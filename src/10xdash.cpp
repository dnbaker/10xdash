#include "bonsai/bonsai/include/bitmap.h"
#include "bonsai/bonsai/include/setcmp.h"
#include "distmat/distmat.h"
#include "htslib/htslib/sam.h"
#include "klib/ketopt.h"
#include "flat_hash_map/flat_hash_map.hpp"
#include <new>
#include <list>
#include <omp.h>
#include "include/sketch_help.h"
#include "include/util.h"
#if __cplusplus >= 201703L
#  include <memory>
#  if __has_include(<execution>)
#    define USE_PAR_EX
#    include <execution>
#  endif
#endif
using namespace bns;
using namespace tenx;

static const char *ex = nullptr;
static uint32_t bbnbits = 16;
static const char *suf = ".unset_suffix";

int bam_usage() {
    std::fprintf(stderr, "Usage:\n%s <flags> in.{sb}am\n%s currently does not assume any sorted ordering or alignment and compares by HLLs by default.\n"
                         "-B\tUse Bloom Filters for sketches. These become more accurate with very large sketches. Default: HLL\n"
                         "-K\tUse full hash sets instead of sketches. These are exact, but expensive to compare and require a relatively larger amount of memory. Default: HLL\n"
                         "-M\tUse bottom-k minhash. Default: HLL\n"
                         "-C\tUse counting bottom-k minhash histogram intersection. Default: HLL\n"
                         "-S\tSketch size in bytes, log2. Default: 8 (256 bytes per sketch)\n"
                         "-k\tKmer length. Default: 31\n"
                         "-f\tFail all reads without all bits in argument set. This can be specified multiple times. Default: 0\n"
                         "-F\tFail all reads with any bits in argument set. This can be specified multiple times. Default: 0\n"
                         "-o\tOutput path. Default: in.{sb}am.distmat\n"
                         "-p\tNumber of threads to use. Default: 1\n"
                         "-b\tCheck 'UB' tag. Default: 'CB'\n"
                         "-r\tCheck 'UR' tag. Default: 'CB'\n"
                         "-s\tIncrement in log2 size of full experiment sketch in bytes. Default: 0\n"
                         "-R\tMap reserve size. Pre-allocate this much space in the hash table. Default: 1 << 16\n"
                         "-z\tOutput zlib compression level. Set to 0 for uncompressed. Default: 0\n"
                         "-w\tWrite sketches to disk. This will be done in one file per barcode\n"
                         "-D\tDo not perform distance calculations. (This should only be done is -w is specified.)\n"
                          "-h/-?\tEmit this usage menu.\n"
                 , ex, ex);
    return EXIT_FAILURE;
}


struct CLIArgs {
    int nthreads = 1;
    int sketch_size_l2 = 8;
    int full_sketch_size_l2_diff = 0; // Defaults to sketch_size_l2
    unsigned fail_flags = 0;
    unsigned required_flags = 0;
    int k = 31;
    int compression_level = 0;
    bool skip_full = false;
    bool write_sketches = false;
    bool skip_distance = false;
    bool write_human_readable = true;
    std::string write_binary_bc_sketch_pairs;
    size_t map_reserve_size = 1 << 16;
    samFile *fp = nullptr;
    bam_hdr_t *hdr = nullptr;
    BCType tag = CB;
    Sketch sketch_type = HLL;
    size_t parallel_chunk_size = 1 << 10;
    const char *omatpath = nullptr;
    static size_t bytesl2_to_arg(int nblog2, Sketch sketch) {
        switch(sketch) {
            case HLL: return nblog2;
            case BLOOM_FILTER: return nblog2 + 3; // 8 bits per byte
            case HYPERMINHASH16: return nblog2 - 1; // Assuming 16-bit HMH
            case HYPERMINHASH32: return nblog2 - 2; // Assuming 16-bit HMH
            case RANGE_MINHASH: return size_t(1) << (nblog2 - 3); // 8 bytes per minimizer
            case FULL_KHASH_SET: return size_t(1) << nblog2; // No real reason
            case COUNTING_RANGE_MINHASH: return size_t(1) << (nblog2) / (sizeof(uint64_t) + sizeof(uint32_t));
            case BB_MINHASH: bbnbits; // since b comes before p for this constructor. TODO: consider allowing different values.
            default: {
                char buf[128];
                std::sprintf(buf, "Sketch %s not yet supported.\n", (sketch >= (sizeof(sketch_names) / sizeof(char *)) ? "Not such sketch": sketch_names[sketch]));
                RUNTIME_ERROR(buf);
                return -1337;
            }
        }

    }
    int next_rec(bam1_t *b) {return sam_read1(fp, hdr, b);}
};


template<typename SketchType, typename... Args>
int core(CLIArgs &args, dm::DistanceMatrix<float> *distmat, uint32_t **bcs, Args &&... sketchargs) {
    const char *const tag = tags[args.tag];
    suf = SketchFileSuffix<SketchType>::suffix;
    bam1_t *b = bam_init1();
    ska::flat_hash_map<u32, SketchType> map;
    map.reserve(args.map_reserve_size); // why not?
    sketch::common::WangHash hasher;
    const uint64_t kmer_mask = UINT64_C(-1) >> (64 - (args.k * 2));
    SketchType *full_set =
        args.skip_full ? nullptr
                       : new SketchType(args.bytesl2_to_arg(
                                            args.sketch_size_l2 +
                                                args.full_sketch_size_l2_diff, args.sketch_type),
                                          std::forward<Args>(sketchargs)...);
    int rc;
    while((rc = args.next_rec(b)) >= 0) {
        if(b->core.flag & (args.fail_flags) || (b->core.flag & args.required_flags) != args.required_flags) continue;
        uint8_t *data = bam_aux_get(b, tag);
        if(unlikely(data == nullptr)) RUNTIME_ERROR(std::string("Missing ") + tag + " tag");
        data = reinterpret_cast<uint8_t *>(bam_aux2Z(data));
        u32 bcbin = encode_bc(reinterpret_cast<char *>(data));
        auto it = map.find(bcbin);
        if(it == map.end())
            it = map.emplace(bcbin,
                             SketchType(args.bytesl2_to_arg(
                                            args.sketch_size_l2, args.sketch_type),
                                        std::forward<Args>(sketchargs)...)).first;
        data = bam_get_seq(b);
        int i = 0, len = b->core.l_qseq, nfilled;
        u64 kmer = lut4b[bam_seqi(data, 0)];
        if(kmer == BF)
            kmer = 0, nfilled = 0;
        else nfilled = 1;
        while(i < len) {
            kmer <<= 2;
            if((kmer |= lut4b[bam_seqi(data, i)]) == BF) {
                kmer = nfilled = 0;
            } else if(++nfilled == args.k) {
                kmer &= kmer_mask;
                it->second.addh(kmer);
                if(full_set) full_set->addh(kmer);
                --nfilled;
            }
            ++i;
        }
    }
    bam_destroy1(b);
    std::future<void> full_writer;
    if(full_set) {
        std::string opath = args.fp->fn;
        opath += SketchFileSuffix<SketchType>::suffix;
        full_writer = std::async(std::launch::async, [&]() {
            full_set->write(opath.data());
            delete full_set;
        });
    }
    std::list<std::future<void>> q;
    if(args.write_sketches) {
        for(const auto &pair: map) {
            q.emplace_back(std::async(std::launch::async, [&]() {
                std::string opath = std::to_string(pair.first) + SketchFileSuffix<SketchType>::suffix;
                pair.second.write(opath.data());
            }));
            while(q.size() >= args.nthreads) {
                for(auto it = q.begin(); it != q.end();) {
                    if(it->wait_for(std::chrono::seconds(0)) == std::future_status::ready)
                        q.erase(it++);
                    else ++it;
                }
            }
        }
    }
    if(args.write_binary_bc_sketch_pairs.size()) {
        auto names = args.write_binary_bc_sketch_pairs + SketchFileSuffix<SketchType>::suffix;
        gzFile fp = gzopen(names.data(), "wb");
        uint64_t nsketches = map.size();
        gzwrite(fp, &nsketches, sizeof(nsketches));
        for(const auto &pair: map) {
            gzwrite(fp, &pair.first, sizeof(pair.first)); // Write 32-bit int.
            pair.second.write(fp);
        }
        gzclose(fp);
    }
    if(args.skip_distance)
        return EXIT_SUCCESS;
    const size_t map_size = map.size();
    using FinalType = typename FinalSketch<SketchType>::final_type;
    FinalType *ptrs;
    if(ptrs == nullptr) throw std::bad_alloc();
    uint32_t *barcodes = static_cast<uint32_t *>(std::malloc(sizeof(uint32_t) * map_size)), *bcp = barcodes;
    if(barcodes == nullptr) throw std::bad_alloc();
    *bcs = barcodes;
    auto pp = ptrs;
    for(auto &&pair: map) {
        *bcp++ = pair.first;
        new(pp++) FinalType(std::move(pair.second));
    }
    {decltype(map) tmap(std::move(map));} // Free map now that it's not needed
    distmat->resize(map_size);
    omp_set_num_threads(args.nthreads);
    if(map_size < args.parallel_chunk_size) {
        for(size_t i = 0; i < map_size; ++i) {
            auto &cmpsketch = ptrs[i];
            auto span = distmat->row_span(i);
            #pragma omp parallel for
            for(size_t j = i + 1; j < map_size; ++j) {
                assert(j >= i - 1 && (j - i - 1 <= span.second));
                span.first[j - i - 1] = similarity(cmpsketch, ptrs[j]);
            }
        }
    }
    // Core distance code
#if __cplusplus >= 201703L
    std::destroy_n(
#  ifdef USE_PAR_EX
        std::execution::par_unseq,
#  endif
        ptrs, map_size);
#else
    std::for_each(ptrs, pp, [](auto &sketch) {sketch.~FinalType();});
#endif
    std::free(ptrs);
    if(rc != EOF) std::fprintf(stderr, "Warning: Wrong error code from rc: %d\n", rc);
    return rc;
}

int bam_main(int argc, char *argv[]) {
    if(argc == 1) return bam_usage();
    CLIArgs args;
    ketopt_t opt = KETOPT_INIT;
    // TODO: Add count{,min}-sketch filtering for errors in sequencing
    // TODO: Perform richer comparisons based on regions within the genome
    int rc;
    while((rc = ketopt(&opt, argc, argv, 0, "W:s:o:k:R:p:S:z:f:F:N:DPKCdwdBbrh?", nullptr)) >= 0) {
        switch(rc) {
            case '3': args.sketch_type = HYPERMINHASH32;
            case 'B': args.sketch_type = BLOOM_FILTER; break;
            case 'C': args.sketch_type = COUNTING_RANGE_MINHASH; break;
            case 'D': args.skip_distance = true; break;
            case 'F': args.required_flags |= std::atoi(optarg); break;
            case 'H': args.sketch_type = HYPERMINHASH16;
            case 'K': args.sketch_type = FULL_KHASH_SET; break;
            case 'N': bbnbits = std::atoi(optarg); break;
            case 'P': args.skip_full = true; break;
            case 'R': args.map_reserve_size = std::strtoull(opt.arg, nullptr, 10); break;
            case 'W': args.write_binary_bc_sketch_pairs = opt.arg; break;
            case 'S': args.sketch_size_l2 = std::atoi(opt.arg); break;
            case 'b': args.tag = UB; break;
            case 'd': args.write_human_readable = false;
            case 'f': args.fail_flags |= std::atoi(optarg); break;
            case 'k': args.k = std::atoi(opt.arg); break;
            case 'o': args.omatpath = optarg; break;
            case 'p': args.nthreads = std::atoi(opt.arg); assert(args.nthreads > 0); break;
            case 'r': args.tag = UR; break;
            case 's': args.full_sketch_size_l2_diff = std::atoi(opt.arg); break;
            case 'w': args.write_sketches = true; break;
            case 'z': args.compression_level = std::atoi(optarg) % 10; break; break;
            case 'h': case '?': return bam_usage();
        }
    }
    if(args.sketch_type == HYPERMINHASH16 || args.sketch_type == HYPERMINHASH32) {
        std::fprintf(stderr, "Warning: something seems wrong with this HyperMinHash draft implementation. It is experimental and not as accurate as other, more vetted methods.\n");
    }
    samFile *fp;
    bam_hdr_t *hdr;
    if(!(fp = sam_open(argv[opt.ind], "r"))) RUNTIME_ERROR(std::string("Could not open sam at ") + argv[opt.ind]);
    if(!(hdr = sam_hdr_read(fp)))            RUNTIME_ERROR(std::string("Could not parse header from file at ") + argv[opt.ind]);
    dm::DistanceMatrix<float> distances;
    uint32_t *bcs = nullptr;
    switch(args.sketch_type) {
        case HLL: core<hll::hll_t>(args, &distances, &bcs); break;
        case BLOOM_FILTER: core<bf::bf_t>(args, &distances, &bcs); break;
        case RANGE_MINHASH: core<mh::RangeMinHash<uint64_t>>(args, &distances, &bcs); break;
        case FULL_KHASH_SET: core<khset64_t>(args, &distances, &bcs); break;
        case COUNTING_RANGE_MINHASH: core<mh::CountingRangeMinHash<uint64_t>>(args, &distances, &bcs); break;
        case HYPERMINHASH16: core<mh::HyperMinHash<uint64_t>>(args, &distances, &bcs, 10); break;
        case HYPERMINHASH32: core<mh::HyperMinHash<uint64_t>>(args, &distances, &bcs, 26); break;
        case BB_MINHASH:    core<mh::BBitMinHasher<uint64_t>>(args, &distances, &bcs, args.sketch_size_l2 - std::ceil(std::log2(bbnbits / 8))); break;
        default: {
            char buf[128];
            std::sprintf(buf, "Sketch %s not yet supported.\n", args.sketch_type >= sizeof(sketch_names) / sizeof(char *) ? "Not such sketch": sketch_names[args.sketch_type]);
            RUNTIME_ERROR(buf);
        }
    }
    if(distances.size()) {
        std::string opath = args.omatpath ? args.omatpath: (std::string(argv[opt.ind]) + ".distmat").data();
        gzFile ofp = gzopen(opath.data(), (std::string("wb") + std::to_string(args.compression_level)).data());
        if(!ofp) RUNTIME_ERROR(std::string("Could not open file for writing at ") + opath.data());
        if(args.write_human_readable) {
            gzputs(ofp, "#Labels");
            std::for_each(bcs, bcs + distances.size(), [ofp](auto bc) {gzprintf(ofp, "\t%u", bc);});
            gzputc(ofp, '\n');
            distances.printf(ofp, true); // Always emit scientific fmt for now
        } else {
            std::FILE *cfp = std::fopen((opath + ".labels").data(), "w");
            if(!cfp) RUNTIME_ERROR(std::string("Could not open file at ") + opath + ".labels");
            std::fprintf(cfp, "#Barcodes\n");
            std::for_each(bcs, bcs + distances.size(), [cfp](auto bc) {std::fprintf(cfp, "%u\n", bc);});
            distances.write(ofp);
            std::fclose(cfp);
        }
        gzclose(ofp);
    }
    std::free(bcs);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    ex = argv[0];
    return bam_main(argc, argv);
}
