#include "bonsai/bonsai/include/util.h"
#include "bonsai/bonsai/include/bitmap.h"
#include "bonsai/bonsai/include/setcmp.h"
#include "bonsai/hll/hll.h"
#include "bonsai/hll/ccm.h"
#include "distmat/distmat.h"
#include "htslib/htslib/sam.h"
#include "klib/ketopt.h"
#include "flat_hash_map/flat_hash_map.hpp"
#include <new>
#include <omp.h>
#if __cplusplus >= 201703L
#  include <memory>
#  if __has_include(<execution>)
#    define USE_PAR_EX
#    include <execution>
#  endif
#endif
using namespace sketch;
using namespace bns;

static const char *ex = nullptr;

u32 encode_bc(const char *s) {
    // Be willing to consider parsing the last integer describing the GEM group
    const char *e = s + 16;
    u32 ret = cstr_lut[*s++];
    while(s != e) ret <<= 2, ret |= cstr_lut[*s++];
    return ret;
}

enum BCType {
    CB = 0,
    UB = 1,
    UR = 2
};

static const char *tags[] {
    "CB", "UB", "UR"
};

int bam_usage() {
    std::fprintf(stderr, "TODO: write usage.\nUsage: %s <opts>\n", ex);
    return EXIT_FAILURE;
}

static const int8_t lut [] {
    -1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1
    // For 4-bit nucleotides
    // We ignore anything that isn't A, C, G, or T
    // the |= -1 converts the kmer to UINT64_C(-1), which we reserve as an error code.
};

enum Sketch {
    HLL = 0,
    BLOOM_FILTER = 1,
    RANGE_MIN_HASH = 2,
    HYPER_MIN_HASH = 3,
    FULL_KHASH_SET = 4 // Not yet supported
};

struct CLIArgs {
    int nthreads = 1;
    int sketch_size_l2 = 8;
    int full_sketch_size_l2_diff = 0; // Defaults to sketch_size_l2
    int k = 31;
    int compression_level = 0;
    bool skip_full = false;
    bool write_sketches_only = false;
    bool write_human_readable = true;
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
            case HYPER_MIN_HASH: return nblog2 - 1; // Assuming 16-bit HMH
            case RANGE_MIN_HASH: return size_t(1) << nblog2;
            default: RUNTIME_ERROR("Impocerous!"); return -1337;
        }
        
    }
    int next_rec(bam1_t *b) {
        return sam_read1(fp, hdr, b);
    }
};

template<typename SketchType, typename... Args>
int core(CLIArgs &args, dm::DistanceMatrix<float> *distmat, uint32_t **bcs, Args &&... sketchargs) {
    const char *const tag = tags[args.tag];
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
        if(b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY)) continue;
        uint8_t *data = bam_aux_get(b, tag);
        if(unlikely(data == nullptr)) RUNTIME_ERROR(std::string("Missing ") + tag + " tag");
        data = reinterpret_cast<uint8_t *>(bam_aux2Z(data));
        u32 bcbin = encode_bc(reinterpret_cast<char *>(data));
        auto it = map.find(bcbin);
        if(it == map.end()) it = map.emplace(bcbin, SketchType(args.bytesl2_to_arg(args.sketch_size_l2, args.sketch_type), std::forward<Args>(sketchargs)...)).first;
        data = bam_get_seq(b);
        int i = 0, len = b->core.l_qseq, nfilled;
        u64 kmer = lut[bam_seqi(data, 0)];
        if(kmer == BF)
            kmer = 0, nfilled = 0; // Skip 'k', as 
        else nfilled = 1;
        while(i < len) {
            kmer <<= 2;
            if((kmer |= lut[bam_seqi(data, i)]) == BF)
                kmer = nfilled = 0;
            else {
                kmer &= kmer_mask;
                if(++nfilled == args.k) {
                    it->second.addh(kmer);
                    if(full_set) full_set->addh(kmer);
                    --nfilled;
                }
            }
            ++i;
        }
    }
    bam_destroy1(b);
    std::string opath;
    if(full_set) {
        opath = args.fp->fn;
        opath += ".hll";
        full_set->write(opath.data());
    }
    if(args.write_sketches_only) {
        for(const auto &pair: map) {
            opath = pair.first;
            opath += ".hll";
            pair.second.write(opath.data());
        }
        return EXIT_SUCCESS;
    }
    const size_t map_size = map.size();
    SketchType *ptrs = static_cast<SketchType *>(std::malloc(sizeof(SketchType) * map_size));
    if(ptrs == nullptr) throw std::bad_alloc();
    uint32_t *barcodes = static_cast<uint32_t *>(std::malloc(sizeof(uint32_t) * map_size)), *bcp = barcodes;
    *bcs = barcodes;
    auto pp = ptrs;
    for(auto &&pair: map) {
        *bcp++ = pair.first;
        new(pp++) SketchType(std::move(pair.second));
    }
    distmat->resize(map_size);
    omp_set_num_threads(args.nthreads);
    if(map_size < args.parallel_chunk_size) {
        for(size_t i = 0; i < map_size; ++i) {
            auto &cmpsketch = ptrs[i];
            auto span = distmat->row_span(i);
            #pragma omp parallel for
            for(size_t j = i + 1; j < map_size; ++j) {
                assert(j >= i - 1 && (j - i - 1 <= span.second));
                span.first[j - i - 1] = jaccard_index(cmpsketch, ptrs[j]);
            }
        }
    }
    // Core distance code
#if __cplusplus >= 201703L
    std::destroy_n(
#ifdef USE_PAR_EX
        std::execution::par,
#endif
        ptrs, map_size);
#else
    std::for_each(ptrs, pp, [](auto &sketch) {sketch.~SketchType();});
#endif
    std::free(ptrs);
    if(rc != EOF) std::fprintf(stderr, "Warning: Wrong error code from rc: %d\n", rc);
    delete full_set;
    return rc;
}


int bam_main(int argc, char *argv[]) {
    if(argc == 1) return bam_usage();
    CLIArgs args;
    const char *path = argv[1];
    ketopt_t opt = KETOPT_INIT;
    int rc;
    while((rc = ketopt(&opt, argc, argv, 0, "o:k:R:p:dBbrh?", nullptr)) >= 0) {
        switch(rc) {
            case 'B': args.sketch_type = BLOOM_FILTER; break;
            case 'p': args.nthreads = std::atoi(opt.arg); assert(args.nthreads > 0); break;
            case 'b': args.tag = UB; break;
            case 'd': args.write_human_readable = false;
            case 'r': args.tag = UR; break;
            case 'S': args.sketch_size_l2 = std::atoi(opt.arg); break;
            case 'k': args.k = std::atoi(opt.arg); break;
            case 'f': args.full_sketch_size_l2_diff = std::atoi(opt.arg); break;
            case 'P': args.skip_full = true; break;
            case 'w': args.write_sketches_only = true; break;
            case 'R': args.map_reserve_size = std::strtoull(opt.arg, nullptr, 10); break;
            case 'o': args.omatpath = optarg; break;
            case 'z': args.compression_level = 6; break;
            case 'h': return bam_usage();
        }
    }
    samFile *fp = sam_open(path, "r");
    if(!fp) RUNTIME_ERROR(std::string("Could not open sam at ") + path);
    auto hdr = sam_hdr_read(fp);
    if(!hdr) RUNTIME_ERROR(std::string("Could not parse header from file at ") + path);
    dm::DistanceMatrix<float> distances;
    uint32_t *bcs = nullptr;
    switch(args.sketch_type) {
        case HLL: core<hll::hll_t>(args, &distances, &bcs); break;
        case BLOOM_FILTER: core<bf::bf_t>(args, &distances, &bcs); break;
        default: throw std::runtime_error(std::string("NotImplemented sketch type ") + std::to_string(args.sketch_type));
    }
    if(distances.size()) {
        std::string opath = args.omatpath ? args.omatpath: (std::string(path) + ".distmat").data();
        gzFile ofp = gzopen(opath.data(), (std::string("wb") + std::to_string(args.compression_level)).data());
        if(!ofp) RUNTIME_ERROR(std::string("Could not open file for writing at ") + opath.data());
        if(args.write_human_readable) {
            gzputs(ofp, "#Labels");
            std::for_each(bcs, bcs + distances.size(), [ofp](auto bc) {gzprintf(ofp, "\t%u", bc);});
            gzputc(ofp, '\n');
            distances.printf(ofp, true); // Always emit scientific fmt for now
        }
        else {
            std::FILE *cfp = std::fopen((opath + ".labels").data(), "w");
            if(!cfp) RUNTIME_ERROR(std::string("Could not open file at ") + opath + ".labels");
            std::fprintf(cfp, "#Barcodes\n");
            std::for_each(bcs, bcs + distances.size(), [cfp](auto bc) {std::fprintf(cfp, "%u\n", bc);});
            distances.write(ofp);
            std::fclose(cfp);
        }
        gzclose(ofp);
    }
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    ex = argv[0];
    return bam_main(argc, argv);
}
