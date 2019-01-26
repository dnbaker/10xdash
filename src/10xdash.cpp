#include "bonsai/bonsai/include/util.h"
#include "bonsai/bonsai/include/database.h"
#include "bonsai/bonsai/include/bitmap.h"
#include "bonsai/bonsai/include/setcmp.h"
#include "bonsai/hll/hll.h"
#include "bonsai/hll/ccm.h"
#include "bonsai/pdqsort/pdqsort.h"
#include "distmat/distmat.h"
#include <sstream>
#include "htslib/htslib/sam.h"
#include "klib/ketopt.h"
#include "flat_hash_map/flat_hash_map.hpp"
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
};

int bam_main(int argc, char *argv[]) {
    if(argc == 1) return bam_usage();
    bam1_t *b = bam_init1();
    const char *path = argv[1];
    samFile *fp = sam_open(path, "r");
    auto hdr = sam_hdr_read(fp);
    ketopt_t opt = KETOPT_INIT;
    int rc, nthreads = 1;
    int p = 8, full_p = 0, k = 31;
    bool skip_full = false;
    size_t map_reserve_size = 1 << 16;
    BCType tagtype = CB;
    while((rc = ketopt(&opt, argc, argv, 0, "k:R:p:brh?", nullptr)) >= 0) {
        switch(rc) {
           case 'p': nthreads = std::atoi(optarg); assert(nthreads > 0); break;
           case 'b': tagtype = UB; break;
           case 'r': tagtype = UR; break;
           case 'S': p = std::atoi(optarg); break;
           case 'k': k = std::atoi(optarg); break;
           case 'f': full_p = std::atoi(optarg); break;
           case 'P': skip_full = true; break;
           case 'R': map_reserve_size = std::strtoull(optarg, nullptr, 10); break;
           case 'h': return bam_usage();
        }
    }
    if(!full_p) full_p = p;
    const char *const tag = tags[tagtype];
    ska::flat_hash_map<u32, hll::hll_t> map;
    map.reserve(map_reserve_size); // why not?
    // ma
    void *data;
    hll::hll_t *full_set = skip_full ? nullptr: new hll::hll_t(full_p);
    sketch::common::WangHash hasher;
    const uint64_t kmer_mask = UINT64_C(-1) >> (64 - (k * 2));
    while((rc = sam_read1(fp, hdr, b)) >= 0) {
        if(b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY)) continue;
        uint8_t *data = bam_aux_get(b, tag);
        if(unlikely(data == nullptr)) throw std::runtime_error(std::string("Missing ") + tag + " tag");
        data = reinterpret_cast<uint8_t *>(bam_aux2Z(data));
        u32 bcbin = encode_bc(reinterpret_cast<char *>(data));
        auto it = map.find(bcbin);
        if(it == map.end()) it = map.emplace(bcbin, hll::hll_t(p)).first;
        data = bam_get_seq(b);
        int i = 0, len = b->core.l_qseq, nfilled;
        u64 kmer = lut[bam_seqi(data, 0)];
        if(kmer == BF)
            i += k, kmer = 0, nfilled = 0;
        else nfilled = 1;
        while(i < len) {
            kmer <<= 2;
            if((kmer |= lut[bam_seqi(data, i)]) == BF)
                kmer = 0, i += k;
            else {
                kmer &= kmer_mask;
                if(++nfilled == k) {
                    uint64_t hv = hasher(kmer);
                    it->second.add(hv);
                    if(full_set) full_set->add(hv);
                    --nfilled;
                }
            }
        }
    }
    const size_t map_size = map.size();
    hll::hll_t *ptrs = static_cast<hll::hll_t *>(std::malloc(sizeof(hll::hll_t) * map_size));
    if(ptrs == nullptr) throw std::bad_alloc();
    auto pp = ptrs;
    for(auto &&pair: map) {
        new(pp++) hll::hll_t(std::move(pair.second));
    }
    // Core distance code
    std::for_each(ptrs, ptrs + map_size, [](auto &h) {h.~hllbase_t();});
    std::free(ptrs);
    if(rc != EOF) std::fprintf(stderr, "Wrong error code from rc: %d\n", rc);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);
    delete full_set;
    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    ex = argv[0];
    return bam_main(argc, argv);
    if(std::strcmp(argv[1], "bam") == 0)
        return bam_main(argc - 1, argv + 1);
    RUNTIME_ERROR("Unsupported command");
}
