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
#include "flat_hash_map/flat_hash_map.hpp"
using namespace sketch;
using namespace bns;

u32 encode_bc(const char *s) {
    const char *e = s + 16;
    u32 ret = cstr_lut[*s++];
    while(s != e) ret <<= 2, ret |= cstr_lut[*s++];
    return ret;
}

int bam_main(int argc, char *argv[]) {
    bam1_t *b = bam_init1();
    const char *path = argv[1];
    samFile *fp = sam_open(path, "r");
    auto hdr = sam_hdr_read(fp);
    int rc;
    const char *tag = "UB"; // TODO: support CB optionally
    ska::flat_hash_map<u32, hll::hll_t> map;
    map.reserve(1 << 16); // why not?
    size_t p = 8;
    // ma
    void *data;
    while((rc = sam_read1(fp, hdr, b)) >= 0) {
        const char *bc(((data = bam_aux_get(b, tag)) != nullptr) ? bam_aux2Z(static_cast<const uint8_t *>(data)): nullptr);
        if(unlikely(bc == nullptr)) throw std::runtime_error(std::string("Missing ") + tag + " tag");
        u32 bcbin = encode_bc(bc);
        auto it = map.find(bcbin);
        if(it == map.end()) it = map.emplace(bcbin, hll::hll_t(p)).first;
    }
    if(rc != EOF) std::fprintf(stderr, "Wrong error code from rc: %d\n", rc);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);
    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    if(std::strcmp(argv[1], "bam") == 0)
        return bam_main(argc - 1, argv + 1);
    RUNTIME_ERROR("Unsupported command");
}
