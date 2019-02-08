#ifndef TENX_DASH_UTIL_H__
#define TENX_DASH_UTIL_H__
#include <cstdint>
#include "bonsai/bonsai/include/util.h"

namespace tenx {
using namespace bns;
static inline u32 encode_bc(const char *s) {
    // Be willing to consider parsing the last integer describing the GEM group
    const char *e = s + 16;
    u32 ret = cstr_lut[*s++];
    while(s != e) ret <<= 2, ret |= cstr_lut[*s++];
    return ret;
}

enum BCType {
    CB = 0,
    UB = 1,
    BOTH = 2
};

static const int8_t lut4b [] {
    -1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1
    // For 4-bit nucleotides
    // We ignore anything that isn't A, C, G, or T
    // the |= -1 converts the kmer to UINT64_C(-1), which we reserve as an error code.
};

} // tenx

#endif
