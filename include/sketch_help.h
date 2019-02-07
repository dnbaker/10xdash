#ifndef SKETCH_HELPERS_H___
#define SKETCH_HELPERS_H___
#include "bonsai/hll/hll.h"
#include "bonsai/hll/ccm.h"
#include "bonsai/hll/mh.h"
#include "bonsai/hll/bbmh.h"
#include "khset/khset.h"

namespace tenx {
using namespace sketch;
using namespace hll;

enum Sketch {
    HLL,
    BLOOM_FILTER,
    RANGE_MINHASH,
    FULL_KHASH_SET,
    COUNTING_RANGE_MINHASH,
    HYPERMINHASH16,
    HYPERMINHASH32,
    TF_IDF_COUNTING_RANGE_MINHASH, // TODO. I'm willing to make two passes through the data for this.
    BB_MINHASH,
    COUNTING_BB_MINHASH, // TODO
};

static const char *sketch_names [] {
    "HLL",
    "BLOOM_FILTER",
    "RANGE_MINHASH",
    "FULL_KHASH_SET",
    "COUNTING_RANGE_MINHASH",
    "HYPERMINHASH16",
    "HYPERMINHASH32",
    "TF_IDF_COUNTING_RANGE_MINHASH",
    "BB_MINHASH",
    "COUNTING_BB_MINHASH",
};
using CRMFinal = mh::FinalCRMinHash<uint64_t, std::greater<uint64_t>, uint32_t>;
template<typename T>
double similarity(const T &a, const T &b) {
    return jaccard_index(a, b);
}
template<>
double similarity<CRMFinal>(const CRMFinal &a, const CRMFinal &b) {
    return a.histogram_intersection(b);
}

struct khset64_t: public kh::khset64_t {
    using super = kh::khset64_t;
    void addh(uint64_t v) {this->insert(v);}
    khset64_t(): kh::khset64_t() {}
    khset64_t(size_t reservesz): kh::khset64_t(reservesz) {}
    khset64_t(const khset64_t &o): super(*reinterpret_cast<const super *>(&o)) {}
    khset64_t(khset64_t &&o): kh::khset64_t() /* memset(0) initialization */ {
        std::swap_ranges((uint8_t *)this, (uint8_t *)this + sizeof(*this), (uint8_t *)std::addressof(o));
    }
    double jaccard_index(const khset64_t &other) const {
        auto p1 = this, p2 = &other;
        if(size() > other.size())
            std::swap(p1, p2);
        size_t olap = 0;
        p1->for_each([&](auto v) {olap += p2->contains(v);});
        return static_cast<double>(olap) / (p1->size() + p2->size() - olap);
    }
};
using CBBMinHashType = mh::CountingBBitMinHasher<uint64_t, uint32_t>; // Is counting to 65536 enough for a transcriptome? Maybe we can use 16...
template<typename SketchType>
struct FinalSketch {
    using final_type = SketchType;
};
#define FINAL_OVERLOAD(type) \
template<> struct FinalSketch<type> { \
    using final_type = typename type::final_type;}
FINAL_OVERLOAD(mh::CountingRangeMinHash<uint64_t>);
FINAL_OVERLOAD(mh::RangeMinHash<uint64_t>);
FINAL_OVERLOAD(mh::BBitMinHasher<uint64_t>);
FINAL_OVERLOAD(CBBMinHashType);

template<typename T>struct SketchFileSuffix {static constexpr const char *suffix = ".sketch";};
#define SSS(type, suf) template<> struct SketchFileSuffix<type> {static constexpr const char *suffix = suf;}
SSS(mh::CountingRangeMinHash<uint64_t>, ".crmh");
SSS(mh::RangeMinHash<uint64_t>, ".rmh");
SSS(khset64_t, ".khs");
SSS(bf::bf_t, ".bf");
SSS(mh::BBitMinHasher<uint64_t>, ".bmh");
SSS(CBBMinHashType, ".cbmh");
SSS(mh::HyperMinHash<uint64_t>, ".hmh");
SSS(hll::hll_t, ".hll");

}
#endif /* SKETCH_HELPERS_H___ */
