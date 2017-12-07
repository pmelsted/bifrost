#ifndef KALLISTO_REPHASH_H
#define KALLISTO_REPHASH_H

#include <stdint.h>
#include <cassert>
#include <time.h>

static const unsigned char twin[32] = {
    0, 20, 2, 7, 4, 5, 6, 3,
    8,  9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 1, 21, 22, 23,
    24, 25, 26, 27, 28, 29, 30, 31
};

static const uint64_t hvals[4] = {
    2053695854357871005ULL, 5073395517033431291ULL,
    10060236952204337488ULL, 7783083932390163561ULL
};

class RepHash {

    public:

        RepHash(int _k) : k(_k) { k = k & 63; }

        RepHash() : k(0) {}

        inline void fastleftshiftk(uint64_t& x) { x = (x << k) | (x >> (64-k)); }

        inline void fastrightshiftk(uint64_t& x) { x = (x >> k) | (x << (64-k)); }

        inline void fastleftshift1(uint64_t& x) { x = (x << 1) | (x >> 63); }

        inline void fastrightshift1(uint64_t& x) { x = (x >> 1) | (x << 63); }

        inline uint64_t hash() const { return (h ^ ht); }

        // hash of _s[0:k]
        void init(const char* _s) {

            h = 0;
            ht = 0;

            const unsigned char* s = (const unsigned char*) _s;

            for (size_t i = 0; i < k; ++i) {

                fastleftshift1(h);

                h ^= hvals[charmask(s[i])];

                fastleftshift1(ht);

                ht ^= hvals[twinmask(s[k-1-i])];
            }
        }

        // hash of s[1:]+in where s[0] == out
        inline void update(const unsigned char out, const unsigned char in) {

            uint64_t z(hvals[charmask(out)]);
            uint64_t zt(hvals[twinmask(in)]);

            fastleftshiftk(z);
            fastleftshiftk(zt);

            fastleftshift1(h);

            h ^= z;
            h ^= hvals[charmask(in)];

            ht ^= zt;
            ht ^= hvals[twinmask(out)];

            fastrightshift1(ht);
        }

        // hash of s+fw
        inline void extendFW(const unsigned char fw) {

            // update h, h(a[0.k]) = p(a[k]) + x*h(a[0..k-1])
            fastleftshift1(h);
            h ^= hvals[charmask(fw)];

            // update ht, ht(a[0..k]) = p(~a[k])*x^k + ht(a[0..k-1])
            uint64_t hvalt(hvals[twinmask(fw)]);

            fastleftshiftk(hvalt);
            ht ^= hvalt;

            // extend k by 1
            increaseK();
        }

        // hash of bw+s
        inline void extendBW(const unsigned char bw) {

            // update h, h(a[-1..k-1]) = p(a[-1])*x^k + h(a[0..k-1])
            uint64_t hval(hvals[charmask(bw)]);

            fastleftshiftk(hval);
            h ^= hval;

            // update ht, ht(a[-1..k-1]) = p(~a[-1]) + x * ht(a[0..k-1])
            fastleftshift1(ht);
            ht ^= hvals[twinmask(bw)];

            // extend k by 1
            increaseK();
        }

        inline void increaseK() { k = (k+1) & 63; }

        inline void setK(const int _k) { h = 0; ht = 0; k = _k; }

        inline uint64_t charmask (const unsigned char x) const { return (x & 6) >> 1; }

        inline uint64_t twinmask (const unsigned char x) const { return ((x ^ 4) & 6) >> 1; }

    private:

        size_t k;
        uint64_t h, ht;
};

#endif // KALLISTO_REPHASH_H
