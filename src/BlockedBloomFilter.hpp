#ifndef BFG_BLOCKEDBLOOMFILTER_HPP
#define BFG_BLOCKEDBLOOMFILTER_HPP

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "hash.hpp"
#include "libdivide.h"

#include "minHashIterator.hpp"
#include <vector>

extern "C" {
    #include "xxhash.h"
}

static const uint64_t mask[8] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};

/* Short description:
 *  - Extended BloomFilter which hashes into 64-bit blocks
 *    that can be accessed very fast from the CPU cache
 * */
class BlockedBloomFilter {

 private:

    uint64_t *table_; //Bit array
    uint64_t blocks_; //Nb blocks (each block is 512 bits -> 8 * sizeof(uint64_t))
    uint32_t seed_; //Seed for random number generation
    uint64_t size_; //Size of bit array
    size_t k_; //Nb hash functions
    libdivide::divider<uint64_t> fast_div_; // fast division

    public:

        BlockedBloomFilter() : seed_(0), size_(0), table_(NULL), k_(0), blocks_(0), fast_div_() {}
        BlockedBloomFilter(size_t num, size_t bits, uint32_t seed) : seed_(seed), size_(0), table_(NULL), fast_div_() {

            size_ = rndup512(bits*num);
            blocks_ = size_/512;

            init_table();
            init_k(bits);
        }

        ~BlockedBloomFilter() { clear(); }

        inline uint64_t* getBlock(const uint64_t min_hash) const{

            return table_ + 8 * (min_hash - (min_hash / fast_div_) * blocks_);
        }

        template<typename T> bool contains(T x) const { return (search(x) == 0); }

        template<typename T>
        size_t search(T x) const {

            //uint64_t hash; MurmurHash3_x64_64((const void *) &x, sizeof(T), seed_, &hash);
            uint64_t hash = XXH64((const void *) &x, sizeof(T), seed_);

            uint64_t hash0 = hash / fast_div_; // hash0 = hash / blocks;
            uint64_t block = hash - hash0 * blocks_; // blocks = hash % blocks;

            __builtin_prefetch(table_+8*block,0,1);

            for (uint64_t i = 0; i < k_; i++) {
                // 0 <= bit < 512, which bit to set
                uint64_t bit = (hash0) & 0x1ffULL; // equal to hash % 512;
                hash0 = (hash0 * 48271) % (2147483647ULL);

                // we set bit number (id % 8) in byte (table_[block + id/8]) to 1
                uint64_t maskcheck = 1ULL << (bit & 0x3fULL);
                uint64_t loc = 8*block + (bit>>6);

                if ((table_[loc] &  maskcheck) == 0) return k_;
            }

            return 0;
        }

        template<typename T>
        size_t insert(T x) {

            size_t r = 0;

            //uint64_t hash; MurmurHash3_x64_64((const void *) &x, sizeof(T), seed_, &hash);
            uint64_t hash = XXH64((const void *) &x, sizeof(T), seed_);

            uint64_t hash0 = hash / fast_div_; // hash0 = hash / blocks;
            uint64_t block = hash - hash0 * blocks_; // blocks = hash % blocks;

            for(uint64_t i = 0; i < k_; i++) {
                // 0 <= bit < 512, which bit to set
                uint64_t bit = (hash0) & 0x1ffULL; // equal to hash % 512;
                hash0 = (hash0 * 48271) % (2147483647ULL);

                // we set bit number (id % 8) in byte (table_[block + id/8]) to 1
                uint64_t maskcheck = 1ULL << (bit & 0x3fULL);
                uint64_t loc = 8*block + (bit>>6);

                if ((table_[loc] &  maskcheck) == 0) {

                    uint64_t val = __sync_fetch_and_or(table_ + loc, maskcheck);

                    if ((val & maskcheck) == 0) r++;
                }
            }

            return r;
        }

        // With minimizers

        inline bool contains(uint64_t kmer_hash, const uint64_t min_hash) const { return (search(kmer_hash, min_hash) == 0); }

        size_t search(uint64_t kmer_hash, const uint64_t min_hash) const {

            uint64_t block = min_hash - (min_hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+8*block,0,1);

            for (uint64_t i = 0; i < k_; i++) {
                // 0 <= bit < 512, which bit to set
                uint64_t bit = kmer_hash & 0x1ffULL; // equal to hash % 512;
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);

                // we set bit number (id % 8) in byte (table_[block + id/8]) to 1
                uint64_t maskcheck = 1ULL << (bit & 0x3fULL);
                uint64_t loc = 8*block + (bit>>6);

                if ((table_[loc] &  maskcheck) == 0) return k_;
            }

            return 0;
        }

        inline bool contains(uint64_t kmer_hash, const uint64_t* block) const { return (search(kmer_hash, block) == 0); }

        size_t search(uint64_t kmer_hash, const uint64_t* block) const {

            __builtin_prefetch(block,0,1);

            for (uint64_t i = 0; i < k_; i++) {
                // 0 <= bit < 512, which bit to set
                uint64_t bit = kmer_hash & 0x1ffULL; // equal to hash % 512;
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);

                // we set bit number (id % 8) in byte (table_[block + id/8]) to 1
                uint64_t maskcheck = 1ULL << (bit & 0x3fULL);

                if ((block[bit>>6] & maskcheck) == 0) return k_;
            }

            return 0;
        }

        size_t insert(uint64_t kmer_hash, uint64_t min_hash) {

            size_t r = 0;

            uint64_t block = min_hash - (min_hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+8*block,1,1);

            for(uint64_t i = 0; i < k_; i++) {
                // 0 <= bit < 512, which bit to set
                uint64_t bit = kmer_hash & 0x1ffULL; // equal to hash % 512;
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);

                // we set bit number (id % 8) in byte (table_[block + id/8]) to 1
                uint64_t maskcheck = 1ULL << (bit & 0x3fULL);
                uint64_t loc = 8*block + (bit>>6);

                if ((table_[loc] &  maskcheck) == 0) {

                    uint64_t val = __sync_fetch_and_or(table_ + loc, maskcheck);

                    if ((val & maskcheck) == 0) r++;
                }
            }

            return r;
        }

        size_t search_and_insert(uint64_t kmer_hash, uint64_t min_hash) {

            size_t r = 0;

            uint64_t bit, maskcheck, loc, val;
            uint64_t block = min_hash - (min_hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+8*block,1,1);

            for( uint64_t i = 0; i < k_; i++) {

                bit = kmer_hash & 0x1ffULL; // kmer_hash % 512 -> which bit to set in the block
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);

                // we set bit number (id % 8) in byte (table_[block + id/8]) to 1
                maskcheck = 1ULL << (bit & 0x3fULL);
                loc = 8*block + (bit>>6);

                if ((table_[loc] & maskcheck) == 0) {

                    val = __sync_fetch_and_or(table_ + loc, maskcheck);

                    if ((val & maskcheck) == 0) r = k_;
                }
            }

            return r;
        }

        /*
        bool contains(uint64_t kmer_hash, const minHashResult& min) const { return (search(kmer_hash, min) == 0); }

        size_t search(uint64_t kmer_hash, const minHashResult& min) const {

            uint64_t block = min.hash - (min.hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+8*block,0,1);

            for (uint64_t i = 0; i < k_; i++) {
                // 0 <= bit < 512, which bit to set
                uint64_t bit = kmer_hash & 0x1ffULL; // equal to hash % 512;
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);

                // we set bit number (id % 8) in byte (table_[block + id/8]) to 1
                uint64_t maskcheck = 1ULL << (bit & 0x3fULL);
                uint64_t loc = 8*block + (bit>>6);

                if ((table_[loc] &  maskcheck) == 0) return k_;
            }

            return 0;
        }

        size_t insert(uint64_t kmer_hash, const minHashResult& min) {

            size_t r = 0;

            uint64_t block = min.hash - (min.hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+8*block,1,1);

            for(uint64_t i = 0; i < k_; i++) {
                // 0 <= bit < 512, which bit to set
                uint64_t bit = kmer_hash & 0x1ffULL; // equal to hash % 512;
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);

                // we set bit number (id % 8) in byte (table_[block + id/8]) to 1
                uint64_t maskcheck = 1ULL << (bit & 0x3fULL);
                uint64_t loc = 8*block + (bit>>6);

                if ((table_[loc] &  maskcheck) == 0) {

                    uint64_t val = __sync_fetch_and_or(table_ + loc, maskcheck);

                    if ((val & maskcheck) == 0) r++;
                }
            }

            return r;
        }

        size_t search_and_insert(uint64_t kmer_hash, const minHashResult& min) {

            size_t r = 0;

            uint64_t bit, maskcheck, loc, val;
            uint64_t block = min.hash - (min.hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+8*block,1,1);

            for( uint64_t i = 0; i < k_; i++) {

                bit = kmer_hash & 0x1ffULL; // kmer_hash % 512 -> which bit to set in the block
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);

                // we set bit number (id % 8) in byte (table_[block + id/8]) to 1
                maskcheck = 1ULL << (bit & 0x3fULL);
                loc = 8*block + (bit>>6);

                if ((table_[loc] & maskcheck) == 0) {

                    val = __sync_fetch_and_or(table_ + loc, maskcheck);

                    if ((val & maskcheck) == 0) r = k_;
                }
            }

            return r;
        }
        */

        bool WriteBloomFilter(FILE *fp) {

            if (fwrite(&size_,   sizeof(size_),   1, fp) != 1) return false;
            if (fwrite(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
            if (fwrite(&seed_,   sizeof(seed_),   1, fp) != 1) return false;
            if (fwrite(&k_,      sizeof(k_),      1, fp) != 1) return false;

            if (fwrite(table_, sizeof(uint64_t), 8*blocks_, fp) != (8*blocks_)) return false;

            return true;
        }

        bool ReadBloomFilter(FILE *fp) {

            clear();

            if (fread(&size_, sizeof(size_), 1, fp) != 1) return false;
            if (fread(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
            if (fread(&seed_, sizeof(seed_), 1, fp) != 1) return false;
            if (fread(&k_,    sizeof(k_),    1, fp) != 1) return false;

            init_table();

            if (fread(table_, sizeof(uint64_t), 8*blocks_, fp) != (8*blocks_)) return false;

            return true;
        }

        size_t count() const {

            unsigned char *t = (unsigned char *) table_;
            size_t c = 0;

            for (size_t i = 0; i < 64*blocks_; i++) {

                unsigned char u = t[i];

                for (size_t j = 128; j != 0; j = j>>1) {

                    if ((u & j) != 0) c++;
                }
            }

            //std::cout << c << " bits set out of " << size_ << " with k = " << k_ << std::endl;
            if (c != 0) {

                double n = size_*(-log(1.0-((double)c)/size_))/k_;
                //cout << "estimate =" << (size_t)n  << endl;
                return (size_t) n;
            }
            else return 0;
        }

        size_t memory() const {

            size_t m = sizeof(BlockedBloomFilter) + (blocks_ / 64 );
            fprintf(stderr, "BlockedBloomFilter:\t\t%zuMB\n",  m >> 20);
            return m;
        }

        void clear() {

            if (table_ != NULL) free(table_);

            table_ = NULL;
            size_ = 0;
            blocks_ = 0;
        }

    private:

        void init_table() {

            fast_div_ = libdivide::divider<uint64_t>(blocks_);
            posix_memalign((void**)&table_, 64, 8*blocks_*sizeof(table_[0]));
            memset(table_, 0, 8*blocks_*sizeof(table_[0]));
        }

        void init_k(size_t bits) {

            size_t k = (size_t) (bits*log(2));

            if (fpp(bits,k) < fpp(bits,k+1)) k_ = k;
            else k_ = k+1;

            //std::cerr << "k="<<k_<<", fpp="<<fpp(bits,k_) <<  std::endl;
        }

        double fpp(size_t bits, size_t k) const {
            //std::cout << bits<<","<<k<<","<<(-((double)k)/((double)bits)) << std::endl;
            return pow(1-exp(-((double)k)/((double)bits)),(double)k);
        }

        uint64_t rndup512(uint64_t x) const {
            return ((x+511)/512)*512;
        }
};

#endif // BFG_BLOCKEDBLOOMFILTER_HPP
