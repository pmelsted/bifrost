#ifndef BIFROST_BLOCKEDBLOOMFILTER_HPP
#define BIFROST_BLOCKEDBLOOMFILTER_HPP

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <atomic>

#include "libdivide.h"
#include "libpopcnt.h"

#define NB_BITS_BLOCK (0x800ULL)
#define MASK_BITS_BLOCK (0x7ffULL)
#define NB_ELEM_BLOCK (32)

/* Short description:
 *  - Extended BloomFilter which hashes into 64-bit blocks
 *    that can be accessed very fast from the CPU cache
 * */

#if defined(__AVX2__)
#include <x86intrin.h>
#endif

class BlockedBloomFilter {

    public:

        typedef std::pair<uint64_t, uint64_t> BBF_Blocks;

        BlockedBloomFilter();
        BlockedBloomFilter(size_t nb_elem, size_t bits_per_elem);
        BlockedBloomFilter(const BlockedBloomFilter& o);
        BlockedBloomFilter(BlockedBloomFilter&& o);

        ~BlockedBloomFilter();

        BlockedBloomFilter& operator=(const BlockedBloomFilter& o);
        BlockedBloomFilter& operator=(BlockedBloomFilter&& o);

        inline BBF_Blocks getBlock(uint64_t min_hash) const{

            uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

            min_hash -= (min_hash / fast_div_) * blocks_;
            min_hash_2 -= (min_hash_2 / fast_div_) * blocks_;

            return std::make_pair(min_hash, min_hash_2);
        }

        int contains(const uint64_t (&kmer_hash)[4], const uint64_t min_hash, bool (&pres)[4], const int limit) const;
        bool contains(const uint64_t kmer_hash, const uint64_t min_hash) const;

        inline bool contains(const uint64_t kmer_hash, const uint64_t min_hash, const BBF_Blocks blockIDs) const {

            return (contains_block(kmer_hash, min_hash, blockIDs) != 0);
        }

        size_t contains_block(const uint64_t kmer_hash, const uint64_t min_hash, const BBF_Blocks blockIDs) const;

        inline bool insert(const uint64_t kmer_hash, const uint64_t min_hash, const bool multi_threaded = false) {

            return (multi_threaded ? insert_par(kmer_hash, min_hash) : insert_unpar(kmer_hash, min_hash));
        }

        bool WriteBloomFilter(FILE *fp) const;
        bool ReadBloomFilter(FILE *fp);

        void clear();

        inline uint64_t getNbBlocks() const { return blocks_; }

    private:

        struct BBF_Block {

            BBF_Block() {

                lck.clear();
                memset(block, 0, NB_ELEM_BLOCK * sizeof(uint64_t));
            }

            inline void lock() { while (lck.test_and_set(std::memory_order_acquire)); }
            inline void unlock() { lck.clear(std::memory_order_release); }

            uint64_t block[NB_ELEM_BLOCK];
            std::atomic_flag lck;
        };

        BBF_Block* table_; //Bit array

        uint64_t blocks_; //Nb blocks

        int k_; //Nb hash functions

        libdivide::divider<uint64_t> fast_div_; // fast division

        #if defined(__AVX2__)
        //__m256i mask_h; // Mask for hash functions that we do not want (>k);

        uint64_t hashes_mask[4];

        static const __m256i mask_and_div; // All MASK_BITS_BLOCK LSB of each 16 bits word
        static const __m256i mask_and_mod; // All 4 LSB of each 16 bits word
        static const __m256i one2shift_lsb; // All 1 LSB of each 32 bits word
        static const __m256i one2shift_msb; // Set the 17th LSB bit of each 32 bits word
        static const __m256i mask_lsb; // All 16 LSB bits of each 32 bits word
        #endif

        void init_table();

        inline double fpp(size_t bits, int k) const {

            return pow(1-exp(-((double)k)/((double)bits)),(double)k);
        }

        bool insert_par(const uint64_t kmer_hash, const uint64_t min_hash);
        bool insert_unpar(const uint64_t kmer_hash, const uint64_t min_hash);
};

#endif // BFG_BLOCKEDBLOOMFILTER_HPP
