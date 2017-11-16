#ifndef BFG_BLOCKEDBLOOMFILTER_HPP
#define BFG_BLOCKEDBLOOMFILTER_HPP

#include <cmath>
#include <iostream>
#include <cstdlib>

//#include "hash.hpp"
#include "libdivide.h"

#include <vector>


// ------ TEST ------
#include <algorithm>
#include <fstream>
#include <random>

#include "KmerHashTable.h"
#include "Kmer.hpp"
#include "RepHash.hpp"

#include "libpopcnt.h"

#define NB_BITS_BLOCK (0x800ULL)
#define MASK_BITS_BLOCK (0x7ffULL)
#define NB_ELEM_BLOCK (32)

/* Short description:
 *  - Extended BloomFilter which hashes into 64-bit blocks
 *    that can be accessed very fast from the CPU cache
 * */
class BlockedBloomFilter {

    private:

        uint64_t* table_; //Bit array

        uint64_t size_table_; //Size of bit array (in bits)
        uint64_t blocks_; //Nb blocks
        int k_; //Nb hash functions

        libdivide::divider<uint64_t> fast_div_; // fast division

    public:

        BlockedBloomFilter() : table_(NULL), size_table_(0), blocks_(0), k_(0), fast_div_() {}

        BlockedBloomFilter(size_t nb_elem, size_t bits_per_elem) : table_(NULL), size_table_(0), blocks_(0), k_(0), fast_div_() {

            size_table_ = ((bits_per_elem * nb_elem + MASK_BITS_BLOCK) / NB_BITS_BLOCK) * NB_BITS_BLOCK;
            blocks_ = size_table_ / NB_BITS_BLOCK;

            init_table();

            k_ = (int) (bits_per_elem * log(2));
            if (fpp(bits_per_elem, k_) >= fpp(bits_per_elem, k_+1)) k_++;
        }

        ~BlockedBloomFilter() {

            clear();
        }

        inline std::pair<uint64_t*,uint64_t*> getBlock(uint64_t min_hash) const{

            uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

            min_hash -= (min_hash / fast_div_) * blocks_;
            min_hash_2 -= (min_hash_2 / fast_div_) * blocks_;

            return std::make_pair(table_ + NB_ELEM_BLOCK * min_hash, table_ + NB_ELEM_BLOCK * min_hash_2);
        }

        bool contains(uint64_t kmer_hash, const uint64_t min_hash) const {

            int i = 0;

            const int k = k_;

            uint64_t kmer_hash_2 = kmer_hash;

            uint64_t* table = table_ + ((min_hash - (min_hash / fast_div_) * blocks_) * NB_ELEM_BLOCK);

            __builtin_prefetch(table, 0, 1);

            for (; i < k; ++i) {

                if ((table[(kmer_hash & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash & 0x3fULL))) == 0) break;
                kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
            }

            if (i != k){

                const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

                table = table_ + ((min_hash_2 - (min_hash_2 / fast_div_) * blocks_) * NB_ELEM_BLOCK);

                __builtin_prefetch(table, 0, 1);

                for (i = 0; i < k; ++i) {

                    if ((table[(kmer_hash_2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_2 & 0x3fULL))) == 0) break;
                    kmer_hash_2 = (kmer_hash_2 * 49157) % (1610612741ULL);
                }
            }

            return i == k;
        }

        inline bool contains(uint64_t kmer_hash, const std::pair<uint64_t*, uint64_t*> block_ptr) const {

            return (contains_block(kmer_hash, block_ptr) != 0);
        }

        size_t contains_block(uint64_t kmer_hash, const std::pair<const uint64_t* const, const uint64_t* const> block_ptr) const {

            uint64_t kmer_hash_2 = kmer_hash;

            int i = 0;

            const int k = k_;

            __builtin_prefetch(block_ptr.first, 0, 1);

            for (; i != k; ++i) {

                if ((block_ptr.first[(kmer_hash & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash & 0x3fULL))) == 0) break;
                kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
            }

            if (i != k){

                __builtin_prefetch(block_ptr.second, 0, 1);

                for (i = 0; i != k; ++i) {

                    if ((block_ptr.second[(kmer_hash_2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_2 & 0x3fULL))) == 0) break;
                    kmer_hash_2 = (kmer_hash_2 * 49157) % (1610612741ULL);
                }

                return (i == k ? 2 : 0);
            }

            return 1;
        }

        bool search_and_insert(uint64_t kmer_hash, const uint64_t min_hash, const bool multi_threaded = false) {

            int i = 0, j = 0;

            const int k = k_;

            uint64_t kmer_hash_2 = kmer_hash;

            uint64_t* table = table_ + ((min_hash - (min_hash / fast_div_) * blocks_) * NB_ELEM_BLOCK);

            __builtin_prefetch(table, 0, 1);

            for (; i != k; ++i) {

                if ((table[(kmer_hash & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash & 0x3fULL))) == 0) break;
                kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
            }

            if (i != k){

                const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

                uint64_t* table2 = table_ + ((min_hash_2 - (min_hash_2 / fast_div_) * blocks_) * NB_ELEM_BLOCK);

                __builtin_prefetch(table2, 0, 1);

                for (; j != k; ++j) {

                    if ((table2[(kmer_hash_2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_2 & 0x3fULL))) == 0) break;
                    kmer_hash_2 = (kmer_hash_2 * 49157) % (1610612741ULL);
                }

                if (j != k){

                    if (!multi_threaded){

                        if (popcnt(table2, NB_ELEM_BLOCK * sizeof(uint64_t)) < popcnt(table, NB_ELEM_BLOCK * sizeof(uint64_t))){

                            i = j;
                            table = table2;
                            kmer_hash = kmer_hash_2;
                        }

                        __builtin_prefetch(table, 1, 1);

                        for (; i != k; ++i) {

                            //__sync_fetch_and_or(table + ((kmer_hash & MASK_BITS_BLOCK) >> 6), 1ULL << (kmer_hash & 0x3fULL));
                            table[(kmer_hash & MASK_BITS_BLOCK) >> 6] |= 1ULL << (kmer_hash & 0x3fULL);
                            kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
                        }
                    }
                    else {

                        if (popcnt(table2, NB_ELEM_BLOCK * sizeof(uint64_t)) < popcnt(table, NB_ELEM_BLOCK * sizeof(uint64_t))){

                            int tmp = i;
                            i = j;
                            j = tmp;

                            uint64_t tmp_size_t = kmer_hash;
                            kmer_hash = kmer_hash_2;
                            kmer_hash_2 = tmp_size_t;

                            uint64_t* tmp_ptr = table;
                            table = table2;
                            table2 = tmp_ptr;
                        }

                        __builtin_prefetch(table, 1, 1);

                        for (; i != k; ++i) {

                            __sync_fetch_and_or(table + ((kmer_hash & MASK_BITS_BLOCK) >> 6), 1ULL << (kmer_hash & 0x3fULL));
                            //table[(kmer_hash & MASK_BITS_BLOCK) >> 6] |= 1ULL << (kmer_hash & 0x3fULL);
                            kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
                        }

                        __builtin_prefetch(table2, 0, 1);

                        for (; j != k; ++j) {

                            if ((table2[(kmer_hash_2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_2 & 0x3fULL))) == 0) break;
                            kmer_hash_2 = (kmer_hash_2 * 49157) % (1610612741ULL);
                        }
                    }

                    return j != k;
                }
            }

            return false;
        }

        inline void insert(uint64_t kmer_hash, const uint64_t min_hash){

            search_and_insert(kmer_hash, min_hash, false);
        }

        bool WriteBloomFilter(FILE *fp) {

            if (fwrite(&size_table_, sizeof(size_table_), 1, fp) != 1) return false;
            if (fwrite(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
            if (fwrite(&k_, sizeof(k_), 1, fp) != 1) return false;

            if (fwrite(table_, sizeof(uint64_t), NB_ELEM_BLOCK * blocks_, fp) != (NB_ELEM_BLOCK * blocks_)) return false;

            return true;
        }

        bool ReadBloomFilter(FILE *fp) {

            clear();

            if (fread(&size_table_, sizeof(size_table_), 1, fp) != 1) return false;
            if (fread(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
            if (fread(&k_, sizeof(k_), 1, fp) != 1) return false;

            init_table();

            if (fread(table_, sizeof(uint64_t), NB_ELEM_BLOCK * blocks_, fp) != (NB_ELEM_BLOCK * blocks_)) return false;

            return true;
        }

        void clear() {

            if (table_ != NULL){

                free(table_);
                table_ = NULL;
            }

            size_table_ = 0;
            blocks_ = 0;
            k_ = 0;
        }

        void get(BlockedBloomFilter& bf) {

            clear();

            table_ = bf.table_;
            size_table_ = bf.size_table_;
            blocks_ = bf.blocks_;
            k_ = bf.k_;
            fast_div_ = bf.fast_div_;

            bf.table_ = NULL;
        }

        inline uint64_t getNbBlocks() const { return blocks_; }

        inline const uint64_t* getTable_ptr() const { return table_; }

    private:

        void init_table(){

            fast_div_ = libdivide::divider<uint64_t>(blocks_);

            posix_memalign((void**)&table_, 64, NB_ELEM_BLOCK * blocks_* sizeof(table_[0]));
            memset(table_, 0, NB_ELEM_BLOCK * blocks_ * sizeof(table_[0]));
        }

        inline double fpp(size_t bits, int k) const {

            return pow(1-exp(-((double)k)/((double)bits)),(double)k);
        }
};

#endif // BFG_BLOCKEDBLOOMFILTER_HPP
