#ifndef BFG_BLOCKEDBLOOMFILTER_HPP
#define BFG_BLOCKEDBLOOMFILTER_HPP

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstring>

#include <vector>

#include <algorithm>
#include <fstream>
#include <random>

#include "KmerHashTable.h"
#include "Kmer.hpp"
#include "RepHash.hpp"

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

class BlockedBloomFilter {

    public:

        typedef std::pair<uint64_t*, uint64_t*> BBF_Blocks;

        BlockedBloomFilter() : table_(NULL), size_table_(0), blocks_(0), k_(0), fast_div_() {}

        BlockedBloomFilter(size_t nb_elem, size_t bits_per_elem) : table_(NULL), size_table_(0), blocks_(0), k_(0), fast_div_() {

            size_table_ = ((bits_per_elem * nb_elem + MASK_BITS_BLOCK) / NB_BITS_BLOCK) * NB_BITS_BLOCK;
            blocks_ = size_table_ / NB_BITS_BLOCK;

            init_table();

            k_ = (int) (bits_per_elem * log(2));
            if (fpp(bits_per_elem, k_) >= fpp(bits_per_elem, k_+1)) k_++;

            if (k_ > 16){

                std::cerr << "BlockedBloomFilter(): The AVX2 Blocked Bloom filter does not support more than 16 hash functions." << std::endl;
                std::cerr << "Either use less bits per element to insert or recompile code of Bifrost without AVX2 activated." << std::endl;

                clear();
            }
            else {

                uint64_t hashes_mask[4] __attribute__((aligned(32))) = {0, 0, 0, 0};

                for (int k = 0; k != k_; ++k) hashes_mask[k/4] = (hashes_mask[k/4] << 16) | 0xffff;

                mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);
            }
        }

        ~BlockedBloomFilter() {

            clear();
        }

        inline BBF_Blocks getBlock(uint64_t min_hash) const{

            uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

            min_hash -= (min_hash / fast_div_) * blocks_;
            min_hash_2 -= (min_hash_2 / fast_div_) * blocks_;

            return std::make_pair(table_ + NB_ELEM_BLOCK * min_hash, table_ + NB_ELEM_BLOCK * min_hash_2);
        }

        int contains(const uint64_t (&kmer_hash)[4], const uint64_t min_hash, bool (&pres)[4], const int limit) const {

            int cpt = 0;

            __m256i table_gather;

            //Gather and compare
            const uint16_t* table = reinterpret_cast<const uint16_t*>(table_ + ((min_hash - (min_hash / fast_div_) * blocks_) * NB_ELEM_BLOCK));

            for (uint8_t j = 0; (j != 4) && (cpt != limit); ++j){

                const uint64_t hashes_div[4] __attribute__((aligned(32))) { kmer_hash[j], kmer_hash[j] + min_hash, kmer_hash[j] + min_hash + min_hash,
                                                                            kmer_hash[j] + min_hash + min_hash + min_hash};

                const __m256i km_hashes = _mm256_and_si256(_mm256_set_epi64x(hashes_div[0], hashes_div[1], hashes_div[2], hashes_div[3]), mask_and_div);
                const __m256i h_shift = _mm256_and_si256(km_hashes, mask_and_mod);

                const __m256i hash_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
                const __m256i hash_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

                const __m256i hash_gather = _mm256_and_si256(mask_h, _mm256_or_si256(hash_gather_lsb, hash_gather_msb));

                _mm256_stream_si256((__m256i*)hashes_div, _mm256_srli_epi16(km_hashes, 4));

                const uint16_t* hashes_div16 = reinterpret_cast<const uint16_t*>(hashes_div);

                if (k_ <= 4){

                    table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                    table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
                }

                if ((k_ > 4) && (k_ <= 8)){

                    table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0,
                                                    table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                                    table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
                }

                if ((k_ > 8) && (k_ <= 12)){

                   table_gather = _mm256_set_epi16(0, 0, 0, 0,
                                                   table[hashes_div16[11]], table[hashes_div16[10]], table[hashes_div16[9]], table[hashes_div16[8]],
                                                   table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                                   table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
                }

                if ((k_ > 12) && (k_ <= 16)){

                   table_gather = _mm256_set_epi16(table[hashes_div16[15]], table[hashes_div16[14]], table[hashes_div16[13]], table[hashes_div16[12]],
                                                   table[hashes_div16[11]], table[hashes_div16[10]], table[hashes_div16[9]], table[hashes_div16[8]],
                                                   table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                                   table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
                }

                const __m256i xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather, hash_gather), hash_gather);

                cpt += (pres[j] = (_mm256_testz_si256(xor_m256i, xor_m256i) != 0));
            }

            if (cpt != limit){

                const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

                const uint16_t* table2 = reinterpret_cast<const uint16_t*>(table_ + ((min_hash_2 - (min_hash_2 / fast_div_) * blocks_) * NB_ELEM_BLOCK));

                for (uint8_t j = 0; j != 4; ++j){

                    if (!pres[j]){

                        const uint64_t hashes_div[4] __attribute__((aligned(32))) { kmer_hash[j], kmer_hash[j] + min_hash, kmer_hash[j] + min_hash + min_hash,
                                                                                    kmer_hash[j] + min_hash + min_hash + min_hash};

                        const __m256i km_hashes = _mm256_and_si256(_mm256_set_epi64x(hashes_div[0], hashes_div[1], hashes_div[2], hashes_div[3]), mask_and_div);
                        const __m256i h_shift = _mm256_and_si256(km_hashes, mask_and_mod);

                        const __m256i hash_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
                        const __m256i hash_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

                        const __m256i hash_gather = _mm256_and_si256(mask_h, _mm256_or_si256(hash_gather_lsb, hash_gather_msb));

                        _mm256_stream_si256((__m256i*)hashes_div, _mm256_srli_epi16(km_hashes, 4));

                        const uint16_t* hashes_div16 = reinterpret_cast<const uint16_t*>(hashes_div);

                        if (k_ <= 4){

                            table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                            table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
                        }

                        if ((k_ > 4) && (k_ <= 8)){

                            table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0,
                                                            table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                                            table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
                        }

                        if ((k_ > 8) && (k_ <= 12)){

                           table_gather = _mm256_set_epi16(0, 0, 0, 0,
                                                           table2[hashes_div16[11]], table2[hashes_div16[10]], table2[hashes_div16[9]], table2[hashes_div16[8]],
                                                           table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                                           table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
                        }

                        if ((k_ > 12) && (k_ <= 16)){

                           table_gather = _mm256_set_epi16(table2[hashes_div16[15]], table2[hashes_div16[14]], table2[hashes_div16[13]], table2[hashes_div16[12]],
                                                           table2[hashes_div16[11]], table2[hashes_div16[10]], table2[hashes_div16[9]], table2[hashes_div16[8]],
                                                           table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                                           table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
                        }

                        const __m256i xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather, hash_gather), hash_gather);

                        cpt += (pres[j] = (_mm256_testz_si256(xor_m256i, xor_m256i) != 0));

                        if (cpt == limit) break;
                    }
                }
            }

            return cpt;
        }

        bool contains(const uint64_t kmer_hash, const uint64_t min_hash) const {

            uint64_t hashes_div[4] __attribute__((aligned(32)));

            hashes_div[0] = kmer_hash;
            hashes_div[1] = hashes_div[0] + min_hash;
            hashes_div[2] = hashes_div[1] + min_hash;
            hashes_div[3] = hashes_div[2] + min_hash;

            const __m256i km_hashes = _mm256_and_si256(_mm256_set_epi64x(hashes_div[0], hashes_div[1], hashes_div[2], hashes_div[3]), mask_and_div);
            const __m256i h_shift = _mm256_and_si256(km_hashes, mask_and_mod);

            const __m256i hash_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
            const __m256i hash_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));
            const __m256i hash_gather = _mm256_and_si256(mask_h, _mm256_or_si256(hash_gather_lsb, hash_gather_msb));

            _mm256_stream_si256((__m256i*)hashes_div, _mm256_srli_epi16(km_hashes, 4));

            const uint16_t* hashes_div16 = reinterpret_cast<const uint16_t*>(hashes_div);

            //Gather and compare
            const uint16_t* table = reinterpret_cast<const uint16_t*>(table_ + ((min_hash - (min_hash / fast_div_) * blocks_) * NB_ELEM_BLOCK));

            __m256i table_gather;

            if (k_ <= 4){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            if ((k_ > 4) && (k_ <= 8)){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0,
                                                table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                                table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            if ((k_ > 8) && (k_ <= 12)){

               table_gather = _mm256_set_epi16(0, 0, 0, 0,
                                               table[hashes_div16[11]], table[hashes_div16[10]], table[hashes_div16[9]], table[hashes_div16[8]],
                                               table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                               table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            if ((k_ > 12) && (k_ <= 16)){

               table_gather = _mm256_set_epi16(table[hashes_div16[15]], table[hashes_div16[14]], table[hashes_div16[13]], table[hashes_div16[12]],
                                               table[hashes_div16[11]], table[hashes_div16[10]], table[hashes_div16[9]], table[hashes_div16[8]],
                                               table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                               table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            const __m256i xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather, hash_gather), hash_gather);

            if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0) return true;

            const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

            const uint16_t* table2 = reinterpret_cast<const uint16_t*>(table_ + ((min_hash_2 - (min_hash_2 / fast_div_) * blocks_) * NB_ELEM_BLOCK));

            if (k_ <= 4){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            if ((k_ > 4) && (k_ <= 8)){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0,
                                                table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                                table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            if ((k_ > 8) && (k_ <= 12)){

               table_gather = _mm256_set_epi16(0, 0, 0, 0,
                                               table2[hashes_div16[11]], table2[hashes_div16[10]], table2[hashes_div16[9]], table2[hashes_div16[8]],
                                               table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                               table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            if ((k_ > 12) && (k_ <= 16)){

               table_gather = _mm256_set_epi16(table2[hashes_div16[15]], table2[hashes_div16[14]], table2[hashes_div16[13]], table2[hashes_div16[12]],
                                               table2[hashes_div16[11]], table2[hashes_div16[10]], table2[hashes_div16[9]], table2[hashes_div16[8]],
                                               table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                               table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            const __m256i xor_m256i_2 = _mm256_xor_si256(_mm256_and_si256(table_gather, hash_gather), hash_gather);

            if (_mm256_testz_si256(xor_m256i_2, xor_m256i_2) != 0) return true;

            return false;
        }

        inline bool contains(const uint64_t kmer_hash, const uint64_t min_hash, const BBF_Blocks block_ptr) const {

            return (contains_block(kmer_hash, min_hash, block_ptr) != 0);
        }

        size_t contains_block(const uint64_t kmer_hash, const uint64_t min_hash, const BBF_Blocks block_ptr) const {

            uint64_t hashes_div[4] __attribute__((aligned(32)));

            hashes_div[0] = kmer_hash;
            hashes_div[1] = hashes_div[0] + min_hash;
            hashes_div[2] = hashes_div[1] + min_hash;
            hashes_div[3] = hashes_div[2] + min_hash;

            const __m256i km_hashes = _mm256_and_si256(_mm256_set_epi64x(hashes_div[0], hashes_div[1], hashes_div[2], hashes_div[3]), mask_and_div);
            const __m256i h_shift = _mm256_and_si256(km_hashes, mask_and_mod);

            const __m256i hash_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
            const __m256i hash_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));
            const __m256i hash_gather = _mm256_and_si256(mask_h, _mm256_or_si256(hash_gather_lsb, hash_gather_msb));

            _mm256_stream_si256((__m256i*)hashes_div, _mm256_srli_epi16(km_hashes, 4));

            const uint16_t* hashes_div16 = reinterpret_cast<const uint16_t*>(hashes_div);

            //Gather and compare
            const uint16_t* table = reinterpret_cast<const uint16_t*>(block_ptr.first);

            __m256i table_gather;

            if (k_ <= 4){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            if ((k_ > 4) && (k_ <= 8)){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0,
                                                table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                                table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            if ((k_ > 8) && (k_ <= 12)){

               table_gather = _mm256_set_epi16(0, 0, 0, 0,
                                               table[hashes_div16[11]], table[hashes_div16[10]], table[hashes_div16[9]], table[hashes_div16[8]],
                                               table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                               table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            if ((k_ > 12) && (k_ <= 16)){

               table_gather = _mm256_set_epi16(table[hashes_div16[15]], table[hashes_div16[14]], table[hashes_div16[13]], table[hashes_div16[12]],
                                               table[hashes_div16[11]], table[hashes_div16[10]], table[hashes_div16[9]], table[hashes_div16[8]],
                                               table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                               table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            const __m256i xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather, hash_gather), hash_gather);

            if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0) return 1;

            const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

            const uint16_t* table2 = reinterpret_cast<const uint16_t*>(block_ptr.second);

            if (k_ <= 4){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            if ((k_ > 4) && (k_ <= 8)){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0,
                                                table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                                table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            if ((k_ > 8) && (k_ <= 12)){

               table_gather = _mm256_set_epi16(0, 0, 0, 0,
                                               table2[hashes_div16[11]], table2[hashes_div16[10]], table2[hashes_div16[9]], table2[hashes_div16[8]],
                                               table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                               table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            if ((k_ > 12) && (k_ <= 16)){

               table_gather = _mm256_set_epi16(table2[hashes_div16[15]], table2[hashes_div16[14]], table2[hashes_div16[13]], table2[hashes_div16[12]],
                                               table2[hashes_div16[11]], table2[hashes_div16[10]], table2[hashes_div16[9]], table2[hashes_div16[8]],
                                               table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                               table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            const __m256i xor_m256i_2 = _mm256_xor_si256(_mm256_and_si256(table_gather, hash_gather), hash_gather);

            if (_mm256_testz_si256(xor_m256i_2, xor_m256i_2) != 0) return 2;

            return 0;
        }

        bool search_and_insert(const uint64_t kmer_hash, const uint64_t min_hash, const bool multi_threaded = false) {

            uint64_t hashes_div[4] __attribute__((aligned(32)));

            hashes_div[0] = kmer_hash;
            hashes_div[1] = hashes_div[0] + min_hash;
            hashes_div[2] = hashes_div[1] + min_hash;
            hashes_div[3] = hashes_div[2] + min_hash;

            const __m256i km_hashes = _mm256_and_si256(_mm256_set_epi64x(hashes_div[0], hashes_div[1], hashes_div[2], hashes_div[3]), mask_and_div);
            const __m256i h_shift = _mm256_and_si256(km_hashes, mask_and_mod);

            const __m256i hash_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
            const __m256i hash_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));
            const __m256i hash_gather = _mm256_and_si256(mask_h, _mm256_or_si256(hash_gather_lsb, hash_gather_msb));

            _mm256_stream_si256((__m256i*)hashes_div, _mm256_srli_epi16(km_hashes, 4));

            const uint16_t* hashes_div16 = reinterpret_cast<const uint16_t*>(hashes_div);

            //Gather and compare
            uint16_t* table = reinterpret_cast<uint16_t*>(table_ + ((min_hash - (min_hash / fast_div_) * blocks_) * NB_ELEM_BLOCK));

            __m256i table_gather;

            if (k_ <= 4){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            if ((k_ > 4) && (k_ <= 8)){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0,
                                                table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                                table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            if ((k_ > 8) && (k_ <= 12)){

               table_gather = _mm256_set_epi16(0, 0, 0, 0,
                                               table[hashes_div16[11]], table[hashes_div16[10]], table[hashes_div16[9]], table[hashes_div16[8]],
                                               table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                               table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            if ((k_ > 12) && (k_ <= 16)){

               table_gather = _mm256_set_epi16(table[hashes_div16[15]], table[hashes_div16[14]], table[hashes_div16[13]], table[hashes_div16[12]],
                                               table[hashes_div16[11]], table[hashes_div16[10]], table[hashes_div16[9]], table[hashes_div16[8]],
                                               table[hashes_div16[7]], table[hashes_div16[6]], table[hashes_div16[5]], table[hashes_div16[4]],
                                               table[hashes_div16[3]], table[hashes_div16[2]], table[hashes_div16[1]], table[hashes_div16[0]]);
            }

            const __m256i xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather, hash_gather), hash_gather);
            if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0) return false;

            const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

            uint16_t* table2 = reinterpret_cast<uint16_t*>(table_ + ((min_hash_2 - (min_hash_2 / fast_div_) * blocks_) * NB_ELEM_BLOCK));

            if (k_ <= 4){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            if ((k_ > 4) && (k_ <= 8)){

                table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0,
                                                table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                                table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            if ((k_ > 8) && (k_ <= 12)){

               table_gather = _mm256_set_epi16(0, 0, 0, 0,
                                               table2[hashes_div16[11]], table2[hashes_div16[10]], table2[hashes_div16[9]], table2[hashes_div16[8]],
                                               table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                               table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            if ((k_ > 12) && (k_ <= 16)){

               table_gather = _mm256_set_epi16(table2[hashes_div16[15]], table2[hashes_div16[14]], table2[hashes_div16[13]], table2[hashes_div16[12]],
                                               table2[hashes_div16[11]], table2[hashes_div16[10]], table2[hashes_div16[9]], table2[hashes_div16[8]],
                                               table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                               table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
            }

            const __m256i xor_m256i_2 = _mm256_xor_si256(_mm256_and_si256(table_gather, hash_gather), hash_gather);
            if (_mm256_testz_si256(xor_m256i_2, xor_m256i_2) != 0) return false;

            uint16_t hashes_mod[16] __attribute__((aligned(32)));
            _mm256_stream_si256(reinterpret_cast<__m256i*>(hashes_mod), hash_gather);

            if (!multi_threaded){

                if (popcnt(table2, NB_ELEM_BLOCK * sizeof(uint64_t)) <= popcnt(table, NB_ELEM_BLOCK * sizeof(uint64_t))) table = table2;

                switch(k_){
                    case 16: table[hashes_div16[15]] |= hashes_mod[15];
                    case 15: table[hashes_div16[14]] |= hashes_mod[14];
                    case 14: table[hashes_div16[13]] |= hashes_mod[13];
                    case 13: table[hashes_div16[12]] |= hashes_mod[12];
                    case 12: table[hashes_div16[11]] |= hashes_mod[11];
                    case 11: table[hashes_div16[10]] |= hashes_mod[10];
                    case 10: table[hashes_div16[9]] |= hashes_mod[9];
                    case 9: table[hashes_div16[8]] |= hashes_mod[8];
                    case 8: table[hashes_div16[7]] |= hashes_mod[7];
                    case 7: table[hashes_div16[6]] |= hashes_mod[6];
                    case 6: table[hashes_div16[5]] |= hashes_mod[5];
                    case 5: table[hashes_div16[4]] |= hashes_mod[4];
                    case 4: table[hashes_div16[3]] |= hashes_mod[3];
                    case 3: table[hashes_div16[2]] |= hashes_mod[2];
                    case 2: table[hashes_div16[1]] |= hashes_mod[1];
                    case 1: table[hashes_div16[0]] |= hashes_mod[0];
                }
            }
            else {

                if (popcnt(table2, NB_ELEM_BLOCK * sizeof(uint64_t)) <= popcnt(table, NB_ELEM_BLOCK * sizeof(uint64_t))){

                    uint16_t* tmp_ptr = table;
                    table = table2;
                    table2 = tmp_ptr;
                }

                switch(k_){
                    case 16: __sync_fetch_and_or(&table[hashes_div16[15]], hashes_mod[15]);
                    case 15: __sync_fetch_and_or(&table[hashes_div16[14]], hashes_mod[14]);
                    case 14: __sync_fetch_and_or(&table[hashes_div16[13]], hashes_mod[13]);
                    case 13: __sync_fetch_and_or(&table[hashes_div16[12]], hashes_mod[12]);
                    case 12: __sync_fetch_and_or(&table[hashes_div16[11]], hashes_mod[11]);
                    case 11: __sync_fetch_and_or(&table[hashes_div16[10]], hashes_mod[10]);
                    case 10: __sync_fetch_and_or(&table[hashes_div16[9]], hashes_mod[9]);
                    case 9: __sync_fetch_and_or(&table[hashes_div16[8]], hashes_mod[8]);
                    case 8: __sync_fetch_and_or(&table[hashes_div16[7]], hashes_mod[7]);
                    case 7: __sync_fetch_and_or(&table[hashes_div16[6]], hashes_mod[6]);
                    case 6: __sync_fetch_and_or(&table[hashes_div16[5]], hashes_mod[5]);
                    case 5: __sync_fetch_and_or(&table[hashes_div16[4]], hashes_mod[4]);
                    case 4: __sync_fetch_and_or(&table[hashes_div16[3]], hashes_mod[3]);
                    case 3: __sync_fetch_and_or(&table[hashes_div16[2]], hashes_mod[2]);
                    case 2: __sync_fetch_and_or(&table[hashes_div16[1]], hashes_mod[1]);
                    case 1: __sync_fetch_and_or(&table[hashes_div16[0]], hashes_mod[0]);
                }

                if (k_ <= 4){

                    table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                    table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
                }

                if ((k_ > 4) && (k_ <= 8)){

                    table_gather = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0,
                                                    table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                                    table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
                }

                if ((k_ > 8) && (k_ <= 12)){

                   table_gather = _mm256_set_epi16(0, 0, 0, 0,
                                                   table2[hashes_div16[11]], table2[hashes_div16[10]], table2[hashes_div16[9]], table2[hashes_div16[8]],
                                                   table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                                   table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
                }

                if ((k_ > 12) && (k_ <= 16)){

                   table_gather = _mm256_set_epi16(table2[hashes_div16[15]], table2[hashes_div16[14]], table2[hashes_div16[13]], table2[hashes_div16[12]],
                                                   table2[hashes_div16[11]], table2[hashes_div16[10]], table2[hashes_div16[9]], table2[hashes_div16[8]],
                                                   table2[hashes_div16[7]], table2[hashes_div16[6]], table2[hashes_div16[5]], table2[hashes_div16[4]],
                                                   table2[hashes_div16[3]], table2[hashes_div16[2]], table2[hashes_div16[1]], table2[hashes_div16[0]]);
                }

                const __m256i xor_m256i_3 = _mm256_xor_si256(_mm256_and_si256(table_gather, hash_gather), hash_gather);

                if (_mm256_testz_si256(xor_m256i_3, xor_m256i_3) != 0) return false;
            }

            return true;
        }

        inline void insert(const uint64_t kmer_hash, const uint64_t min_hash){

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

            uint64_t hashes_mask[4] __attribute__((aligned(32))) = {0, 0, 0, 0};

            for (int k = 0; k != k_; ++k) hashes_mask[k/4] = (hashes_mask[k/4] << 16) | 0xffff;

            mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

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
            mask_h = _mm256_setzero_si256();
        }

        void get(BlockedBloomFilter& bf) {

            clear();

            table_ = bf.table_;
            size_table_ = bf.size_table_;
            blocks_ = bf.blocks_;
            k_ = bf.k_;
            mask_h = bf.mask_h;
            fast_div_ = bf.fast_div_;

            bf.table_ = NULL;
        }

        inline uint64_t getNbBlocks() const { return blocks_; }

        inline const uint64_t* getTable_ptr() const { return table_; }

    private:

        uint64_t* table_; //Bit array

        uint64_t size_table_; //Size of bit array (in bits)
        uint64_t blocks_; //Nb blocks

        int k_; //Nb hash functions

        __m256i mask_h; // Mask for hash functions that we do not want (>k);

        static const __m256i mask_and_div; // All MASK_BITS_BLOCK LSB of each 16 bits word
        static const __m256i mask_and_mod; // All 4 LSB of each 16 bits word
        static const __m256i one2shift_lsb; // All 1 LSB of each 32 bits word
        static const __m256i one2shift_msb; // Set the 17th LSB bit of each 32 bits word
        static const __m256i mask_lsb; // All 16 LSB bits of each 32 bits word

        libdivide::divider<uint64_t> fast_div_; // fast division

        void init_table(){

            fast_div_ = libdivide::divider<uint64_t>(blocks_);

            posix_memalign((void**)&table_, 8, NB_ELEM_BLOCK * blocks_* sizeof(uint64_t));
            memset(table_, 0, NB_ELEM_BLOCK * blocks_ * sizeof(uint64_t));
        }

        inline double fpp(size_t bits, int k) const {

            return pow(1-exp(-((double)k)/((double)bits)),(double)k);
        }
};

#else

class BlockedBloomFilter {

    public:

        typedef std::pair<uint64_t*, uint64_t*> BBF_Blocks;

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

        inline BBF_Blocks getBlock(uint64_t min_hash) const{

            uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

            min_hash -= (min_hash / fast_div_) * blocks_;
            min_hash_2 -= (min_hash_2 / fast_div_) * blocks_;

            return std::make_pair(table_ + NB_ELEM_BLOCK * min_hash, table_ + NB_ELEM_BLOCK * min_hash_2);
        }

        int contains(const uint64_t (&kmer_hash)[4], const uint64_t min_hash, bool (&pres)[4], const int limit) const{

            int cpt = 0;

            const int k = k_;

            uint64_t* table = table_ + ((min_hash - (min_hash / fast_div_) * blocks_) * NB_ELEM_BLOCK);

            __builtin_prefetch(table, 0, 1);

            for (size_t i, j = 0; (j != 4) && (cpt != limit); ++j){

                uint64_t kmer_hash_tmp = kmer_hash[j];

                for (i = 0; i != k; ++i) {

                    if ((table[(kmer_hash_tmp & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_tmp & 0x3fULL))) == 0) break;
                    kmer_hash_tmp = (kmer_hash_tmp * 49157) % (1610612741ULL);
                }

                cpt += (pres[j] = (i == k));
            }

            if (cpt != limit){

                const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

                table = table_ + ((min_hash_2 - (min_hash_2 / fast_div_) * blocks_) * NB_ELEM_BLOCK);

                __builtin_prefetch(table, 0, 1);

                for (size_t i, j = 0; (j != 4) && (cpt != limit); ++j){

                    if (!pres[j]){

                        uint64_t kmer_hash_tmp = kmer_hash[j];

                        for (i = 0; i != k; ++i) {

                            if ((table[(kmer_hash_tmp & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_tmp & 0x3fULL))) == 0) break;
                            kmer_hash_tmp = (kmer_hash_tmp * 49157) % (1610612741ULL);
                        }

                        cpt += (pres[j] = (i == k));
                    }
                }
            }

            return cpt;
        }

        bool contains(uint64_t kmer_hash, const uint64_t min_hash) const {

            int i = 0;

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

                table = table_ + ((min_hash_2 - (min_hash_2 / fast_div_) * blocks_) * NB_ELEM_BLOCK);

                __builtin_prefetch(table, 0, 1);

                for (i = 0; i != k; ++i) {

                    if ((table[(kmer_hash_2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_2 & 0x3fULL))) == 0) return false;
                    kmer_hash_2 = (kmer_hash_2 * 49157) % (1610612741ULL);
                }

                return true;
            }

            return true;
        }

        inline bool contains(const uint64_t kmer_hash, const uint64_t min_hash, const BBF_Blocks block_ptr) const {

            return (contains_block(kmer_hash, min_hash, block_ptr) != 0);
        }

        size_t contains_block(uint64_t kmer_hash, const uint64_t min_hash, const BBF_Blocks block_ptr) const {

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

                    if ((block_ptr.second[(kmer_hash_2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_2 & 0x3fULL))) == 0) return 0;
                    kmer_hash_2 = (kmer_hash_2 * 49157) % (1610612741ULL);
                }

                return 2;
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

                        if (popcnt(table2, NB_ELEM_BLOCK * sizeof(uint64_t)) <= popcnt(table, NB_ELEM_BLOCK * sizeof(uint64_t))){

                            i = j;
                            table = table2;
                            kmer_hash = kmer_hash_2;
                        }

                        __builtin_prefetch(table, 1, 1);

                        for (; i != k; ++i) {

                            table[(kmer_hash & MASK_BITS_BLOCK) >> 6] |= 1ULL << (kmer_hash & 0x3fULL);
                            kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
                        }

                        return true;
                    }
                    else {

                        if (popcnt(table2, NB_ELEM_BLOCK * sizeof(uint64_t)) <= popcnt(table, NB_ELEM_BLOCK * sizeof(uint64_t))){

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
                            kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
                        }

                        __builtin_prefetch(table2, 0, 1);

                        for (; j != k; ++j) {

                            if ((table2[(kmer_hash_2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_2 & 0x3fULL))) == 0) return true;
                            kmer_hash_2 = (kmer_hash_2 * 49157) % (1610612741ULL);
                        }

                        return false;
                    }
                }
            }

            return false;
        }

        inline void insert(const uint64_t kmer_hash, const uint64_t min_hash){

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

        uint64_t* table_; //Bit array

        uint64_t size_table_; //Size of bit array (in bits)
        uint64_t blocks_; //Nb blocks
        int k_; //Nb hash functions

        libdivide::divider<uint64_t> fast_div_; // fast division

        void init_table(){

            fast_div_ = libdivide::divider<uint64_t>(blocks_);

            posix_memalign((void**)&table_, 8, NB_ELEM_BLOCK * blocks_* sizeof(table_[0]));
            memset(table_, 0, NB_ELEM_BLOCK * blocks_ * sizeof(table_[0]));
        }

        inline double fpp(size_t bits, int k) const {

            return pow(1-exp(-((double)k)/((double)bits)),(double)k);
        }
};

#endif

#endif // BFG_BLOCKEDBLOOMFILTER_HPP
