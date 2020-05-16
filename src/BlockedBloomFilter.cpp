#include "BlockedBloomFilter.hpp"

#if defined(__AVX2__)

BlockedBloomFilter::BlockedBloomFilter() : table_(nullptr), blocks_(0), k_(0) {

    hashes_mask[0] = 0;
    hashes_mask[1] = 0;
    hashes_mask[2] = 0;
    hashes_mask[3] = 0;
}

BlockedBloomFilter::BlockedBloomFilter(size_t nb_elem, size_t bits_per_elem) : table_(nullptr), blocks_(0), k_(0) {

    if ((nb_elem != 0) && (bits_per_elem != 0)){

        blocks_ = (bits_per_elem * nb_elem + MASK_BITS_BLOCK) / NB_BITS_BLOCK;
        k_ = (int) (bits_per_elem * log(2));

        hashes_mask[0] = 0;
        hashes_mask[1] = 0;
        hashes_mask[2] = 0;
        hashes_mask[3] = 0;

        if (fpp(bits_per_elem, k_) >= fpp(bits_per_elem, k_+1)) ++k_;

        if (k_ > 16){

            std::cerr << "BlockedBloomFilter(): The AVX2 Blocked Bloom filter does not support more than 16 hash functions." << std::endl;
            std::cerr << "Either use less bits per element to insert or recompile code of Bifrost with AVX2 deactivated." << std::endl;

            clear();
        }
        else {

            //uint64_t hashes_mask[4] __attribute__((aligned(32))) = {0, 0, 0, 0};

            for (int k = 0; k != k_; ++k) hashes_mask[k/4] = (hashes_mask[k/4] << 16) | 0xffff;

            //mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);
        }

        init_table();
    }
}

BlockedBloomFilter::BlockedBloomFilter(const BlockedBloomFilter& o) : table_(nullptr), blocks_(o.blocks_), k_(o.k_), fast_div_(o.fast_div_) {

    std::memcpy(hashes_mask, o.hashes_mask, 4 * sizeof(uint64_t));

    if (blocks_ != 0){

        init_table();

        for (uint64_t i = 0; i != blocks_; ++i){

            memcpy(&(table_[i].block), &(o.table_[i].block), NB_ELEM_BLOCK * sizeof(uint64_t));
        }
    }
}

BlockedBloomFilter::BlockedBloomFilter(BlockedBloomFilter&& o) : table_(o.table_), blocks_(o.blocks_), k_(o.k_), fast_div_(o.fast_div_) {

    std::memcpy(hashes_mask, o.hashes_mask, 4 * sizeof(uint64_t));

    o.table_ = nullptr;

    o.clear();
}

BlockedBloomFilter::~BlockedBloomFilter() {

    clear();
}

BlockedBloomFilter& BlockedBloomFilter::operator=(const BlockedBloomFilter& o) {

    clear();

    table_ = nullptr;
    blocks_ = o.blocks_;
    k_ = o.k_;
    fast_div_ = o.fast_div_;

    std::memcpy(hashes_mask, o.hashes_mask, 4 * sizeof(uint64_t));

    init_table();

    for (uint64_t i = 0; i != blocks_; ++i){

        memcpy(&(table_[i].block), &(o.table_[i].block), NB_ELEM_BLOCK * sizeof(uint64_t));
    }

    return *this;
}

BlockedBloomFilter& BlockedBloomFilter::operator=(BlockedBloomFilter&& o) {

    if (this != &o) {

        clear();

        table_ = o.table_;
        blocks_ = o.blocks_;
        k_ = o.k_;
        fast_div_ = o.fast_div_;

        std::memcpy(hashes_mask, o.hashes_mask, 4 * sizeof(uint64_t));

        o.table_ = nullptr;

        o.clear();
    }

    return *this;
}

int BlockedBloomFilter::contains(const uint64_t (&kmer_hash)[4], const uint64_t min_hash, bool (&pres)[4], const int limit) const {

    int cpt = 0;

    __m256i table_gather;

    const __m256i mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

    //Gather and compare
    const uint16_t* table = reinterpret_cast<const uint16_t*>(table_[min_hash - (min_hash / fast_div_) * blocks_].block);

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

        const uint16_t* table2 = reinterpret_cast<const uint16_t*>(table_[min_hash_2 - (min_hash_2 / fast_div_) * blocks_].block);

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

bool BlockedBloomFilter::contains(const uint64_t kmer_hash, const uint64_t min_hash) const {

    uint64_t hashes_div[4] __attribute__((aligned(32)));

    hashes_div[0] = kmer_hash;
    hashes_div[1] = hashes_div[0] + min_hash;
    hashes_div[2] = hashes_div[1] + min_hash;
    hashes_div[3] = hashes_div[2] + min_hash;

    const __m256i km_hashes = _mm256_and_si256(_mm256_set_epi64x(hashes_div[0], hashes_div[1], hashes_div[2], hashes_div[3]), mask_and_div);
    const __m256i h_shift = _mm256_and_si256(km_hashes, mask_and_mod);

    const __m256i hash_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
    const __m256i hash_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

    const __m256i mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

    const __m256i hash_gather = _mm256_and_si256(mask_h, _mm256_or_si256(hash_gather_lsb, hash_gather_msb));

    _mm256_stream_si256((__m256i*)hashes_div, _mm256_srli_epi16(km_hashes, 4));

    const uint16_t* hashes_div16 = reinterpret_cast<const uint16_t*>(hashes_div);

    //Gather and compare
    const uint16_t* table = reinterpret_cast<const uint16_t*>(table_[min_hash - (min_hash / fast_div_) * blocks_].block);

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

    const uint16_t* table2 = reinterpret_cast<const uint16_t*>(table_[min_hash_2 - (min_hash_2 / fast_div_) * blocks_].block);

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

size_t BlockedBloomFilter::contains_block(const uint64_t kmer_hash, const uint64_t min_hash, const BBF_Blocks blockIDs) const {

    uint64_t hashes_div[4] __attribute__((aligned(32)));

    hashes_div[0] = kmer_hash;
    hashes_div[1] = hashes_div[0] + min_hash;
    hashes_div[2] = hashes_div[1] + min_hash;
    hashes_div[3] = hashes_div[2] + min_hash;

    const __m256i km_hashes = _mm256_and_si256(_mm256_set_epi64x(hashes_div[0], hashes_div[1], hashes_div[2], hashes_div[3]), mask_and_div);
    const __m256i h_shift = _mm256_and_si256(km_hashes, mask_and_mod);

    const __m256i hash_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
    const __m256i hash_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

    const __m256i mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

    const __m256i hash_gather = _mm256_and_si256(mask_h, _mm256_or_si256(hash_gather_lsb, hash_gather_msb));

    _mm256_stream_si256((__m256i*)hashes_div, _mm256_srli_epi16(km_hashes, 4));

    const uint16_t* hashes_div16 = reinterpret_cast<const uint16_t*>(hashes_div);

    //Gather and compare
    const uint16_t* table = reinterpret_cast<const uint16_t*>(table_[blockIDs.first].block);

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

    const uint16_t* table2 = reinterpret_cast<const uint16_t*>(table_[blockIDs.second].block);

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

bool BlockedBloomFilter::WriteBloomFilter(FILE *fp) const {

    if (fwrite(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
    if (fwrite(&k_, sizeof(k_), 1, fp) != 1) return false;

    for (uint64_t i = 0; i != blocks_; ++i){

        if (fwrite(&(table_[i].block), sizeof(uint64_t), NB_ELEM_BLOCK, fp) != NB_ELEM_BLOCK) return false;
    }

    return true;
}

bool BlockedBloomFilter::ReadBloomFilter(FILE *fp) {

    clear();

    if (fread(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
    if (fread(&k_, sizeof(k_), 1, fp) != 1) return false;

    init_table();

    for (uint64_t i = 0; i != blocks_; ++i){

        if (fread(&(table_[i].block), sizeof(uint64_t), NB_ELEM_BLOCK, fp) != NB_ELEM_BLOCK) return false;
    }

    //uint64_t hashes_mask[4] __attribute__((aligned(32))) = {0, 0, 0, 0};

    for (int k = 0; k != k_; ++k) hashes_mask[k/4] = (hashes_mask[k/4] << 16) | 0xffff;

    //mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

    return true;
}

void BlockedBloomFilter::clear() {

    if (table_ != nullptr){

        delete[] table_;
        table_ = nullptr;
    }

    k_ = 0;
    blocks_ = 0;

    //mask_h = _mm256_setzero_si256();

    hashes_mask[0] = 0;
    hashes_mask[1] = 0;
    hashes_mask[2] = 0;
    hashes_mask[3] = 0;
}

void BlockedBloomFilter::init_table(){

    fast_div_ = libdivide::divider<uint64_t>(blocks_);
    table_ = new BBF_Block[blocks_];
}

bool BlockedBloomFilter::insert_par(const uint64_t kmer_hash, const uint64_t min_hash) {

    uint64_t hashes_div[4] __attribute__((aligned(32)));

    hashes_div[0] = kmer_hash;
    hashes_div[1] = hashes_div[0] + min_hash;
    hashes_div[2] = hashes_div[1] + min_hash;
    hashes_div[3] = hashes_div[2] + min_hash;

    const __m256i km_hashes = _mm256_and_si256(_mm256_set_epi64x(hashes_div[0], hashes_div[1], hashes_div[2], hashes_div[3]), mask_and_div);
    const __m256i h_shift = _mm256_and_si256(km_hashes, mask_and_mod);

    const __m256i hash_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
    const __m256i hash_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

    const __m256i mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

    const __m256i hash_gather = _mm256_and_si256(mask_h, _mm256_or_si256(hash_gather_lsb, hash_gather_msb));

    _mm256_stream_si256((__m256i*)hashes_div, _mm256_srli_epi16(km_hashes, 4));

    const uint16_t* hashes_div16 = reinterpret_cast<const uint16_t*>(hashes_div);

    const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

    uint64_t blockID = min_hash - (min_hash / fast_div_) * blocks_;
    uint64_t blockID_2 = min_hash_2 - (min_hash_2 / fast_div_) * blocks_;

    if (blockID_2 < blockID) std::swap(blockID_2, blockID);

    //Gather and compare
    uint16_t* table = reinterpret_cast<uint16_t*>(table_[blockID].block);

    __m256i table_gather;

    table_[blockID].lock();

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

    if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0){

        table_[blockID].unlock();
        return false;
    }

    uint16_t* table2 = reinterpret_cast<uint16_t*>(table_[blockID_2].block);

    if (blockID != blockID_2){

        table_[blockID_2].lock();

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

        if (_mm256_testz_si256(xor_m256i_2, xor_m256i_2) != 0){

            table_[blockID].unlock();
            table_[blockID_2].unlock();

            return false;
        }
    }

    uint16_t hashes_mod[16] __attribute__((aligned(32)));
    _mm256_stream_si256(reinterpret_cast<__m256i*>(hashes_mod), hash_gather);

    if (popcnt(table2, NB_ELEM_BLOCK * sizeof(uint64_t)) <= popcnt(table, NB_ELEM_BLOCK * sizeof(uint64_t))){

        table = table2;

        std::swap(blockID, blockID_2);
    }

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

    table_[blockID].unlock();
    if (blockID != blockID_2) table_[blockID_2].unlock();

    return true;
}

bool BlockedBloomFilter::insert_unpar(const uint64_t kmer_hash, const uint64_t min_hash) {

    uint64_t hashes_div[4] __attribute__((aligned(32)));

    hashes_div[0] = kmer_hash;
    hashes_div[1] = hashes_div[0] + min_hash;
    hashes_div[2] = hashes_div[1] + min_hash;
    hashes_div[3] = hashes_div[2] + min_hash;

    const __m256i km_hashes = _mm256_and_si256(_mm256_set_epi64x(hashes_div[0], hashes_div[1], hashes_div[2], hashes_div[3]), mask_and_div);
    const __m256i h_shift = _mm256_and_si256(km_hashes, mask_and_mod);

    const __m256i hash_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
    const __m256i hash_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

    const __m256i mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

    const __m256i hash_gather = _mm256_and_si256(mask_h, _mm256_or_si256(hash_gather_lsb, hash_gather_msb));

    _mm256_stream_si256((__m256i*)hashes_div, _mm256_srli_epi16(km_hashes, 4));

    const uint16_t* hashes_div16 = reinterpret_cast<const uint16_t*>(hashes_div);

    const uint64_t blockID = min_hash - (min_hash / fast_div_) * blocks_;

    //Gather and compare
    uint16_t* table = reinterpret_cast<uint16_t*>(table_[blockID].block);

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
    const uint64_t blockID_2 = min_hash_2 - (min_hash_2 / fast_div_) * blocks_;

    uint16_t* table2 = reinterpret_cast<uint16_t*>(table_[blockID_2].block);

    if (blockID != blockID_2){

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
    }

    uint16_t hashes_mod[16] __attribute__((aligned(32)));
    _mm256_stream_si256(reinterpret_cast<__m256i*>(hashes_mod), hash_gather);

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

    return true;
}

const __m256i BlockedBloomFilter::mask_and_div = _mm256_set1_epi16(MASK_BITS_BLOCK);
// All 4 LSB of each 16 bits word
const __m256i BlockedBloomFilter::mask_and_mod = _mm256_set1_epi16(0xf);
// All 1 LSB of each 32 bits word
const __m256i BlockedBloomFilter::one2shift_lsb = _mm256_set1_epi32(0x1);
// Set the 17th LSB bit of each 32 bits word
const __m256i BlockedBloomFilter::one2shift_msb = _mm256_set1_epi32(0x10000);
// All 16 LSB bits of each 32 bits word
const __m256i BlockedBloomFilter::mask_lsb = _mm256_set1_epi64x(0x0000ffff0000ffff);

#else

BlockedBloomFilter::BlockedBloomFilter() : table_(nullptr), blocks_(0), k_(0) {}

BlockedBloomFilter::BlockedBloomFilter(size_t nb_elem, size_t bits_per_elem) : table_(nullptr), blocks_(0), k_(0) {

    if ((nb_elem != 0) && (bits_per_elem != 0)){

        blocks_ = (bits_per_elem * nb_elem + MASK_BITS_BLOCK) / NB_BITS_BLOCK;
        k_ = (int) (bits_per_elem * log(2));

        if (fpp(bits_per_elem, k_) >= fpp(bits_per_elem, k_+1)) ++k_;

        init_table();
    }
}

BlockedBloomFilter::BlockedBloomFilter(const BlockedBloomFilter& o) : table_(nullptr), blocks_(o.blocks_), k_(o.k_), fast_div_(o.fast_div_) {

    if (blocks_ != 0){

        init_table();

        for (uint64_t i = 0; i != blocks_; ++i){

            memcpy(&(table_[i].block), &(o.table_[i].block), NB_ELEM_BLOCK * sizeof(uint64_t));
        }
    }
}

BlockedBloomFilter::BlockedBloomFilter(BlockedBloomFilter&& o) : table_(o.table_), blocks_(o.blocks_), k_(o.k_), fast_div_(o.fast_div_) {

    o.table_ = nullptr;

    o.clear();
}

BlockedBloomFilter::~BlockedBloomFilter() {

    clear();
}

BlockedBloomFilter& BlockedBloomFilter::operator=(const BlockedBloomFilter& o) {

    clear();

    blocks_ = o.blocks_;
    k_ = o.k_;
    fast_div_ = o.fast_div_;

    init_table();

    for (uint64_t i = 0; i != blocks_; ++i){

        memcpy(&(table_[i].block), &(o.table_[i].block), NB_ELEM_BLOCK * sizeof(uint64_t));
    }

    return *this;
}

BlockedBloomFilter& BlockedBloomFilter::operator=(BlockedBloomFilter&& o) {

    if (this != &o) {

        clear();

        table_ = o.table_;
        blocks_ = o.blocks_;
        k_ = o.k_;
        fast_div_ = o.fast_div_;

        o.table_ = nullptr;

        o.clear();
    }

    return *this;
}

int BlockedBloomFilter::contains(const uint64_t (&kmer_hash)[4], const uint64_t min_hash, bool (&pres)[4], const int limit) const{

    int cpt = 0;

    const uint64_t blockID = min_hash - (min_hash / fast_div_) * blocks_;

    for (size_t i, j = 0; (j != 4) && (cpt != limit); ++j){

        uint64_t kmer_hash_tmp = kmer_hash[j];

        for (i = 0; i != k_; ++i) {

            if ((table_[blockID].block[(kmer_hash_tmp & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_tmp & 0x3fULL))) == 0) break;
            kmer_hash_tmp = (kmer_hash_tmp * 49157) % (1610612741ULL);
        }

        cpt += (pres[j] = (i == k_));
    }

    if (cpt != limit){

        const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);
        const uint64_t blockID_2 = min_hash_2 - (min_hash_2 / fast_div_) * blocks_;

        for (size_t i, j = 0; (j != 4) && (cpt != limit); ++j){

            if (!pres[j]){

                uint64_t kmer_hash_tmp = kmer_hash[j];

                for (i = 0; i != k_; ++i) {

                    if ((table_[blockID_2].block[(kmer_hash_tmp & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_tmp & 0x3fULL))) == 0) break;
                    kmer_hash_tmp = (kmer_hash_tmp * 49157) % (1610612741ULL);
                }

                cpt += (pres[j] = (i == k_));
            }
        }
    }

    return cpt;
}

bool BlockedBloomFilter::contains(uint64_t kmer_hash, const uint64_t min_hash) const {

    int i = 0;

    uint64_t kmer_hash_2 = kmer_hash;

    const uint64_t blockID = min_hash - (min_hash / fast_div_) * blocks_;

    for (; i != k_; ++i) {

        if ((table_[blockID].block[(kmer_hash & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash & 0x3fULL))) == 0) break;
        kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
    }

    if (i != k_){

        const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);
        const uint64_t blockID_2 = min_hash_2 - (min_hash_2 / fast_div_) * blocks_;

        for (i = 0; i != k_; ++i) {

            if ((table_[blockID_2].block[(kmer_hash_2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_2 & 0x3fULL))) == 0) return false;
            kmer_hash_2 = (kmer_hash_2 * 49157) % (1610612741ULL);
        }

        return true;
    }

    return true;
}

size_t BlockedBloomFilter::contains_block(uint64_t kmer_hash, const uint64_t min_hash, const BBF_Blocks blockIDs) const {

    int i = 0;

    uint64_t kmer_hash_2 = kmer_hash;

    for (; i != k_; ++i) {

        if ((table_[blockIDs.first].block[(kmer_hash & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash & 0x3fULL))) == 0) break;
        kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
    }

    if (i != k_){

        for (i = 0; i != k_; ++i) {

            if ((table_[blockIDs.second].block[(kmer_hash_2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_2 & 0x3fULL))) == 0) return 0;
            kmer_hash_2 = (kmer_hash_2 * 49157) % (1610612741ULL);
        }

        return 2;
    }

    return 1;
}

bool BlockedBloomFilter::WriteBloomFilter(FILE *fp) const {

    if (fwrite(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
    if (fwrite(&k_, sizeof(k_), 1, fp) != 1) return false;

    for (uint64_t i = 0; i != blocks_; ++i){

        if (fwrite(&(table_[i].block), sizeof(uint64_t), NB_ELEM_BLOCK, fp) != NB_ELEM_BLOCK) return false;
    }

    return true;
}

bool BlockedBloomFilter::ReadBloomFilter(FILE *fp) {

    clear();

    if (fread(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
    if (fread(&k_, sizeof(k_), 1, fp) != 1) return false;

    init_table();

    for (uint64_t i = 0; i != blocks_; ++i){

        if (fread(&(table_[i].block), sizeof(uint64_t), NB_ELEM_BLOCK, fp) != NB_ELEM_BLOCK) return false;
    }

    return true;
}

void BlockedBloomFilter::clear() {

    if (table_ != nullptr){

        delete[] table_;
        table_ = nullptr;
    }

    blocks_ = 0;
    k_ = 0;
}

void BlockedBloomFilter::init_table(){

    fast_div_ = libdivide::divider<uint64_t>(blocks_);
    table_ = new BBF_Block[blocks_];
}

bool BlockedBloomFilter::insert_par(uint64_t kmer_hash, const uint64_t min_hash) {

    int i = 0;

    const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

    uint64_t kmer_hash_2 = kmer_hash;

    uint64_t blockID = min_hash - (min_hash / fast_div_) * blocks_;
    uint64_t blockID_2 = min_hash_2 - (min_hash_2 / fast_div_) * blocks_;

    if (blockID_2 < blockID) std::swap(blockID, blockID_2);

    table_[blockID].lock();

    for (; i != k_; ++i) {

        if ((table_[blockID].block[(kmer_hash & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash & 0x3fULL))) == 0) break;
        kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
    }

    if (i != k_){

        int j = 0;

        if (blockID_2 == blockID){

            j = i;
            kmer_hash_2 = kmer_hash;
        }
        else {

            table_[blockID_2].lock();

            for (; j != k_; ++j) {

                if ((table_[blockID_2].block[(kmer_hash_2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_2 & 0x3fULL))) == 0) break;
                kmer_hash_2 = (kmer_hash_2 * 49157) % (1610612741ULL);
            }
        }

        if (j != k_){

            if (popcnt(table_[blockID_2].block, NB_ELEM_BLOCK * sizeof(uint64_t)) <= popcnt(table_[blockID].block, NB_ELEM_BLOCK * sizeof(uint64_t))){

                i = j;
                kmer_hash = kmer_hash_2;

                std::swap(blockID, blockID_2);
            }

            for (; i != k_; ++i) {

                table_[blockID].block[(kmer_hash & MASK_BITS_BLOCK) >> 6] |= 1ULL << (kmer_hash & 0x3fULL);
                kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
            }

            table_[blockID].unlock();
            if (blockID_2 != blockID) table_[blockID_2].unlock();

            return true;
        }

        if (blockID_2 != blockID) table_[blockID_2].unlock();
    }

    table_[blockID].unlock();

    return false;
}

bool BlockedBloomFilter::insert_unpar(uint64_t kmer_hash, const uint64_t min_hash) {

    int i = 0;

    uint64_t kmer_hash_2 = kmer_hash;

    uint64_t blockID = (min_hash - (min_hash / fast_div_) * blocks_);

    for (; i != k_; ++i) {

        if ((table_[blockID].block[(kmer_hash & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash & 0x3fULL))) == 0) break;
        kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
    }

    if (i != k_){

        int j = 0;

        const uint64_t min_hash_2 = (min_hash * 49157) % (1610612741ULL);

        uint64_t blockID_2 = (min_hash_2 - (min_hash_2 / fast_div_) * blocks_);

        for (; j != k_; ++j) {

            if ((table_[blockID_2].block[(kmer_hash_2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmer_hash_2 & 0x3fULL))) == 0) break;
            kmer_hash_2 = (kmer_hash_2 * 49157) % (1610612741ULL);
        }

        if (j != k_){

            if (popcnt(table_[blockID_2].block, NB_ELEM_BLOCK * sizeof(uint64_t)) <= popcnt(table_[blockID].block, NB_ELEM_BLOCK * sizeof(uint64_t))){

                i = j;
                blockID = blockID_2;
                kmer_hash = kmer_hash_2;
            }

            for (; i != k_; ++i) {

                table_[blockID].block[(kmer_hash & MASK_BITS_BLOCK) >> 6] |= 1ULL << (kmer_hash & 0x3fULL);
                kmer_hash = (kmer_hash * 49157) % (1610612741ULL);
            }

            return true;
        }
    }

    return false;
}

#endif
