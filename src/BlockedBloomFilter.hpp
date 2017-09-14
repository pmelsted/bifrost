#ifndef BFG_BLOCKEDBLOOMFILTER_HPP
#define BFG_BLOCKEDBLOOMFILTER_HPP

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "hash.hpp"
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

static const uint64_t mask[8] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};

/* Short description:
 *  - Extended BloomFilter which hashes into 64-bit blocks
 *    that can be accessed very fast from the CPU cache
 * */
/*class BlockedBloomFilter {

    private:

        uint64_t* table_; //Bit array

        uint64_t size_table_; //Size of bit array (in bits)
        uint64_t size_blocks_; //Size of each block (in bits)
        uint64_t blocks_; //Nb blocks
        uint64_t nb_table_t_block; //Number of elements of type table_t in each block
        uint64_t k_; //Nb hash functions

        libdivide::divider<uint64_t> fast_div_; // fast division




        uint64_t* counts;
        std::vector<Minimizer>* v_min;

        typedef MinimizerHashTable<std::vector<std::pair<Kmer, size_t>>> hmap_min_t;
        hmap_min_t hmap_min;

    public:

        BlockedBloomFilter() : table_(NULL), size_table_(0), size_blocks_(0), blocks_(0), k_(0), nb_table_t_block(0), fast_div_() {}

        BlockedBloomFilter(size_t nb_elem, size_t bits_per_elem, size_t bits_per_block = 512) : table_(NULL), size_table_(0), size_blocks_(0),
                                                                                                blocks_(0), k_(0), nb_table_t_block(0), fast_div_() {

            size_blocks_ = ((bits_per_block + 63) / 64) * 64;
            size_table_ = ((bits_per_elem * nb_elem + size_blocks_ - 1) / size_blocks_) * size_blocks_;
            blocks_ = size_table_ / size_blocks_;
            nb_table_t_block = size_blocks_ / 64;
            k_ = (uint64_t) (bits_per_elem * log(2));

            if (fpp(bits_per_elem, k_) >= fpp(bits_per_elem, k_+1)) k_++;

            init_table();

            init_counts();
        }

        ~BlockedBloomFilter() {

            clear();

            clear_counts();
        }

        inline uint64_t* getBlock(const uint64_t min_hash) const{

            return table_ + nb_table_t_block * (min_hash - (min_hash / fast_div_) * blocks_);
        }

        inline bool contains(uint64_t kmer_hash, const uint64_t min_hash) const { return search(kmer_hash, min_hash); }

        inline bool contains(uint64_t kmer_hash, const uint64_t* block_ptr) const { return search(kmer_hash, block_ptr); }

        bool search(uint64_t kmer_hash, const uint64_t min_hash) const {

            uint64_t block = min_hash - (min_hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+ nb_table_t_block * block, 0, 1);

            for (uint64_t bit, maskcheck, loc, i = 0; i < k_; i++) {

                bit = kmer_hash & (size_blocks_ - 1); // position of bit to look for in block;
                maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                loc = nb_table_t_block * block + (bit >> 6); // position of elem uint64_t containing bit to check

                if ((table_[loc] &  maskcheck) == 0) return false;

                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);
            }

            return true;
        }

        bool search(uint64_t kmer_hash, const uint64_t* block_ptr) const {

            __builtin_prefetch(block_ptr, 0, 1);

            for (uint64_t bit, maskcheck, i = 0; i < k_; i++) {

                bit = kmer_hash & (size_blocks_ - 1); // position of bit to look for in block;
                maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit

                if ((block_ptr[bit >> 6] & maskcheck) == 0) return false;

                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);
            }

            return true;
        }

        void insert(uint64_t kmer_hash, uint64_t min_hash) {

            uint64_t block = min_hash - (min_hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+ nb_table_t_block * block, 1, 1);

            for (uint64_t val, bit, maskcheck, loc, i = 0; i < k_; i++) {

                bit = kmer_hash & (size_blocks_ - 1); // position of bit to look for in block;
                maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                loc = nb_table_t_block * block + (bit >> 6); // position of elem uint64_t containing bit to check
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL); //Next hash value for k-mer

                if ((table_[loc] &  maskcheck) == 0) __sync_fetch_and_or(table_ + loc, maskcheck);
            }

            return;
        }

        bool search_and_insert(uint64_t kmer_hash, uint64_t min_hash) {

            bool r = false;

            uint64_t block = min_hash - (min_hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+ nb_table_t_block * block, 1, 1);

            for (uint64_t bit, maskcheck, loc, val, i = 0; i < k_; i++) {

                bit = kmer_hash & (size_blocks_ - 1); // position of bit to look for in block;
                maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                loc = nb_table_t_block * block + (bit >> 6); // position of elem uint64_t containing bit to check
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL); //Next hash value for k-mer

                if ((table_[loc] & maskcheck) == 0) {

                    val = __sync_fetch_and_or(table_ + loc, maskcheck);

                    if (!r && ((val & maskcheck) == 0)) r = true;
                }
            }

            return r;
        }

        bool WriteBloomFilter(FILE *fp) {

            if (fwrite(&nb_table_t_block, sizeof(nb_table_t_block), 1, fp) != 1) return false;
            if (fwrite(&size_table_, sizeof(size_table_), 1, fp) != 1) return false;
            if (fwrite(&size_blocks_, sizeof(size_blocks_), 1, fp) != 1) return false;
            if (fwrite(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
            if (fwrite(&k_, sizeof(k_), 1, fp) != 1) return false;

            if (fwrite(table_, sizeof(uint64_t), nb_table_t_block * blocks_, fp) != (nb_table_t_block * blocks_)) return false;

            return true;
        }

        bool ReadBloomFilter(FILE *fp) {

            clear();

            if (fread(&nb_table_t_block, sizeof(nb_table_t_block), 1, fp) != 1) return false;
            if (fread(&size_table_, sizeof(size_table_), 1, fp) != 1) return false;
            if (fread(&size_blocks_, sizeof(size_blocks_), 1, fp) != 1) return false;
            if (fread(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
            if (fread(&k_, sizeof(k_), 1, fp) != 1) return false;

            init_table();

            if (fread(table_, sizeof(uint64_t), nb_table_t_block * blocks_, fp) != (nb_table_t_block * blocks_)) return false;

            return true;
        }

        void clear() {

            if (table_ != NULL){

                free(table_);
                table_ = NULL;
            }

            nb_table_t_block = 0;
            size_table_ = 0;
            size_blocks_ = 0;
            blocks_ = 0;
            k_ = 0;
        }

        void clear_counts(){

            if (counts != NULL){

                delete[] counts;
                counts = NULL;
            }

            if (v_min != NULL){

                delete[] v_min;
                v_min = NULL;
            }
        }

        inline void inc_count_block(uint64_t min_hash, Minimizer minz, Kmer km){

            uint64_t block = min_hash - (min_hash / fast_div_) * blocks_; //Block ID

            Minimizer minz_rep = minz.rep();
            Kmer km_rep = km.rep();

            std::vector<Minimizer>& v_min_b = v_min[block]; //If minz not present for this block, add it

            if (std::find(v_min_b.begin(), v_min_b.end(), minz_rep) == v_min_b.end()) v_min_b.push_back(minz_rep);

            std::pair<hmap_min_t::iterator, bool> it = hmap_min.insert(make_pair(minz_rep, std::vector<std::pair<Kmer, size_t>>()));
            std::vector<std::pair<Kmer, size_t>>& v = it.first->second;

            bool add = true;

            if (v.size() != 0){

                for (auto& p : v){

                    if (p.first == km_rep){

                        p.second++;
                        add = false;
                        break;
                    }
                }
            }

            if (add){

                v.push_back(std::make_pair(km_rep, (size_t)1));
                counts[block]++;
            }
        }

        void print_count_blocks(std::string filename, std::string suffix){

            std::cerr << "hmap_min.size(): " << hmap_min.size() << ", blocks_: " << blocks_ << std::endl;

            std::string real_filename = filename + "_kmer_counts_" + suffix + ".csv";

            std::ofstream ofile;
            std::ostream sfile(0);

            ofile.open(real_filename.c_str());
            assert(!ofile.fail());
            sfile.rdbuf(ofile.rdbuf());

            uint64_t min_count_kmer = 0xffffffffffffffff, max_count_kmer = 0, tot_count_kmers = 0;
            uint64_t min_count_minz = 0xffffffffffffffff, max_count_minz = 0, tot_count_minz = 0;
            uint64_t min_count_kmer_per_minz = 0xffffffffffffffff, max_count_kmer_per_minz = 0;

            for (size_t i = 0; i < blocks_; i++){

                if (counts[i] < min_count_kmer) min_count_kmer = counts[i];
                if (counts[i] > max_count_kmer) max_count_kmer = counts[i];

                if (v_min[i].size() < min_count_minz) min_count_minz = v_min[i].size();
                if (v_min[i].size() > max_count_minz) max_count_minz = v_min[i].size();

                tot_count_kmers += counts[i];
                tot_count_minz += v_min[i].size();
            }

            for (hmap_min_t::iterator it = hmap_min.begin(); it != hmap_min.end(); it++){

                //std::cerr << it->first.toString() << std::endl;

                std::vector<std::pair<Kmer, size_t>>& v = it->second;

                if (v.size() < min_count_kmer_per_minz) min_count_kmer_per_minz = v.size();
                if (v.size() > max_count_kmer_per_minz) max_count_kmer_per_minz = v.size();
            }

            uint64_t nb_blocks_with_kmer_count[max_count_kmer - min_count_kmer + 1] = {0};
            uint64_t nb_blocks_with_minz_count[max_count_minz - min_count_minz + 1] = {0};
            uint64_t nb_minz_with_kmer_count[max_count_kmer_per_minz - min_count_kmer_per_minz + 1] = {0};

            for (size_t i = 0; i < blocks_; i++){

                nb_blocks_with_kmer_count[counts[i] - min_count_kmer]++;
                nb_blocks_with_minz_count[v_min[i].size() - min_count_minz]++;
            }

            RepHash rep_h;
            uint64_t avg_hash_value = 0, cnt_hash = 0, max_hash_v = 0, min_hash_v = 0xffffffffffffffff;

            rep_h.setK(Minimizer::g);

            for (hmap_min_t::iterator it = hmap_min.begin(); it != hmap_min.end(); it++){

                std::string what = it->first.toString();

                rep_h.init(what.c_str());
                avg_hash_value += rep_h.hash();
                cnt_hash++;

                if (rep_h.hash() < min_hash_v) min_hash_v = rep_h.hash();
                if (rep_h.hash() > max_hash_v) max_hash_v = rep_h.hash();

                nb_minz_with_kmer_count[it->second.size() - min_count_kmer_per_minz]++;
            }

            std::cerr << "Average hash value = " << (avg_hash_value/cnt_hash) << std::endl;
            std::cerr << "Max hash value = " << max_hash_v << std::endl;
            std::cerr << "Min hash value = " << min_hash_v << std::endl;

            size_t avg_nb_kmer_per_block = ceil(((double)tot_count_kmers) / ((double)blocks_));

            uint64_t nb_blocks_overloaded = 0;

            for (size_t i = avg_nb_kmer_per_block + 1 - min_count_kmer; i <= max_count_kmer - min_count_kmer; i++)
                nb_blocks_overloaded += nb_blocks_with_kmer_count[i];

            std::cerr << std::endl << "Total number of unique k-mer inserted is " << tot_count_kmers << std::endl;
            std::cerr << "Total nb blocks is " << blocks_ << std::endl;
            //cerr << "Average number of k-mer per block should be " << avg_nb_kmer_per_block << endl;
            //cerr << "Number of blocks overloaded compared to average is " << nb_blocks_overloaded <<
            //" (" << ((((double)nb_blocks_overloaded) / ((double)blocks_)) * 100) << "%)" << endl << endl;

            for (size_t i = 0; i < min_count_kmer; i++) sfile << i << "\t0" << std::endl;

            for (size_t i = 0; i <= max_count_kmer - min_count_kmer; i++){
                //cerr << "Nb blocks with " << i << " k-mers is " << nb_blocks_with_kmer_count[i] << "("
                //<< ((((double)nb_blocks_with_kmer_count[i])/((double)blocks_)) * 100) << "%)" << endl;
                sfile << (i + min_count_kmer) << "\t" << nb_blocks_with_kmer_count[i] << std::endl;
            }

            ofile.close();

            real_filename = filename + "_minimizer_counts_" + suffix + ".csv";

            ofile.open(real_filename.c_str());
            assert(!ofile.fail());
            sfile.rdbuf(ofile.rdbuf());

            for (size_t i = 0; i < min_count_minz; i++) sfile << i << "\t0" << std::endl;

            for (size_t i = 0; i <= max_count_minz - min_count_minz; i++){
                //cerr << "Nb blocks with " << i << " minimizers is " << nb_blocks_with_minz_count[i] << "("
                //<< ((((double)nb_blocks_with_minz_count[i])/((double)blocks_)) * 100) << "%)" << endl;
                sfile << (i + min_count_minz) << "\t" << nb_blocks_with_minz_count[i] << std::endl;
            }

            ofile.close();

            real_filename = filename + "_minimizer2kmer_counts_" + suffix + ".csv";

            ofile.open(real_filename.c_str());
            assert(!ofile.fail());
            sfile.rdbuf(ofile.rdbuf());

            for (size_t i = 0; i < min_count_kmer_per_minz; i++) sfile << i << "\t0" << std::endl;

            for (size_t i = 0; i <= max_count_kmer_per_minz - min_count_kmer_per_minz; i++){
                sfile << (i + min_count_kmer_per_minz) << "\t" << nb_minz_with_kmer_count[i] << std::endl;
            }

            ofile.close();

            //for (size_t i = 1; i < blocks_; i++) counts[i] += counts[i-1];

            real_filename = filename + "_distrib_minz" + suffix + ".csv";

            ofile.open(real_filename.c_str());
            assert(!ofile.fail());
            sfile.rdbuf(ofile.rdbuf());

            for (size_t i = 0; i < blocks_; i++){
                sfile << i << "\t" << counts[i] << std::endl;
            }

            ofile.close();
        }

    private:

        void init_table(){

            fast_div_ = libdivide::divider<uint64_t>(blocks_);

            posix_memalign((void**)&table_, 64, nb_table_t_block * blocks_* sizeof(table_[0]));
            memset(table_, 0, nb_table_t_block * blocks_ * sizeof(table_[0]));
        }

        inline double fpp(size_t bits, size_t k) const {

            return pow(1-exp(-((double)k)/((double)bits)),(double)k);
        }

        void init_counts(){

            counts = new uint64_t[blocks_];
            memset(counts, 0, blocks_ * sizeof(counts[0]));

            v_min = new std::vector<Minimizer>[blocks_];
        }
};*/

/*class BlockedBloomFilter {

    private:

        uint64_t* table_; //Bit array

        uint64_t size_table_; //Size of bit array (in bits)
        uint64_t blocks_; //Nb blocks
        uint64_t k_; //Nb hash functions

        libdivide::divider<uint64_t> fast_div_; // fast division

    public:

        BlockedBloomFilter() : table_(NULL), size_table_(0), blocks_(0), k_(0), fast_div_() {}

        BlockedBloomFilter(size_t nb_elem, size_t bits_per_elem) : table_(NULL), size_table_(0), blocks_(0), k_(0), fast_div_() {

            size_table_ = ((bits_per_elem * nb_elem + MASK_BITS_BLOCK) / NB_BITS_BLOCK) * NB_BITS_BLOCK;
            blocks_ = size_table_ / NB_BITS_BLOCK;

            init_table();

            k_ = (uint64_t) (bits_per_elem * log(2));
            if (fpp(bits_per_elem, k_) >= fpp(bits_per_elem, k_+1)) k_++;
        }

        ~BlockedBloomFilter() {

            clear();
        }

        inline uint64_t* getBlock(const uint64_t min_hash) const{

            return table_ + NB_ELEM_BLOCK * (min_hash - (min_hash / fast_div_) * blocks_);
        }

        inline bool contains(uint64_t kmer_hash, const uint64_t min_hash) const { return search(kmer_hash, min_hash); }

        inline bool contains(uint64_t kmer_hash, const uint64_t* block_ptr) const { return search(kmer_hash, block_ptr); }

        bool search(uint64_t kmer_hash, const uint64_t min_hash) const {

            uint64_t block = min_hash - (min_hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+ NB_ELEM_BLOCK * block, 0, 1);

            for (uint64_t bit, maskcheck, loc, i = 0; i < k_; i++) {

                bit = kmer_hash & MASK_BITS_BLOCK; // position of bit to look for in block;
                maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                loc = NB_ELEM_BLOCK * block + (bit >> 6); // position of elem uint64_t containing bit to check

                if ((table_[loc] &  maskcheck) == 0) return false;

                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);
            }

            return true;
        }

        bool search(uint64_t kmer_hash, const uint64_t* block_ptr) const {

            __builtin_prefetch(block_ptr, 0, 1);

            for (uint64_t bit, maskcheck, i = 0; i < k_; i++) {

                bit = kmer_hash & MASK_BITS_BLOCK; // position of bit to look for in block;
                maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit

                if ((block_ptr[bit >> 6] & maskcheck) == 0) return false;

                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);
            }

            return true;
        }

        void insert(uint64_t kmer_hash, uint64_t min_hash) {

            uint64_t block = min_hash - (min_hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+ NB_ELEM_BLOCK * block, 1, 1);

            for (uint64_t val, bit, maskcheck, loc, i = 0; i < k_; i++) {

                bit = kmer_hash & MASK_BITS_BLOCK; // position of bit to look for in block;
                maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                loc = NB_ELEM_BLOCK * block + (bit >> 6); // position of elem uint64_t containing bit to check
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL); //Next hash value for k-mer

                if ((table_[loc] &  maskcheck) == 0) __sync_fetch_and_or(table_ + loc, maskcheck);
            }

            return;
        }

        bool search_and_insert(uint64_t kmer_hash, uint64_t min_hash) {

            bool r = false;

            uint64_t block = min_hash - (min_hash / fast_div_) * blocks_;

            __builtin_prefetch(table_+ NB_ELEM_BLOCK * block, 1, 1);

            for (uint64_t bit, maskcheck, loc, val, i = 0; i < k_; i++) {

                bit = kmer_hash & MASK_BITS_BLOCK; // position of bit to look for in block;
                maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                loc = NB_ELEM_BLOCK * block + (bit >> 6); // position of elem uint64_t containing bit to check
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL); //Next hash value for k-mer

                if ((table_[loc] & maskcheck) == 0) {

                    val = __sync_fetch_and_or(table_ + loc, maskcheck);

                    if (!r && ((val & maskcheck) == 0)) r = true;
                }
            }

            return r;
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

    private:

        void init_table(){

            fast_div_ = libdivide::divider<uint64_t>(blocks_);

            posix_memalign((void**)&table_, 64, NB_ELEM_BLOCK * blocks_* sizeof(table_[0]));
            memset(table_, 0, NB_ELEM_BLOCK * blocks_ * sizeof(table_[0]));
        }

        inline double fpp(size_t bits, size_t k) const {

            return pow(1-exp(-((double)k)/((double)bits)),(double)k);
        }
};*/

class BlockedBloomFilter {

    private:

        uint64_t* table_; //Bit array

        uint64_t size_table_; //Size of bit array (in bits)
        uint64_t blocks_; //Nb blocks
        uint64_t k_; //Nb hash functions

        libdivide::divider<uint64_t> fast_div_; // fast division

    public:

        BlockedBloomFilter() : table_(NULL), size_table_(0), blocks_(0), k_(0), fast_div_() {}

        BlockedBloomFilter(size_t nb_elem, size_t bits_per_elem) : table_(NULL), size_table_(0), blocks_(0), k_(0), fast_div_() {

            size_table_ = ((bits_per_elem * nb_elem + MASK_BITS_BLOCK) / NB_BITS_BLOCK) * NB_BITS_BLOCK;
            blocks_ = size_table_ / NB_BITS_BLOCK;

            init_table();

            k_ = (uint64_t) (bits_per_elem * log(2));
            if (fpp(bits_per_elem, k_) >= fpp(bits_per_elem, k_+1)) k_++;
        }

        ~BlockedBloomFilter() {

            clear();
        }

        inline std::pair<uint64_t*,uint64_t*> getBlock(uint64_t min_hash) const{

            uint64_t min_hash_2 = (min_hash * 48271) % (2147483647ULL);

            min_hash -= (min_hash / fast_div_) * blocks_;
            min_hash_2 -= (min_hash_2 / fast_div_) * blocks_;

            return std::make_pair(table_ + NB_ELEM_BLOCK * min_hash, table_ + NB_ELEM_BLOCK * min_hash_2);
        }

        bool contains(uint64_t kmer_hash, const uint64_t min_hash) const {

            bool r = true;

            uint64_t kmer_hash_2 = kmer_hash;
            uint64_t block = (min_hash - (min_hash / fast_div_) * blocks_) * NB_ELEM_BLOCK;

            __builtin_prefetch(table_ + block, 0, 1);

            for (uint64_t bit, maskcheck, loc, i = 0; (i < k_) && r; i++) {

                bit = kmer_hash & MASK_BITS_BLOCK; // position of bit to look for in block;
                maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                loc = block + (bit >> 6); // position of elem uint64_t containing bit to check
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);

                if ((table_[loc] & maskcheck) == 0) r = false;
            }

            if (!r){

                uint64_t min_hash_2 = (min_hash * 48271) % (2147483647ULL);
                block = (min_hash_2 - (min_hash_2 / fast_div_) * blocks_) * NB_ELEM_BLOCK;

                r = true;

                __builtin_prefetch(table_ + block, 0, 1);

                for (uint64_t bit, maskcheck, loc, i = 0; (i < k_) && r; i++) {

                    bit = kmer_hash_2 & MASK_BITS_BLOCK; // position of bit to look for in block;
                    maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                    loc = block + (bit >> 6); // position of elem uint64_t containing bit to check
                    kmer_hash_2 = (kmer_hash_2 * 48271) % (2147483647ULL);

                    if ((table_[loc] & maskcheck) == 0) r = false;
                }
            }

            return r;
        }

        inline bool contains(uint64_t kmer_hash, const std::pair<uint64_t*, uint64_t*>& block_ptr) const {

            return (contains_block(kmer_hash, block_ptr) != 0);
        }

        size_t contains_block(uint64_t kmer_hash, const std::pair<uint64_t*, uint64_t*>& block_ptr) const {

            size_t r = 1;

            uint64_t kmer_hash_2 = kmer_hash;

            __builtin_prefetch(block_ptr.first, 0, 1);

            for (uint64_t bit, maskcheck, i = 0; (i < k_) && (r != 0); i++) {

                bit = kmer_hash & MASK_BITS_BLOCK; // position of bit to look for in block;
                maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);

                if ((block_ptr.first[bit >> 6] & maskcheck) == 0) r = 0;
            }

            if (!r){

                r = 2;

                __builtin_prefetch(block_ptr.second, 0, 1);

                for (uint64_t bit, maskcheck, loc, i = 0; (i < k_) && (r != 0); i++) {

                    bit = kmer_hash_2 & MASK_BITS_BLOCK; // position of bit to look for in block;
                    maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                    kmer_hash_2 = (kmer_hash_2 * 48271) % (2147483647ULL);

                    if ((block_ptr.second[bit >> 6] & maskcheck) == 0) r = 0;
                }
            }

            return r;
        }

        bool search_and_insert(uint64_t kmer_hash, const uint64_t min_hash, bool multi_threaded = false) {

            bool r = true;

            uint64_t kmer_hash_2 = kmer_hash;
            uint64_t block = (min_hash - (min_hash / fast_div_) * blocks_) * NB_ELEM_BLOCK;

            uint64_t i, j;
            uint64_t bit, maskcheck, loc;

            __builtin_prefetch(table_ + block, 0, 1);

            for (i = 0; i < k_; i++) {

                bit = kmer_hash & MASK_BITS_BLOCK; // position of bit to look for in block;
                maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                loc = block + (bit >> 6); // position of elem uint64_t containing bit to check

                if ((table_[loc] & maskcheck) == 0){

                    r = false;
                    break;
                }

                kmer_hash = (kmer_hash * 48271) % (2147483647ULL);
            }

            if (!r){

                r = true;

                uint64_t min_hash_2 = (min_hash * 48271) % (2147483647ULL);
                uint64_t block_2 = (min_hash_2 - (min_hash_2 / fast_div_) * blocks_) * NB_ELEM_BLOCK;

                __builtin_prefetch(table_ + block_2, 0, 1);

                for (j = 0; j < k_; j++) {

                    bit = kmer_hash_2 & MASK_BITS_BLOCK; // position of bit to look for in block;
                    maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                    loc = block_2 + (bit >> 6); // position of elem uint64_t containing bit to check

                    if ((table_[loc] & maskcheck) == 0){

                        r = false;
                        break;
                    }

                    kmer_hash_2 = (kmer_hash_2 * 48271) % (2147483647ULL);
                }

                if (!r){

                    if (!multi_threaded){

                        if (popcnt(table_ + block_2, NB_ELEM_BLOCK * sizeof(uint64_t)) < popcnt(table_ + block, NB_ELEM_BLOCK * sizeof(uint64_t))){

                            i = j;
                            block = block_2;
                            kmer_hash = kmer_hash_2;
                        }

                        __builtin_prefetch(table_ + block, 1, 1);

                        for (; i < k_; i++) {

                            bit = kmer_hash & MASK_BITS_BLOCK; // position of bit to look for in block;
                            maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                            loc = block + (bit >> 6); // position of elem uint64_t containing bit to check
                            kmer_hash = (kmer_hash * 48271) % (2147483647ULL); //Next hash value for k-mer

                            __sync_fetch_and_or(table_ + loc, maskcheck);
                        }
                    }
                    else {

                        if (popcnt(table_ + block_2, NB_ELEM_BLOCK * sizeof(uint64_t)) < popcnt(table_ + block, NB_ELEM_BLOCK * sizeof(uint64_t))){

                            uint64_t tmp = i;

                            i = j;
                            j = tmp;

                            tmp = block;
                            block = block_2;
                            block_2 = tmp;

                            tmp = kmer_hash;
                            kmer_hash = kmer_hash_2;
                            kmer_hash_2 = kmer_hash;
                        }

                        __builtin_prefetch(table_ + block, 1, 1);

                        for (; i < k_; i++) {

                            bit = kmer_hash & MASK_BITS_BLOCK; // position of bit to look for in block;
                            maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                            loc = block + (bit >> 6); // position of elem uint64_t containing bit to check
                            kmer_hash = (kmer_hash * 48271) % (2147483647ULL); //Next hash value for k-mer

                            __sync_fetch_and_or(table_ + loc, maskcheck);
                        }

                        __builtin_prefetch(table_ + block_2, 0, 1);

                        for (; j < k_; j++) {

                            bit = kmer_hash_2 & MASK_BITS_BLOCK; // position of bit to look for in block;
                            maskcheck = 1ULL << (bit & 0x3fULL); // mask to get this bit
                            loc = block_2 + (bit >> 6); // position of elem uint64_t containing bit to check

                            if ((table_[loc] & maskcheck) == 0) break;

                            kmer_hash_2 = (kmer_hash_2 * 48271) % (2147483647ULL);
                        }

                        if (j == k_) r = true;
                    }
                }
            }

            return !r;
        }

        inline void insert(uint64_t kmer_hash, const uint64_t min_hash){

            search_and_insert(kmer_hash, min_hash, false);
            return;
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

        inline uint64_t getNbBlocks() const { return blocks_; }

        inline uint64_t* getTable_ptr() const { return table_; }

    private:

        void init_table(){

            fast_div_ = libdivide::divider<uint64_t>(blocks_);

            posix_memalign((void**)&table_, 64, NB_ELEM_BLOCK * blocks_* sizeof(table_[0]));
            memset(table_, 0, NB_ELEM_BLOCK * blocks_ * sizeof(table_[0]));
        }

        inline double fpp(size_t bits, size_t k) const {

            return pow(1-exp(-((double)k)/((double)bits)),(double)k);
        }
};

#endif // BFG_BLOCKEDBLOOMFILTER_HPP
