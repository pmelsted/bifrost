#ifndef COMPACTED_DBG_HPP
#define COMPACTED_DBG_HPP

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <vector>

#include <thread>
#include <atomic>

#include "BlockedBloomFilter.hpp"
#include "Common.hpp"
#include "CompressedSequence.hpp"
#include "Contig.hpp"
#include "fastq.hpp"
#include "Kmer.hpp"
#include "KmerHashTable.h"
#include "KmerIterator.hpp"
#include "minHashIterator.hpp"
#include "RepHash.hpp"
#include "TinyVector.hpp"
#include "UnitigMap.hpp"

#define MASK_CONTIG_ID (0xffffffff00000000)
#define MASK_CONTIG_TYPE (0x80000000)
#define MASK_CONTIG_POS (0x7fffffff)

#define RESERVED_ID (0xffffffff)

using namespace std;

template<bool is_const>
class unitigIterator;

class CompactedDBG {

    public:

        friend struct UnitigMap;

        template<bool is_const>
        friend class neighborIterator;

        template<bool is_const>
        friend class unitigIterator;

        CompactedDBG(int kmer_length = 31, int minimizer_length = 23);
        ~CompactedDBG();

        void clear();
        void empty();

        inline size_t size() const {

            return v_unitigs.size() + v_kmers.size() + h_kmers_ccov.size();
        }

        bool build(const vector<string>& fastx_filename_in, const size_t nb_unique_kmers, const size_t nb_non_unique_kmers,
                    const bool reference_mode = false, const size_t nb_threads = 1, const size_t nb_bits_unique_kmers_bf = 14,
                    const size_t nb_bits_non_unique_kmers_bf = 14, const string& inFilenameBBF = "", const string& outFilenameBBF = "",
                    const size_t read_chunksize = 10000, const size_t unitig_size = 1000000, const bool verbose = false);

        bool simplify(const bool delete_short_isolated_unitigs = true, const bool clip_short_tips = true, const bool verbose = false);

        bool write(const string output_filename, const bool verbose = false) const;

        UnitigMap find(const Kmer& km, bool extremities_only = false) const;

        bool add(const string& seq, const bool verbose);

        bool remove(const UnitigMap& um, const bool verbose);

        inline int getK() const {

            return k_;
        }

        typedef unitigIterator<true> iterator;
        typedef unitigIterator<false> const_iterator;

        iterator begin();
        iterator end();

        const_iterator begin() const;
        const_iterator end() const;

    private:

        bool join(const UnitigMap& um, const bool verbose);
        bool join(const bool verbose);

        bool filter(const vector<string>& fastx_filename_in, const size_t nb_unique_kmers, const size_t nb_non_unique_kmers,
                    const bool reference_mode = false, const size_t nb_threads = 1, const size_t nb_bits_unique_kmers_bf = 14,
                    const size_t nb_bits_non_unique_kmers_bf = 14, const size_t read_chunksize = 10000, const size_t unitig_size = 1000000,
                    const bool verbose = false);

        // -------- BBF must be non empty for these functions -----------
        bool construct(const vector<string>& fastx_filename_in, const size_t nb_threads = 1, const size_t read_chunksize = 1000,
                       const bool verbose = false);

        bool addUnitigSequenceBBF(Kmer km, const string& read, size_t pos, const string& seq, vector<Kmer>& l_ignored_km_tip);
        size_t findUnitigSequenceBBF(Kmer km, string& s, bool& selfLoop, bool& isIsolated, vector<Kmer>& l_ignored_km_tip);

        bool fwStepBBF(Kmer km, Kmer& end, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;
        bool bwStepBBF(Kmer km, Kmer& front, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;
        // ----------------------------------------------------------------

        UnitigMap findUnitig(const Kmer& km, const string& s, size_t pos) const;
        UnitigMap findUnitig(const Kmer& km, const string& s, size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h) const;

        bool addUnitig(const string& str_unitig, const size_t id_unitig);
        void deleteUnitig(const bool isShort, const bool isAbundant, const size_t id_unitig);
        void swapUnitigs(const bool isShort, const size_t id_a, const size_t id_b);

        bool splitUnitig(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs, size_t& v_unitigs_sz, size_t& v_kmers_sz,
                         const vector<pair<int,int>>& sp);

        UnitigMap find(const Kmer& km, const preAllocMinHashIterator<RepHash>& it_min_h) const;

        inline size_t find(const preAllocMinHashIterator<RepHash>& it_min_h) const {

            const int pos = it_min_h.getPosition();

            return (hmap_min_unitigs.find(Minimizer(&it_min_h.s[pos]).rep()) != hmap_min_unitigs.end() ? 0 : pos - it_min_h.p);
        }

        pair<size_t, size_t> splitAllUnitigs();
        size_t joinAllUnitigs(vector<Kmer>* v_joins = NULL);
        bool checkJoin(const Kmer& a, const UnitigMap& cm_a, Kmer& b);

        void writeGFA(string graphfilename) const;

        size_t unitigCount() const;

        void check_fp_tips(KmerHashTable<bool>& ignored_km_tips);

        size_t removeUnitigs(bool rmIsolated, bool clipTips, vector<Kmer>& v);

        void mapRead(const UnitigMap& cc);

        size_t cstrMatch(const char* a, const char* b) const;

        inline size_t stringMatch(const string& a, const string& b, size_t pos) {

            return distance(a.begin(), mismatch(a.begin(), a.end(), b.begin() + pos).first);
        }


        int k_;
        int g_;

        bool invalid;

        static const int tiny_vector_sz = 2;
        static const int min_abundance_lim = 15;
        static const int max_abundance_lim = 15;

        typedef KmerHashTable<CompressedCoverage> h_kmers_ccov_t;
        typedef MinimizerHashTable<tiny_vector<size_t,tiny_vector_sz>> hmap_min_unitigs_t;

        vector<Unitig*> v_unitigs;
        vector<pair<Kmer, CompressedCoverage>> v_kmers;

        hmap_min_unitigs_t hmap_min_unitigs;
        h_kmers_ccov_t h_kmers_ccov;

        BlockedBloomFilter bf;
};

template<bool is_const = true>
class unitigIterator : public std::iterator<std::input_iterator_tag, UnitigMap, int> {

    public:

        typedef typename std::conditional<is_const, const UnitigMap&, UnitigMap&>::type UnitigMap_ref_t;
        typedef typename std::conditional<is_const, const UnitigMap*, UnitigMap*>::type UnitigMap_ptr_t;

        unitigIterator() :  i_(0), v_unitigs_sz_(0), v_kmers_sz_(0), h_kmers_ccov_sz_(0), sz_(0), invalid_(true), dbg_(NULL) {}

        unitigIterator(const CompactedDBG* dbg) :
                        i_(0), v_unitigs_sz_(0), v_kmers_sz_(0), h_kmers_ccov_sz_(0), sz_(0), invalid_(true), dbg_(dbg),
                        it_h_kmers_ccov((dbg == NULL) || dbg->invalid ? KmerHashTable<CompressedCoverage>::const_iterator() : dbg->h_kmers_ccov.begin()){

            if ((dbg_ != NULL) && !dbg_->invalid && (dbg_->size() != 0)){

                invalid_ = false;

                v_unitigs_sz_ = dbg_->v_unitigs.size();
                v_kmers_sz_ = dbg_->v_kmers.size();
                h_kmers_ccov_sz_ = dbg_->h_kmers_ccov.size();

                sz_ = v_unitigs_sz_ + v_kmers_sz_ + h_kmers_ccov_sz_;
            }
        }

        unitigIterator(const unitigIterator& o) :   i_(o.i_), v_unitigs_sz_(o.v_unitigs_sz_), v_kmers_sz_(o.v_kmers_sz_),
                                                    h_kmers_ccov_sz_(o.h_kmers_ccov_sz_), sz_(o.sz_), invalid_(o.invalid_),
                                                    um_(o.um_), dbg_(o.dbg_) {}

        unitigIterator& operator++() {

            if (invalid_) return *this;

            if ((dbg_ == NULL) || dbg_->invalid || (i_ >= sz_)){

                invalid_ = true;
                return *this;
            }

            if (i_ < v_unitigs_sz_) um_ = UnitigMap(i_, 0, 1, dbg_->v_unitigs[i_]->seq.size(), false, false, true, dbg_);
            else if (i_ < (v_unitigs_sz_ + v_kmers_sz_)){

                um_ = UnitigMap(i_ - v_unitigs_sz_, 0, 1, dbg_->getK(), true, false, true, dbg_);
            }
            else {

                um_ = UnitigMap(it_h_kmers_ccov.getHash(), 0, 1, dbg_->getK(), false, true, true, dbg_);

                it_h_kmers_ccov++;
            }

            i_++;

            return *this;
        }

        unitigIterator operator++(int) {

            unitigIterator tmp(*this);
            operator++();

            return tmp;
        }

        bool operator==(const unitigIterator& o) {

            if (invalid_ || o.invalid_) return invalid_ && o.invalid_;
            return  (i_ == o.i_) && (v_unitigs_sz_ == o.v_unitigs_sz_) && (v_kmers_sz_ == o.v_kmers_sz_) &&
                    (h_kmers_ccov_sz_ == o.h_kmers_ccov_sz_) && (sz_ == o.sz_) && (it_h_kmers_ccov == o.it_h_kmers_ccov) &&
                    (dbg_ == o.dbg_) && (um_ == o.um_);
        }

        bool operator!=(const unitigIterator& o) { return !operator==(o); }

        UnitigMap_ref_t operator*() { return um_; }
        UnitigMap_ptr_t operator->() { return &um_; }

    protected:

        size_t i_;

        size_t v_unitigs_sz_;
        size_t v_kmers_sz_;
        size_t h_kmers_ccov_sz_;
        size_t sz_;

        bool invalid_;

        KmerHashTable<CompressedCoverage>::const_iterator it_h_kmers_ccov;

        UnitigMap um_;

        const CompactedDBG* dbg_;
};

template<typename T>
class superCompactedDBG {

    public:

        superCompactedDBG(int kmer_length = 31, int minimizer_length = 23) : data_unitigs(NULL), data_kmers(NULL) {

            dbg = CompactedDBG(kmer_length, minimizer_length);
            has_data = typeid(T) != typeid(void);

            if (has_data){

                data_unitigs = new T[1];
                data_kmers = new T[1];
            }
        }

        ~superCompactedDBG(){

            if (data_unitigs != NULL) delete[] data_unitigs;
            if (data_kmers != NULL) delete[] data_kmers;
        }

        inline void clear(){ dbg.clear(); }
        inline void empty(){ dbg.empty(); }

        inline size_t size() const { return dbg.size(); }

        inline bool build(const vector<string>& fastx_filename_in, const size_t nb_unique_kmers, const size_t nb_non_unique_kmers,
                    const bool reference_mode = false, const size_t nb_threads = 1, const size_t nb_bits_unique_kmers_bf = 14,
                    const size_t nb_bits_non_unique_kmers_bf = 14, const string& inFilenameBBF = "", const string& outFilenameBBF = "",
                    const size_t read_chunksize = 10000, const size_t unitig_size = 1000000, const bool verbose = false){

            return dbg.build(fastx_filename_in, nb_unique_kmers, nb_non_unique_kmers, reference_mode, nb_threads,
                             nb_bits_unique_kmers_bf, nb_bits_non_unique_kmers_bf, inFilenameBBF, outFilenameBBF,
                             read_chunksize, unitig_size, verbose);
        }

        inline bool simplify(const bool delete_short_isolated_unitigs = true, const bool clip_short_tips = true, const bool verbose = false){

            return dbg.simplify(delete_short_isolated_unitigs, clip_short_tips, verbose);
        }

        inline bool write(const string output_filename, const bool verbose = false) const { return dbg.write(output_filename, verbose); }

        inline UnitigMap find(const Kmer& km, bool extremities_only = false) const { return dbg.find(km, extremities_only); }

        inline bool add(const string& seq, const bool verbose) { return dbg.add(seq, verbose); }

        inline bool remove(const UnitigMap& um, const bool verbose){ return dbg.remove(um, verbose); }

        inline int getK() const { return dbg.getK(); }

        typedef unitigIterator<true> iterator;
        typedef unitigIterator<false> const_iterator;

        inline iterator begin() { return dbg.begin(); }
        inline iterator end() { return dbg.end(); }

        inline const_iterator begin() const { return dbg.begin(); }
        inline const_iterator end() const { return dbg.end(); }

    private:

        CompactedDBG dbg;

        bool has_data;

        T* data_unitigs;
        T* data_kmers;
};

#endif // MINHASHITERATOR_H
