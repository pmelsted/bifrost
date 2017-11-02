#ifndef SUPER_COMPACTED_DBG_HPP
#define SUPER_COMPACTED_DBG_HPP

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
#include "fastq.hpp"
#include "Kmer.hpp"
#include "KmerHashTable.h"
#include "KmerIterator.hpp"
#include "minHashIterator.hpp"
#include "NeighborIterator.hpp"
#include "RepHash.hpp"
#include "TinyVector.hpp"
#include "Unitig.hpp"
#include "UnitigIterator.hpp"
#include "UnitigMap.hpp"

#define MASK_CONTIG_ID (0xffffffff00000000)
#define MASK_CONTIG_TYPE (0x80000000)
#define MASK_CONTIG_POS (0x7fffffff)
#define RESERVED_ID (0xffffffff)

#define DEFAULT_K 31
#define DEFAULT_G 23

using namespace std;

template<typename T> class CompactedDBG;

template<typename T>
using f_join = T* (*)(const UnitigMap&, const UnitigMap&, const CompactedDBG<T>&);

template<typename T>
using f_split = vector<T*> (*)(const UnitigMap&, const vector<pair<int,int>>&, const CompactedDBG<T>&);

template<typename T>
using f_add = bool (*)(const UnitigMap&, const string& filename, const size_t, const size_t, const size_t, const CompactedDBG<T>&);

template<typename T = void>
class CompactedDBG {

    private:

        int k_;
        int g_;

        bool invalid;

        const bool has_data;

        static const int tiny_vector_sz = 2;
        static const int min_abundance_lim = 15;
        static const int max_abundance_lim = 15;

        f_join<T> joinData;
        f_split<T> splitData;

        typedef KmerHashTable<CompressedCoverage_t<T>> h_kmers_ccov_t;
        typedef MinimizerHashTable<tiny_vector<size_t,tiny_vector_sz>> hmap_min_unitigs_t;
        typedef typename hmap_min_unitigs_t::iterator hmap_min_unitigs_iterator;
        typedef typename hmap_min_unitigs_t::const_iterator hmap_min_unitigs_const_iterator;
        

        vector<Unitig<T>*> v_unitigs;
        vector<pair<Kmer, CompressedCoverage_t<T>>> v_kmers;

        hmap_min_unitigs_t hmap_min_unitigs;
        h_kmers_ccov_t h_kmers_ccov;

        BlockedBloomFilter bf;

    public:

        typedef T U;

        template<typename U, bool is_const> friend class unitigIterator;
        template<typename U, bool is_const> friend class neighborIterator;

        typedef unitigIterator<T, false> iterator;
        typedef unitigIterator<T, true> const_iterator;

        typedef neighborIterator<T, false> neighbor_iterator;
        typedef neighborIterator<T, true> const_neighbor_iterator;

        const T* getData(const UnitigMap& um) const;
        void setData(const UnitigMap& um, const T* const data);

        CompactedDBG(f_join<T> joinData_ = nullptr, f_split<T> splitData_ = nullptr);
        ~CompactedDBG();

        void setK(int kmer_length = DEFAULT_K, int minimizer_length = DEFAULT_G);
        inline int getK() const { return k_; }

        void clear();
        void empty();

        inline size_t size() const { return v_unitigs.size() + v_kmers.size() + h_kmers_ccov.size(); }

        bool build(const vector<string>& fastx_filename_in, const size_t nb_unique_kmers, const size_t nb_non_unique_kmers,
                    const bool reference_mode = false, const size_t nb_threads = 1, f_add<T> addData = nullptr,
                    const size_t nb_bits_unique_kmers_bf = 14, const size_t nb_bits_non_unique_kmers_bf = 14,
                    const string& inFilenameBBF = "", const string& outFilenameBBF = "", const size_t read_chunksize = 10000,
                    const size_t unitig_size = 1000000, const bool verbose = false);

        bool simplify(const bool delete_short_isolated_unitigs = true, const bool clip_short_tips = true, const bool verbose = false);
        bool write(const string output_filename, const bool verbose = false) const;

        UnitigMap find(const Kmer& km, bool extremities_only = false) const;

        bool add(const string& seq, const bool verbose = false);
        bool remove(const UnitigMap& um, const bool verbose = false);

        string toString(const UnitigMap& um) const;

        Kmer getHead(const UnitigMap& um) const;
        Kmer getTail(const UnitigMap& um) const;

        iterator begin();
        const_iterator begin() const;

        iterator end();
        const_iterator end() const;

        inline BackwardCDBG<T> getPredecessors(const UnitigMap& um) const { return BackwardCDBG<T>(um, *this); }
        inline ForwardCDBG<T> getSuccessors(const UnitigMap& um) const { return ForwardCDBG<T>(um, *this); }

    private:

        bool join(const UnitigMap& um, const bool verbose);

        bool join(const bool verbose);

        bool filter(const vector<string>& fastx_filename_in, const size_t nb_unique_kmers, const size_t nb_non_unique_kmers,
                    const bool reference_mode = false, const size_t nb_threads = 1, const size_t nb_bits_unique_kmers_bf = 14,
                    const size_t nb_bits_non_unique_kmers_bf = 14, const size_t read_chunksize = 10000, const size_t unitig_size = 1000000,
                    const bool verbose = false);

        bool construct(const vector<string>& fastx_filename_in, const size_t nb_threads = 1, f_add<T> addData = nullptr,
                       const size_t read_chunksize = 1000, const bool verbose = false);

        bool addUnitigSequenceBBF(Kmer km, const string& read, size_t pos, const string& seq, vector<Kmer>& l_ignored_km_tip);
        size_t findUnitigSequenceBBF(Kmer km, string& s, bool& selfLoop, bool& isIsolated, vector<Kmer>& l_ignored_km_tip);
        bool bwStepBBF(Kmer km, Kmer& front, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;
        bool fwStepBBF(Kmer km, Kmer& end, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;

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
        size_t joinAllUnitigs(vector<Kmer>* v_joins = nullptr);

        bool checkJoin(const Kmer& a, const UnitigMap& cm_a, Kmer& b);
        void check_fp_tips(KmerHashTable<bool>& ignored_km_tips);
        size_t removeUnitigs(bool rmIsolated, bool clipTips, vector<Kmer>& v);

        void writeGFA(string graphfilename) const;
        void mapRead(const UnitigMap& cc);

        size_t cstrMatch(const char* a, const char* b) const;
        inline size_t stringMatch(const string& a, const string& b, size_t pos);

        neighbor_iterator bw_begin(const UnitigMap& um);
        const_neighbor_iterator bw_begin(const UnitigMap& um) const;
        neighbor_iterator bw_end();
        const_neighbor_iterator bw_end() const;

        neighbor_iterator fw_begin(const UnitigMap& um);
        const_neighbor_iterator fw_begin(const UnitigMap& um) const;
        neighbor_iterator fw_end();
        const_neighbor_iterator fw_end() const;
};

#include "CompactedDBG.tpp"

#endif
