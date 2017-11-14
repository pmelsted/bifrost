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

struct CDBG_Build_opt {

    bool reference_mode; // Reference mode
    bool verbose; // Print messages during running time
    bool clipTips; // Clip short (<2k) tips
    bool deleteIsolated; // Delete isolated short (<2k) unitigs

    size_t nb_threads; // Number of threads to use for building
    size_t k, g; // Length of k-mers and g-mers (minimizers)
    size_t read_chunksize; // Number of reads shared and processed by x threads at the same time
    size_t unitig_size; // Maximum length of a unitig
    size_t nb_unique_kmers;
    size_t nb_non_unique_kmers;
    size_t nb_bits_unique_kmers_bf;
    size_t nb_bits_non_unique_kmers_bf;

    string prefixFilenameGFA;
    string filenameGFA;
    string inFilenameBBF;
    string outFilenameBBF;

    vector<string> fastx_filename_in;

    CDBG_Build_opt() :  nb_threads(1), k(31), g(23), nb_unique_kmers(0), nb_non_unique_kmers(0), nb_bits_unique_kmers_bf(14),
                        nb_bits_non_unique_kmers_bf(14), read_chunksize(10000), unitig_size(1000000), reference_mode(false),
                        verbose(false), clipTips(false), deleteIsolated(false), inFilenameBBF(""), outFilenameBBF("") {}
};

template<typename T> //Curiously Recurring Template Pattern (CRTP)
class CDBG_Data_t {

    virtual void join(const T& data, CompactedDBG<T>& cdbg) = 0;
    virtual void split(const size_t pos, const size_t len, T& new_data, CompactedDBG<T>& cdbg) const;
};

template<typename T = void>
class CompactedDBG {

    static_assert(is_base_of<CDBG_Data_t<T>, T>::value || is_void<T>::value,
                  "Type of data associated with vertices of class CompactedDBG must be void (no data) or a class extending class CDBG_Data_t");

    private:

        int k_;
        int g_;

        bool invalid;

        const bool has_data;

        static const int tiny_vector_sz = 2;
        static const int min_abundance_lim = 15;
        static const int max_abundance_lim = 15;

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

        template<typename U> friend class UnitigMap;
        template<typename U, bool is_const> friend class unitigIterator;
        template<typename U, bool is_const> friend class neighborIterator;

        typedef unitigIterator<T, false> iterator;
        typedef unitigIterator<T, true> const_iterator;

        CompactedDBG(int kmer_length = DEFAULT_K, int minimizer_length = DEFAULT_G);
        ~CompactedDBG();

        inline int getK() const { return k_; }

        void clear();
        void empty();

        inline size_t size() const { return v_unitigs.size() + v_kmers.size() + h_kmers_ccov.size(); }

        bool build(const CDBG_Build_opt& opt);
        bool simplify(const bool delete_short_isolated_unitigs = true, const bool clip_short_tips = true, const bool verbose = false);
        bool write(const string output_filename, const bool verbose = false);

        UnitigMap<T> find(const Kmer& km, bool extremities_only = false);

        bool add(const string& seq, const bool verbose = false);
        bool remove(const UnitigMap<T>& um, const bool verbose = false);

        iterator begin();
        const_iterator begin() const;

        iterator end();
        const_iterator end() const;

    private:

        bool join(const UnitigMap<T>& um, const bool verbose);

        bool join(const bool verbose);

        bool filter(const CDBG_Build_opt& opt);
        bool construct(const CDBG_Build_opt& opt);

        bool addUnitigSequenceBBF(Kmer km, const string& read, size_t pos, const string& seq, vector<Kmer>& l_ignored_km_tip);
        size_t findUnitigSequenceBBF(Kmer km, string& s, bool& selfLoop, bool& isIsolated, vector<Kmer>& l_ignored_km_tip);
        bool bwStepBBF(Kmer km, Kmer& front, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;
        bool fwStepBBF(Kmer km, Kmer& end, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;

        UnitigMap<T> findUnitig(const Kmer& km, const string& s, size_t pos);
        UnitigMap<T> findUnitig(const Kmer& km, const string& s, size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h);

        bool addUnitig(const string& str_unitig, const size_t id_unitig);
        void deleteUnitig(const bool isShort, const bool isAbundant, const size_t id_unitig);
        void swapUnitigs(const bool isShort, const size_t id_a, const size_t id_b);
        bool splitUnitig(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs, size_t& v_unitigs_sz, size_t& v_kmers_sz,
                        const vector<pair<int,int>>& sp);

        UnitigMap<T> find(const Kmer& km, const preAllocMinHashIterator<RepHash>& it_min_h);

        inline size_t find(const preAllocMinHashIterator<RepHash>& it_min_h) const {

            const int pos = it_min_h.getPosition();
            return (hmap_min_unitigs.find(Minimizer(&it_min_h.s[pos]).rep()) != hmap_min_unitigs.end() ? 0 : pos - it_min_h.p);
        }

        pair<size_t, size_t> splitAllUnitigs();
        size_t joinAllUnitigs(vector<Kmer>* v_joins = nullptr);

        bool checkJoin(const Kmer& a, const UnitigMap<T>& cm_a, Kmer& b);
        void check_fp_tips(KmerHashTable<bool>& ignored_km_tips);
        size_t removeUnitigs(bool rmIsolated, bool clipTips, vector<Kmer>& v);

        void writeGFA(string graphfilename);
        void mapRead(const UnitigMap<T>& cc);

        size_t cstrMatch(const char* a, const char* b) const;
        inline size_t stringMatch(const string& a, const string& b, size_t pos);
};

#include "CompactedDBG.tpp"

#endif
