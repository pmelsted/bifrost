#ifndef BFG_CONTIGMAPPER_HPP
#define BFG_CONTIGMAPPER_HPP

#include <cstring> // for size_t

#include "Kmer.hpp"
#include "BlockedBloomFilter.hpp"
#include "Contig.hpp"
#include "CompressedCoverage.hpp"
#include "ContigMethods.hpp"
#include "KmerHashTable.h"

#include "minHashIterator.hpp"
#include "RepHash.hpp"
#include "TinyVector.hpp"

#define MASK_CONTIG_ID (0xffffffff00000000)
#define MASK_CONTIG_TYPE (0x80000000)
#define MASK_CONTIG_POS (0x7fffffff)

#define RESERVED_ID (0xffffffff)

/*
  Short description:

  This class keeps track of all unitigs with coverage information as
  necessary.

  For a unitig c, we denote the canonical k-mer as the minimum of the
  two representative k-mers at the endpoints.  The class stores shortcut
  information for longer unitigs for efficiency reasons.
 */

class UnitigMapper {

    public:

        UnitigMapper();
        ~UnitigMapper();

        void empty();

        UnitigMap findUnitig(const Kmer& km, const string& s, size_t pos) const;
        UnitigMap findUnitig(const Kmer& km, const string& s, size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h) const;

        UnitigMap find(const Kmer& km, bool extremities_only = false) const;

        inline size_t find(const preAllocMinHashIterator<RepHash>& it_min_h) const {

            const int pos = it_min_h.getPosition();

            return (hmap_min_unitigs.find(Minimizer(&it_min_h.s[pos]).rep()) != hmap_min_unitigs.end() ? 0 : pos - it_min_h.p);
        }

        void mapRead(const UnitigMap& cc);

        bool addUnitigSequence(Kmer km, const string& read, size_t pos, const string& seq, vector<Kmer>& l_ignored_km_tip);
        size_t findUnitigSequence(Kmer km, string& s, bool& selfLoop, bool& isIsolated, vector<Kmer>& l_ignored_km_tip);

        pair<size_t, size_t> splitAllUnitigs();

        size_t joinAllUnitigs(vector<Kmer>* v_joins = NULL);
        bool checkJoin(const Kmer& a, const UnitigMap& cm_a, Kmer& b);

        void writeGFA(string graphfilename) const;

        void printState() const;
        size_t unitigCount() const;
        void printUnitigCount() const;

        void check_fp_tips(KmerHashTable<bool>& ignored_km_tips);

        size_t removeUnitigs(bool rmIsolated, bool clipTips, vector<Kmer>& v);

        void mapBloomFilter(const BlockedBloomFilter *bf);

        const BlockedBloomFilter *bf;

    private:

        UnitigMap find(const Kmer& km, const preAllocMinHashIterator<RepHash>& it_min_h) const;

        bool fwBfStep(Kmer km, Kmer& end, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;
        bool bwBfStep(Kmer km, Kmer& front, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;

        bool addUnitig(const string& str_unitig, const size_t id_unitig);
        void deleteUnitig(const bool isShort, const bool isAbundant, const size_t id_unitig);
        void swapUnitigs(const bool isShort, const size_t id_a, const size_t id_b);

        bool splitUnitig(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs, size_t& v_unitigs_sz, size_t& v_kmers_sz,
                         const vector<pair<int,int>>& sp);

        static const int tiny_vector_sz = 2;
        static const int min_abundance_lim = 15;
        static const int max_abundance_lim = 15;

        typedef KmerHashTable<CompressedCoverage> h_kmers_ccov_t;
        typedef MinimizerHashTable<tiny_vector<size_t,tiny_vector_sz>> hmap_min_unitigs_t;

        vector<Unitig*> v_unitigs;
        vector<pair<Kmer, CompressedCoverage>> v_kmers;

        hmap_min_unitigs_t hmap_min_unitigs;
        h_kmers_ccov_t h_kmers_ccov;
};

#endif //BFG_CONTIGMAPPER_HPP
