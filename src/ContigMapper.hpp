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

  This class keeps track of all contigs with coverage information as
  necessary.

  For a contig c, we denote the canonical k-mer as the minimum of the
  two representative k-mers at the endpoints.  The class stores shortcut
  information for longer contigs for efficiency reasons.
 */

class ContigMapper {

    public:

        ContigMapper();
        ~ContigMapper();

        void mapBloomFilter(const BlockedBloomFilter *bf);

        ContigMap findContig(const Kmer& km, const string& s, size_t pos) const;
        ContigMap findContig(const Kmer& km, const string& s, size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h) const;
        ContigMap findContig(const Kmer& km, const string& s, size_t pos_km_in_s, size_t pos_min_in_s, size_t it_h);

        ContigMap find(const Kmer& km, bool extremities_only = false) const;

        void mapRead(const ContigMap& cc);

        bool addContigSequence(Kmer km, const string& read, size_t pos, const string& seq, vector<Kmer>& l_ignored_km_tip);

        //size_t findContigSequence(Kmer km, string& s, bool& selfLoop, bool& isIsolated);
        size_t findContigSequence(Kmer km, string& s, bool& selfLoop, bool& isIsolated, vector<Kmer>& l_ignored_km_tip);

        pair<size_t, size_t> splitAllContigs();

        size_t joinAllContigs(vector<Kmer>* v_joins = NULL);
        bool checkJoin(const Kmer& a, const ContigMap& cm_a, Kmer& b/*, bool& dir*/);

        void writeGFA(string graphfilename);

        void checkIntegrity();
        void printState() const;
        size_t contigCount() const;
        void printContigCount() const;

        const BlockedBloomFilter *bf;

        inline size_t find(const preAllocMinHashIterator<RepHash>& it_min_h) const {

            int pos = it_min_h.getPosition();

            return (hmap_min_contigs.find(Minimizer(&it_min_h.s[pos]).rep()) != hmap_min_contigs.end() ? 0 : pos - it_min_h.p);
        }

        void check_fp_tips(KmerHashTable<bool>& ignored_km_tips);

        size_t removeUnitigs(bool rmIsolated, bool clipTips, vector<Kmer>& v);

    private:

        ContigMap find(const Kmer& km, const preAllocMinHashIterator<RepHash>& it_min_h) const;
        ContigMap find(const Kmer& km, const size_t pos_min, const size_t it_h, bool extremities_only = false);

        //bool fwBfStep(Kmer km, Kmer& end, char& c, bool& has_no_neighbor) const;
        //bool bwBfStep(Kmer km, Kmer& front, char& c, bool& has_no_neighbor) const;

        bool fwBfStep(Kmer km, Kmer& end, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;
        bool bwBfStep(Kmer km, Kmer& front, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;

        bool addContig(const string& str_contig, const size_t id_contig);
        void deleteContig(const bool isShort, const bool isAbundant, const size_t id_contig);
        void swapContigs(const bool isShort, const size_t id_a, const size_t id_b);
        bool splitContig(size_t& pos_v_contigs, size_t& nxt_pos_insert_v_contigs, size_t& v_contigs_sz, size_t& v_kmers_sz, const vector<pair<int,int>>& sp);

        static const int tiny_vector_sz = 2;
        static const int min_abundance_lim = 15;

        typedef KmerHashTable<CompressedCoverage> h_kmers_ccov_t;
        typedef MinimizerHashTable<tiny_vector<size_t,tiny_vector_sz>> hmap_min_contigs_t;

        vector<Contig*> v_contigs;
        vector<pair<Kmer, CompressedCoverage>> v_kmers;

        hmap_min_contigs_t hmap_min_contigs;
        h_kmers_ccov_t h_kmers_ccov;
};

#endif //BFG_CONTIGMAPPER_HPP
