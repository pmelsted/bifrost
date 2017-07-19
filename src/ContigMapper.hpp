#ifndef BFG_CONTIGMAPPER_HPP
#define BFG_CONTIGMAPPER_HPP

#include <cstring> // for size_t

#include "Kmer.hpp"
#include "BlockedBloomFilter.hpp"
#include "Contig.hpp"
#include "CompressedCoverage.hpp"
#include "ContigMethods.hpp"
#include "KmerHashTable.h"

#include "RepHash.hpp"
#include "TinyVector.hpp"

#define LOWER_32_MASK  (0xffffffff)
#define UPPER_32_MASK  (0xffffffff00000000)

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

        ContigMap findContig(Kmer km, const string& s, size_t pos) const;
        void mapRead(const ContigMap& cc);

        bool addContigSequence(Kmer km, const string& read, size_t pos, const string& seq);
        size_t findContigSequence(Kmer km, string& s, bool& selfLoop);

        pair<size_t, size_t> splitAllContigs();

        size_t joinAllContigs();
        bool checkJoin(Kmer a, Kmer& b, bool& dir);
        bool checkEndKmer(Kmer b, bool& dir);

        void checkIntegrity();
        void printState() const;
        size_t contigCount() const;
        size_t writeGFA(int count1, string graphfilename);

    private:

        ContigMap find(Kmer km) const;

        bool fwBfStep(Kmer km, Kmer& end, char& c, size_t& deg) const;
        bool bwBfStep(Kmer km, Kmer& front, char& c, size_t& deg) const;

        void addContig(const string& str_contig, const size_t id_contig);
        KmerHashTable<CompressedCoverage>::iterator addShortContig(const string& str_contig);
        void deleteContig(const size_t id_contig);
        KmerHashTable<CompressedCoverage>::iterator deleteShortContig(const size_t id_contig);
        void swapContigs(const size_t id_contig_a, const size_t id_contig_b);

        static const int tiny_vector_sz = 2;

        typedef KmerHashTable<CompressedCoverage> hmap_kmer_contigs_t;
        typedef MinimizerHashTable<tiny_vector<size_t,tiny_vector_sz>> hmap_min_contigs_t;

        const BlockedBloomFilter *bf;
        vector<Contig*> v_contigs;

        hmap_kmer_contigs_t hmap_kmer_contigs;
        hmap_min_contigs_t hmap_min_contigs;
};

#endif //BFG_CONTIGMAPPER_HPP
