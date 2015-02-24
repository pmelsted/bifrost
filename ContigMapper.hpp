#ifndef BFG_CONTIGMAPPER_HPP
#define BFG_CONTIGMAPPER_HPP

#include "Kmer.hpp"
#include "BlockedBloomFilter.hpp"
#include <cstring> // for size_t
#include "Contig.hpp"
#include "CompressedCoverage.hpp"
#include "ContigMethods.hpp"
#include "KmerHashTable.h"


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
  ContigMapper(size_t init = 10000);
  ~ContigMapper();
  void mapBloomFilter(const BlockedBloomFilter *bf);


  ContigMap findContig(Kmer km, const string& s, size_t pos) const;
  void mapRead(const ContigMap& cc);

  bool addContig(Kmer km, const string& read, size_t pos, const string& seq);
  void findContigSequence(Kmer km, string& s, bool& selfLoop);

  size_t contigCount() const;


  size_t writeGFA(int count1, string graphfilename, bool debug);
  size_t joinAllContigs();
  pair<size_t, size_t> splitAllContigs();
  size_t clipTips();

  void moveShortContigs();
  void fixShortContigs();
  size_t removeIsolatedContigs();

  bool checkTip(Kmer tip);

  bool checkJoin(Kmer a, Kmer& b, bool& dir);
  bool checkEndKmer(Kmer b, bool& dir);

  bool checkShortcuts();
  void setStride(size_t stride_) { stride = stride_; }
  void printState() const;

 private:
  const BlockedBloomFilter *bf;
  size_t limit;
  size_t stride;

  void removeShortcuts(const string& s);

  ContigMap find(Kmer km) const;
  bool fwBfStep(Kmer km, Kmer& end, char& c, size_t& deg) const;
  bool bwBfStep(Kmer km, Kmer& front, char& c, size_t& deg) const;




  typedef KmerHashTable<CompressedCoverage> hmap_short_contig_t;
  typedef KmerHashTable<Contig *> hmap_long_contig_t;
  typedef KmerHashTable<pair<Kmer, size_t>> hmap_shortcut_t;

  hmap_short_contig_t sContigs;
  hmap_long_contig_t  lContigs;
  hmap_shortcut_t     shortcuts;

};

#endif //BFG_CONTIGMAPPER_HPP
