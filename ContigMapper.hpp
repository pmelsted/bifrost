#ifndef BFG_CONTIGMAPPER_HPP
#define BFG_CONTIGMAPPER_HPP

#include "Kmer.hpp"
#include "BloomFilter.hpp"
#include <cstring> // for size_t
#include "HashTables.hpp"
#include "Contig.hpp"
#include "CompressedCoverage.hpp"
#include "ContigMethods.hpp"


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
  ContigMapper(size_t init = 10000) : bf(NULL), limit(Kmer::MAX_K-1) {};
  ~ContigMapper();
  void mapBloomFilter(const BloomFilter* bf);


  ContigMap findContig(Kmer km, const string& s, size_t pos) const;
  void mapRead(const ContigMap& cc);

  bool addContig(Kmer km, const string& read, size_t pos);
  void findContigSequence(Kmer km, string& s);

  size_t contigCount() const;

  
  size_t writeContigs(int count1, string contigfilename, string graphfilename);
  size_t joinAllContigs();
  pair<size_t, size_t> splitAllContigs(); 
  pair<pair<size_t, size_t>, size_t> splitAndJoinAllContigs(); 

  
 private:
  const BloomFilter *bf;
  size_t limit;

  ContigMap find(Kmer km) const;
  bool fwBfStep(Kmer km, Kmer &end, char &c) const;
  bool bwBfStep(Kmer km, Kmer &front, char &c) const;

  size_t joinTwoContigs(ContigRef a, ContigRef b, int a_direction, int b_direction);
  
  typedef google::sparse_hash_map<Kmer, CompressedCoverage, KmerHash> hmap_short_contig_t;
  typedef google::sparse_hash_map<Kmer, Contig*, KmerHash> hmap_long_contig_t;
  typedef google::sparse_hash_map<Kmer, pair<Kmer, size_t>, KmerHash> hmap_shortcut_t; //TODO: switch to smaller type
  
  hmap_short_contig_t sContigs;
  hmap_long_contig_t  lContigs;
  hmap_shortcut_t     shortcuts;
  
};

#endif //BFG_CONTIGMAPPER_HPP
