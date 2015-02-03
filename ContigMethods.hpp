#ifndef BFG_CONTIGMETHODS_HPP
#define BFG_CONTIGMETHODS_HPP

#include <string>
#include "Common.hpp"
#include "Kmer.hpp"
#include "BlockedBloomFilter.hpp"

/* Structs for Contig and Kmer information */

struct ContigMap {
  /**
     xxxxxxxxxxxxxxxxxxxxxxxxxxxxx contig
    |dist ->         xxxxxxxxxyyyyyyy  read on forward strand                
                     | len ->|
          yyyyXXXXXXXX                 read on reverse strand
    | dist -> | len  | 

   */
  Kmer head; // reference to start of contig
  size_t dist; // 0-based distance from start of contig
  size_t len;  // length of match, >= 1
  size_t size; // length of the contig in k-mers
  bool isEmpty; // true if proper match found
  bool isShort; // true if the contig is short
  bool strand; // true for forward strand
  bool selfLoop; // true if this is a self-loop or hairpin
  bool isIsolated;
  bool isTip;    // true if this is a short tip
  Kmer tipHead;  // only used if isTip is true, points to branching k-mer
  ContigMap(Kmer ref, size_t i, size_t l,  size_t sz, bool eq, bool sh) : dist(i), strand(eq), size(sz), len(l), isShort(sh), head(ref), isEmpty(false), selfLoop(false), isTip(false), isIsolated(false) {}
  ContigMap(size_t l = 1) : isEmpty(true), len(l), isTip(false), isIsolated(false) {}
};


struct NewContig {
  Kmer km;
  string read;
  size_t pos;
  NewContig(Kmer o, string &s, size_t p) : km(o), read(s), pos(p) {}
};


#endif // BFG_CONTIGMETHODS_HPP
