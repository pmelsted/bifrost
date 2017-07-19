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
  size_t pos_contig; // position of the contig in vector v_contigs or in hash table hmap_kmer_contigs
  size_t dist; // 0-based distance from start of contig
  size_t len;  // length of match in k-mers, >= 1
  size_t size; // length of the contig in k-mers
  bool strand; // true for forward strand
  bool isEmpty; // true if proper match found
  bool isShort; // true if the contig has length k
  bool selfLoop; // true if this is a self-loop or hairpin
  bool isIsolated; // true if contig is isolated
  bool isTip;    // true if this is a short tip
  Kmer tipHead;  // only used if isTip is true, points to branching k-mer
  ContigMap(size_t p, size_t i, size_t l, size_t sz, bool short_, bool strd) : pos_contig(p), dist(i), len(l), size(sz), strand(strd), selfLoop(false), isEmpty(false), isShort(short_), isTip(false), isIsolated(false) {}
  ContigMap(size_t l = 1) : len(l), isTip(false), isIsolated(false), isShort(true), isEmpty(true) {}
};


struct NewContig {
  Kmer km;
  string read;
  size_t pos;
  string seq;
  NewContig(const Kmer o, const string& s, size_t p, const string& seq_) : km(o), read(s), pos(p), seq(seq_) {}
};


#endif // BFG_CONTIGMETHODS_HPP
