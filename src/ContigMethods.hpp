#ifndef BFG_CONTIGMETHODS_HPP
#define BFG_CONTIGMETHODS_HPP

#include <string>
#include "Common.hpp"
#include "Kmer.hpp"
#include "BlockedBloomFilter.hpp"

/* Structs for Unitig and Kmer information */

struct UnitigMap {
    /**
     xxxxxxxxxxxxxxxxxxxxxxxxxxxxx unitig
    |dist ->         xxxxxxxxxyyyyyyy  read on forward strand
                     | len ->|
          yyyyXXXXXXXX                 read on reverse strand
    | dist -> | len  |

    */
    size_t pos_unitig; // unitig pos. in v_unitigs or v_kmers or h_kmers
    size_t pos_min; // position in hmap_min_unitigs_t of the first minimizer of the first k-mer of the match
    size_t dist; // 0-based distance from start of unitig
    size_t len;  // length of match in k-mers, >= 1
    size_t size; // length of the unitig
    bool strand; // true for forward strand
    bool selfLoop; // true if this is a self-loop or hairpin
    bool isEmpty; // true if proper match found
    bool isShort; // true if the unitig has length k
    bool isAbundant; // true if isShort=true and the unitig has an abundant minimizer
    bool isIsolated; // true if unitig is isolated
    bool isTip;    // true if this is a short tip

    UnitigMap(size_t p_unitig, size_t p_min, size_t i, size_t l, size_t sz, bool short_, bool abundance, bool strd) :
            pos_unitig(p_unitig), pos_min(p_min), dist(i), len(l), size(sz), strand(strd), isShort(short_), isAbundant(abundance),
            selfLoop(false), isEmpty(false), isTip(false), isIsolated(false) {}

    UnitigMap(size_t p_min = 0, size_t l = 1) : pos_min(p_min), len(l), isTip(false), isIsolated(false), isShort(false), isAbundant(false), isEmpty(true) {}
};


struct NewUnitig {
  Kmer km;
  string read;
  size_t pos;
  string seq;
  NewUnitig(const Kmer o, const string& s, size_t p, const string& seq_) : km(o), read(s), pos(p), seq(seq_) {}
};

#endif // BFG_CONTIGMETHODS_HPP
