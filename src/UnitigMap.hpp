#ifndef SUPER_UNITIGMAP_HPP
#define SUPER_UNITIGMAP_HPP

#include <string>
#include "Common.hpp"
#include "Kmer.hpp"

struct UnitigMap {
    /**
     xxxxxxxxxxxxxxxxxxxxxxxxxxxxx unitig
    |dist ->         xxxxxxxxxyyyyyyy  read on forward strand
                     | len ->|
          yyyyXXXXXXXX                 read on reverse strand
    | dist -> | len  |

    */
    size_t pos_unitig; // unitig pos. in v_unitigs or v_kmers or h_kmers
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

    UnitigMap(size_t p_unitig, size_t i, size_t l, size_t sz, bool short_, bool abundance, bool strd) :
            pos_unitig(p_unitig), dist(i), len(l), size(sz), strand(strd), isShort(short_), isAbundant(abundance),
            selfLoop(false), isEmpty(false), isTip(false), isIsolated(false) {}

    UnitigMap(size_t l = 1) :   len(l), isTip(false), isIsolated(false), isShort(false), isAbundant(false),
                                isEmpty(true), selfLoop(false), strand(true), pos_unitig(0), dist(0), size(0) {}

    bool operator==(const UnitigMap& o) {

        return  (pos_unitig == o.pos_unitig) && (dist == o.dist) && (len == o.len) && (size == o.size) && (strand == o.strand) &&
                (selfLoop == o.selfLoop) && (isEmpty == o.isEmpty) && (isShort == o.isShort) && (isAbundant == o.isAbundant) &&
                (isIsolated == o.isIsolated) && (isTip == o.isTip);
    }

    bool operator!=(const UnitigMap& o) { return !operator==(o); }
};

struct NewUnitig {

    Kmer km;

    string read, seq;

    size_t pos;

    NewUnitig(const Kmer km_, const string& read_, const size_t pos_, const string& seq_) : km(km_), read(read_), pos(pos_), seq(seq_) {}
};

#endif
