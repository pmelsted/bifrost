#ifndef UNITIGMAP_HPP
#define UNITIGMAP_HPP

#include <string>
#include "Common.hpp"
#include "Kmer.hpp"

class CompactedDBG;

class rUnitigMap;

template<bool is_const>
class neighborIterator;

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

    const CompactedDBG* dbg;

    UnitigMap(size_t p_unitig, size_t i, size_t l, size_t sz, bool short_, bool abundance, bool strd, const CompactedDBG* dbg_) :
            pos_unitig(p_unitig), dist(i), len(l), size(sz), strand(strd), isShort(short_), isAbundant(abundance),
            selfLoop(false), isEmpty(false), isTip(false), isIsolated(false), dbg(dbg_) {}

    UnitigMap(size_t l = 1, const CompactedDBG* dbg_ = NULL) :  len(l), isTip(false), isIsolated(false), isShort(false),
                                                                isAbundant(false), isEmpty(true), selfLoop(false), strand(true),
                                                                pos_unitig(0), dist(0), size(0), dbg(dbg_) {}

    bool operator==(const UnitigMap& o) {

        return  (pos_unitig == o.pos_unitig) && (dist == o.dist) && (len == o.len) && (size == o.size) && (strand == o.strand) &&
                (selfLoop == o.selfLoop) && (isEmpty == o.isEmpty) && (isShort == o.isShort) && (isAbundant == o.isAbundant) &&
                (isIsolated == o.isIsolated) && (isTip == o.isTip);
    }

    bool operator!=(const UnitigMap& o) { return !operator==(o); }

    string toString() const;

    Kmer getHead() const;
    Kmer getTail() const;

    typedef neighborIterator<true> iterator;
    typedef neighborIterator<false> const_iterator;

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    iterator rbegin();
    iterator rend();

    const_iterator rbegin() const;
    const_iterator rend() const;

    rUnitigMap backward() const;
};

class rUnitigMap{

    public:

        explicit rUnitigMap(const UnitigMap& um) : um_(um) {}

        UnitigMap::const_iterator begin() const;
        UnitigMap::const_iterator end() const;

    private:

        const UnitigMap um_;
};

struct NewUnitig {

    Kmer km;

    string read, seq;

    size_t pos;

    NewUnitig(const Kmer o, const string& s, size_t p, const string& seq_) : km(o), read(s), pos(p), seq(seq_) {}
};

template<bool is_const = true>
class neighborIterator : public std::iterator<std::input_iterator_tag, UnitigMap, int> {

    public:

        typedef typename std::conditional<is_const, const UnitigMap&, UnitigMap&>::type UnitigMap_ref_t;
        typedef typename std::conditional<is_const, const UnitigMap*, UnitigMap*>::type UnitigMap_ptr_t;

        neighborIterator();
        neighborIterator(const Kmer km, const CompactedDBG* dbg, const bool is_forward);
        neighborIterator(const UnitigMap& um, const CompactedDBG* dbg, const bool is_forward);
        neighborIterator(const neighborIterator& o);

        neighborIterator& operator++();
        neighborIterator operator++(int);

        bool operator==(const neighborIterator& o);
        bool operator!=(const neighborIterator& o);

        UnitigMap_ref_t operator*();
        UnitigMap_ptr_t operator->();

    protected:

        int i_;

        const bool is_fw_;

        UnitigMap um_;

        Kmer km_;

        const CompactedDBG* dbg_;
};

#endif
