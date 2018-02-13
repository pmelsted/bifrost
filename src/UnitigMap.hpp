#ifndef UNITIGMAP_HPP
#define UNITIGMAP_HPP

#include <string>
#include "Common.hpp"
#include "Kmer.hpp"

/** @file src/UnitigMap.hpp
* UnitigMap type interface.
* Code snippets using this interface are provided in snippets.hpp.
*/

template<typename T> class CompactedDBG;
template<typename T> class Unitig;
template<typename T, bool is_const> class BackwardCDBG;
template<typename T, bool is_const> class ForwardCDBG;
template<typename T, bool is_const> class neighborIterator;

/** @struct UnitigMap
* @brief Contain all the information for the mapping of a k-mer or a sequence to a unitig
* of a Compacted de Bruijn graph. Its template parameter indicates the type of data associated with
* the unitig and should be the same as the one specified for CompactedDBG. An example of using such
* a structure is shown in src/snippets.hpp.
* @var UnitigMap<T>::isEmpty
* True if the k-mer or sequence does not match a unitig of the graph, false otherwise
* @var UnitigMap<T>::dist
* 0-based distance of the match from start of the unitig
* @var UnitigMap<T>::len
* Length of the match on the unitig (in k-mers)
* @var UnitigMap<T>::size
* Length of the unitig
* @var UnitigMap<T>::strand
* True if the k-mer or sequence matches the forward strand, false if it matches its reverse-complement.
* @var UnitigMap<T>::cdbg
* Compacted de Bruijn graph containing the unitig associated with the mapping.
*/
template<typename T = void>
struct UnitigMap {

    /*

    |size ->                      |
     xxxxxxxxxxxxxxxxxxxxxxxxxxxxx      unitig
    |dist ->         xxxxxxxxxyyyyyyy   read on forward strand
                     | len ->|
          yyyyXXXXXXXX                  read on reverse strand
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
    bool isAbundant; // true if the unitig has length k and has an abundant minimizer
    bool isIsolated; // true if the unitig is isolated
    bool isTip;    // true if this is a short tip

    CompactedDBG<T>* cdbg;

    typedef neighborIterator<T, false> neighbor_iterator;
    typedef neighborIterator<T, true> const_neighbor_iterator;

    UnitigMap(size_t p_unitig, size_t i, size_t l, size_t sz, bool short_, bool abundance, bool strd, CompactedDBG<T>& cdbg_);
    UnitigMap(size_t l = 1);

    bool operator==(const UnitigMap& o);
    bool operator!=(const UnitigMap& o);

    string toString() const;
    string reverseToString() const;

    size_t lcp(const char* s, const size_t pos_s = 0, const size_t pos_um_seq = 0, const bool um_reversed = false) const;

    Kmer getHead() const;
    Kmer getTail() const;
    Kmer getKmer(const size_t pos) const;

    const T* getData() const;
    T* getData();
    void setData(const T* const data) const;

    void mergeData(const UnitigMap<T>& um);
    Unitig<T> splitData(const bool last_split);

    BackwardCDBG<T, true> getPredecessors() const;
    ForwardCDBG<T, true> getSuccessors() const;

    BackwardCDBG<T, false> getPredecessors();
    ForwardCDBG<T, false> getSuccessors();

    neighbor_iterator bw_begin();
    const_neighbor_iterator bw_begin() const;
    neighbor_iterator bw_end();
    const_neighbor_iterator bw_end() const;

    neighbor_iterator fw_begin();
    const_neighbor_iterator fw_begin() const;
    neighbor_iterator fw_end();
    const_neighbor_iterator fw_end() const;
};

///@cond NO_DOC
struct NewUnitig {

    Kmer km;

    string seq;

    size_t len_match_km; // Match length in k-mers

    bool selfLoop;

    NewUnitig(const Kmer km_, const string& seq_, const size_t len_match_km_, const bool selfLoop_) : km(km_), seq(seq_), len_match_km(len_match_km_), selfLoop(selfLoop_) {}
};
///@endcond

#include "UnitigMap.tpp"

#endif
