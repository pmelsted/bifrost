#ifndef UNITIGMAP_HPP
#define UNITIGMAP_HPP

#include <string>
#include "Common.hpp"
#include "Kmer.hpp"

/** @file src/UnitigMap.hpp
* UnitigMap type interface.
* Code snippets using this interface are provided in snippets.hpp.
*/

template<typename U> class Unitig;
template<typename U, typename G> class CompactedDBG;
template<typename U, typename G, bool is_const> class BackwardCDBG;
template<typename U, typename G, bool is_const> class ForwardCDBG;
template<typename U, typename G, bool is_const> class neighborIterator;

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

template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
class UnitigMap {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    template<typename U, typename G> friend class CompactedDBG;
    template<typename U, typename G, bool is_const_bis> friend class BackwardCDBG;
    template<typename U, typename G, bool is_const_bis> friend class ForwardCDBG;
    template<typename U, typename G, bool is_const_bis> friend class UnitigMap;

    typedef typename std::conditional<is_const, const CompactedDBG<U, G>*, CompactedDBG<U, G>*>::type CompactedDBG_ptr_t;
    typedef typename std::conditional<is_const, const U*, U*>::type Unitig_data_ptr_t;

    public:

        typedef BackwardCDBG<U, G, is_const> UnitigMap_BW;
        typedef ForwardCDBG<U, G, is_const> UnitigMap_FW;

        UnitigMap(size_t p_unitig, size_t i, size_t l, size_t sz, bool short_, bool abundance, bool strd, CompactedDBG_ptr_t cdbg_);
        UnitigMap(size_t l = 1, CompactedDBG_ptr_t cdbg_ = nullptr);

        operator UnitigMap<U, G, true>() const {

            UnitigMap<U, G, true> um(pos_unitig, dist, len, size, isShort, isAbundant, strand, cdbg);

            um.isEmpty = isEmpty;

            return um;
        }

        bool operator==(const UnitigMap& o) const;
        bool operator!=(const UnitigMap& o) const;

        string toString() const;

        size_t lcp(const char* s, const size_t pos_s = 0, const size_t pos_um_seq = 0, const bool um_reversed = false) const;

        Kmer getUnitigHead() const;
        Kmer getUnitigTail() const;
        Kmer getUnitigKmer(const size_t pos) const;

        Kmer getMappedHead() const;
        Kmer getMappedTail() const;
        Kmer getMappedKmer(const size_t pos) const;

        UnitigMap<U, G, is_const> getKmerMapping(const size_t pos) const;

        Unitig_data_ptr_t getData() const;

        UnitigMap_BW getPredecessors() const;
        UnitigMap_FW getSuccessors() const;

        inline CompactedDBG_ptr_t getCompactedDBG() const { return cdbg; }

        size_t pos_unitig; // unitig pos. in v_unitigs or v_kmers or h_kmers
        size_t dist; // 0-based distance from start of unitig
        size_t len;  // length of match in k-mers, >= 1
        size_t size; // length of the unitig

        bool strand; // true for forward strand
        bool isEmpty; // true if proper match found

    private:

        neighborIterator<U, G, is_const> bw_begin() const;
        neighborIterator<U, G, is_const> bw_end() const;

        neighborIterator<U, G, is_const> fw_begin() const;
        neighborIterator<U, G, is_const> fw_end() const;

        template<bool is_void> typename std::enable_if<!is_void, void>::type mergeData_(const UnitigMap<U, G, is_const>& um) const;
        template<bool is_void> typename std::enable_if<is_void, void>::type mergeData_(const UnitigMap<U, G, is_const>& um) const;

        void mergeData(const UnitigMap<U, G, is_const>& um) const;

        template<bool is_void> typename std::enable_if<!is_void, Unitig<U>>::type splitData_(const bool last_split) const;
        template<bool is_void> typename std::enable_if<is_void, Unitig<U>>::type splitData_(const bool last_split) const;

        Unitig<U> splitData(const bool last_split) const;

        template<bool is_void> typename std::enable_if<!is_void, Unitig_data_ptr_t>::type getData_() const;
        template<bool is_void> typename std::enable_if<is_void, Unitig_data_ptr_t>::type getData_() const;

        void partialCopy(const UnitigMap<U, G, is_const>& um);

        bool isShort; // true if the unitig has length k
        bool isAbundant; // true if the unitig has length k and has an abundant minimizer

        CompactedDBG_ptr_t cdbg;
};

template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
struct UnitigMapHash {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    size_t operator()(const UnitigMap<U, G, is_const>& um) const {

        return static_cast<size_t>(XXH64(static_cast<const void*>(&um), sizeof(UnitigMap<U, G, is_const>), 0));
    }
};

#include "UnitigMap.tcc"

#endif
