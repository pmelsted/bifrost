#ifndef BFG_UNITIGMAP_HPP
#define BFG_UNITIGMAP_HPP

#include <string>
#include "Common.hpp"
#include "Kmer.hpp"

/** @file src/UnitigMap.hpp
* UnitigMap type interface.
* Code snippets using this interface are provided in snippets/test.cpp.
*/

template<typename U> class Unitig;
template<typename U, typename G> class CompactedDBG;
template<typename U, typename G, bool is_const> class BackwardCDBG;
template<typename U, typename G, bool is_const> class ForwardCDBG;
template<typename U, typename G, bool is_const> class neighborIterator;

/** @struct UnitigMap
* @brief Contain all the information for the mapping of a k-mer or a sequence to a unitig
* of a Compacted de Bruijn graph. A UnitigMap object has 3 template parameters: the type of data
* associated with the unitigs of the graph, the type of data associated with the graph and a boolean
* indicating if this is a constant UnitigMap (const_UnitigMap) or not. A const_UnitigMap can be
* modified but you can't modify the CompactedDBG you can access using UnitigMap::getCompactedDBG.
* The unitig data and graph data types should be the same as the ones used for the CompactedDBG.
* \code{.cpp}
* UnitigMap<> um_1; // No unitig data, no graph data, NOT constant and its content IS NOT constant
* UnitigMap<void, void, false> um_2; // Equivalent to previous notation
* const UnitigMap<> um_3; // No unitig data, no graph data, constant BUT its content IS NOT constant
* UnitigMap<void, void, true> um_4; // No unitig data, no graph data, NOT constant BUT its content IS constant
* const UnitigMap<void, void, true> um_5; // No unitig data, no graph data, constant AND its content IS constant
* UnitigMap<myUnitigData, myGraphData> um_6; // Unitig data of type myUnitigData for each unitig, graph data of type myGraphData, not constant

* CompactedDBG<>* cdbg_ptr_1 = um_1.getCompactedDBG(); // Associated CompactedDBG can be modified from the UnitigMap
* CompactedDBG<>* cdbg_ptr_2 = um_2.getCompactedDBG(); // Associated CompactedDBG can be modified from the UnitigMap
* CompactedDBG<>* cdbg_ptr_3 = um_3.getCompactedDBG(); // Associated CompactedDBG can be modified from the UnitigMap
* const CompactedDBG<>* cdbg_ptr_4 = um_4.getCompactedDBG(); // Associated CompactedDBG cannot be modified from the UnitigMap
* const CompactedDBG<>* cdbg_ptr_5 = um_5.getCompactedDBG(); // Associated CompactedDBG cannot be modified from the UnitigMap
* CompactedDBG<myUnitigData, myGraphData>* cdbg_ptr_6 = um_6.getCompactedDBG(); // Associated CompactedDBG can be modified from the UnitigMap
* \endcode
* @var UnitigMap<T>::isEmpty
* True if there is no mapping.
* @var UnitigMap<T>::dist
* Start position of the mapping (0-based distance) from the start of the reference unitig.
* @var UnitigMap<T>::len
* Length of the mapping on the reference unitig, in k-mers.
* @var UnitigMap<T>::size
* Length of the reference unitig.
* @var UnitigMap<T>::strand
* True if the mapped k-mer or sequence matches the forward strand, false if it matches its reverse-complement.
*/
template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
class UnitigMap {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    template<typename U, typename G> friend class CompactedDBG;
    template<typename U, typename G, bool C> friend class BackwardCDBG;
    template<typename U, typename G, bool C> friend class ForwardCDBG;
    template<typename U, typename G, bool C> friend class unitigIterator;
    template<typename U, typename G, bool C> friend class UnitigMap;

    typedef typename std::conditional<is_const, const CompactedDBG<U, G>*, CompactedDBG<U, G>*>::type CompactedDBG_ptr_t;
    typedef typename std::conditional<is_const, const U*, U*>::type Unitig_data_ptr_t;

    public:

        typedef BackwardCDBG<U, G, is_const> UnitigMap_BW;
        typedef ForwardCDBG<U, G, is_const> UnitigMap_FW;

        UnitigMap(size_t length = 1, CompactedDBG_ptr_t cdbg_ = nullptr);
        UnitigMap(const size_t start, const size_t length, const size_t unitig_sz, const bool strand);

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

        /** Get a pointer to the CompactedDBG containing the reference unitig used in the mapping.
        * @return a pointer to the CompactedDBG containing the reference unitig used in the mapping.
        * If the mapping is empty, a nullptr pointer is returned. The pointer is a constant pointer
        * if the UnitigMap is constant (UnitigMap<U, G, true>).
        */
        inline CompactedDBG_ptr_t getCompactedDBG() const { return cdbg; }

        size_t dist;
        size_t len;
        size_t size;

        bool strand;
        bool isEmpty;

    private:

        UnitigMap(size_t p_unitig, size_t i, size_t l, size_t sz, bool short_, bool abundance, bool strd, CompactedDBG_ptr_t cdbg_);

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

        size_t pos_unitig; // unitig pos. in v_unitigs or v_kmers or h_kmers

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
