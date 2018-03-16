#ifndef NEIGHBOR_ITERATOR_HPP
#define NEIGHBOR_ITERATOR_HPP

#include "Kmer.hpp"

/** @file src/NeighborIterator.hpp
* The neighborIterator, BackwardCDBG and ForwardCDBG type interfaces.
* Code snippets using these interfaces are provided in snippets.hpp.
*/

template<typename U, typename G> class CompactedDBG;
template<typename U, typename G, bool is_const> class UnitigMap;

/** @class neighborIterator
* @brief Iterator for the neighbors (predecessors or successors) of a mapped unitig object UnitigMap.
* The first template argument, type T, is the type of data associated with the unitigs in the
* Compacted de Bruijn graph. Second template argument indicates whether the iterator is a constant.
* This iterator considers the forward and reverse-complement strand of the mapped unitig so the mapped
* unitig can have more than 4 predecessors/successors (up to 4 in one direction and 4 in the other
* direction). Also, note that no specific order (such as a lexicographic one) is assumed during iteration.
* An example of using such a class is shown in src/snippets.hpp.
*/
template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
class neighborIterator : public std::iterator<std::input_iterator_tag, UnitigMap<Unitig_data_t, Graph_data_t, is_const>, int> {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    public:

        typedef typename std::conditional<is_const, const CompactedDBG<U, G>*, CompactedDBG<U, G>*>::type CompactedDBG_ptr_t;

        neighborIterator();
        neighborIterator(const UnitigMap<U, G, is_const>& um_, const bool is_forward_);
        neighborIterator(const neighborIterator& o);

        neighborIterator& operator++();
        neighborIterator operator++(int);

        bool operator==(const neighborIterator& o) const;
        bool operator!=(const neighborIterator& o) const;

        const UnitigMap<U, G, is_const>& operator*() const;
        const UnitigMap<U, G, is_const>* operator->() const;

    private:

        int i;

        bool is_fw;

        Kmer km_head;
        Kmer km_tail;

        UnitigMap<U, G, is_const> um;

        CompactedDBG_ptr_t cdbg;
};

/** @class BackwardCDBG
* @brief Wrapper for class neighborIterator to iterate over the predecessors of a mapped unitig object UnitigMap.
* An example of using such a class is shown in src/snippets.hpp. The first template argument, type T, is the type
* of data associated with the unitigs in the Compacted de Bruijn graph. Second template argument indicates whether
* the iterator accessible through this wrapper is a constant.
*/
template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
class BackwardCDBG {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    public:

        explicit BackwardCDBG(const UnitigMap<U, G, is_const>& um_);

        neighborIterator<U, G, is_const> begin() const;
        neighborIterator<U, G, is_const> end() const;

    private:

        UnitigMap<U, G, is_const> um;
};

/** @class ForwardCDBG
* @brief Wrapper for class neighborIterator to iterate over the successors of a mapped unitig object UnitigMap.
* An example of using such a class is shown in src/snippets.hpp. The first template argument, type T, is the type
* of data associated with the unitigs in the Compacted de Bruijn graph. Second template argument indicates whether
* the iterator accessible through this wrapper is a constant.
*/
template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
class ForwardCDBG {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    public:

        explicit ForwardCDBG(const UnitigMap<U, G, is_const>& um_);

        neighborIterator<U, G, is_const> begin() const;
        neighborIterator<U, G, is_const> end() const;

    private:

        UnitigMap<U, G, is_const> um;
};

#include "NeighborIterator.tcc"

#endif
