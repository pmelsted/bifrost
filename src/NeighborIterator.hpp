#ifndef NEIGHBOR_ITERATOR_HPP
#define NEIGHBOR_ITERATOR_HPP

#include "Kmer.hpp"

/** @file src/NeighborIterator.hpp
* The neighborIterator, BackwardCDBG and ForwardCDBG type interfaces.
* Code snippets using these interfaces are provided in snippets.hpp.
*/

template<typename T> class CompactedDBG;
template<typename T> class UnitigMap;

/** @class neighborIterator
* @brief Iterator for the neighbors (predecessors or successors) of a mapped unitig object UnitigMap.
* The first template argument, type T, is the type of data associated with the unitigs in the
* Compacted de Bruijn graph. Second template argument indicates whether the iterator is a constant.
* This iterator considers the forward and reverse-complement strand of the mapped unitig so the mapped
* unitig can have more than 4 predecessors/successors (up to 4 in one direction and 4 in the other
* direction). Also, note that no specific order (such as a lexicographic one) is assumed during iteration.
* An example of using such a class is shown in src/snippets.hpp.
*/
template<typename T = void, bool is_const = true>
class neighborIterator : public std::iterator<std::input_iterator_tag, UnitigMap<T>, int> {

    public:

        typedef typename std::conditional<is_const, const UnitigMap<T>&, UnitigMap<T>&>::type UnitigMap_ref_t;
        typedef typename std::conditional<is_const, const UnitigMap<T>*, UnitigMap<T>*>::type UnitigMap_ptr_t;

        neighborIterator();
        neighborIterator(const UnitigMap<T>& um_, const bool is_forward_);
        neighborIterator(const neighborIterator& o);

        neighborIterator& operator++();
        neighborIterator operator++(int);

        bool operator==(const neighborIterator& o);
        bool operator!=(const neighborIterator& o);

        UnitigMap_ref_t operator*();
        UnitigMap_ptr_t operator->();

    private:

        int i;

        const bool is_fw;
        const bool is_km;

        UnitigMap<T> um;

        Kmer km_head;
        Kmer km_tail;

        CompactedDBG<T>* cdbg;
};

/** @class BackwardCDBG
* @brief Wrapper for class neighborIterator to iterate over the predecessors of a mapped unitig object UnitigMap.
* An example of using such a class is shown in src/snippets.hpp. The first template argument, type T, is the type
* of data associated with the unitigs in the Compacted de Bruijn graph. Second template argument indicates whether
* the iterator accessible through this wrapper is a constant.
*/
template<typename T = void, bool is_const = true>
class BackwardCDBG {

    typedef typename std::conditional<is_const, const UnitigMap<T>&, UnitigMap<T>&>::type UnitigMap_ref_t;

    public:

        explicit BackwardCDBG(UnitigMap_ref_t um_);

        typename UnitigMap<T>::neighbor_iterator begin();
        typename UnitigMap<T>::neighbor_iterator end();

        typename UnitigMap<T>::const_neighbor_iterator begin() const;
        typename UnitigMap<T>::const_neighbor_iterator end() const;

    private:

        UnitigMap_ref_t um;
};

/** @class ForwardCDBG
* @brief Wrapper for class neighborIterator to iterate over the successors of a mapped unitig object UnitigMap.
* An example of using such a class is shown in src/snippets.hpp. The first template argument, type T, is the type
* of data associated with the unitigs in the Compacted de Bruijn graph. Second template argument indicates whether
* the iterator accessible through this wrapper is a constant.
*/
template<typename T = void, bool is_const = true>
class ForwardCDBG {

    typedef typename std::conditional<is_const, const UnitigMap<T>&, UnitigMap<T>&>::type UnitigMap_ref_t;

    public:

        explicit ForwardCDBG(UnitigMap_ref_t um_);

        typename UnitigMap<T>::neighbor_iterator begin();
        typename UnitigMap<T>::neighbor_iterator end();

        typename UnitigMap<T>::const_neighbor_iterator begin() const;
        typename UnitigMap<T>::const_neighbor_iterator end() const;

    private:

        UnitigMap_ref_t um;
};

#include "NeighborIterator.tpp"

#endif
