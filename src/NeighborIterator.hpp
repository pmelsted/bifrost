#ifndef NEIGHBOR_ITERATOR_HPP
#define NEIGHBOR_ITERATOR_HPP

#include "Kmer.hpp"

/** @file src/NeighborIterator.hpp
* The neighborIterator, BackwardCDBG and ForwardCDBG type interfaces.
* Code snippets using these interfaces are provided in snippets/test.cpp.
*/

template<typename U, typename G> class CompactedDBG;
template<typename U, typename G, bool is_const> class UnitigMap;

/** @class neighborIterator
* @brief Iterator for the neighbors (predecessors or successors) of a reference unitig used in a UnitigMap object.
* A neighborIterator object has 3 template parameters: the type of data associated with the unitigs of the graph,
* the type of data associated with the graph and a boolean indicating if this is a constant iterator or not.
* Note that you are supposed to use this class as the iterator of a BackwardCDBG or ForwardCDBG object, which can
* be obtained respectively from UnitigMap::getPredecessors() and UnitigMap::getSuccessors(), so you shouldn't
* have to instantiate an object neighborIterator and its template parameters yourself. The unitig data and graph
* data types should be the same as the ones used for the CompactedDBG the iterator is from. No specific order
* (such as a lexicographic one) is assumed during iteration.
* \code{.cpp}
* CompactedDBG<> cdbg;
* ... // Some more code, cdbg construction
* for (const auto& unitig : cdbg){
*   cout << unitig.toString() << endl; // unitig is of type UnitigMap
*   for (const auto& pred : unitig.getPredecessors()) cout << pred.toString() << endl;
*   for (const auto& succ : unitig.getSuccessors()) cout << succ.toString() << endl;
* }
* \endcode
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
* @brief Wrapper for class neighborIterator to iterate over the predecessors of a reference unitig used in a UnitigMap object.
* A BackwardCDBG object has 3 template parameters: the type of data associated with the unitigs of the graph,
* the type of data associated with the graph and a boolean indicating if this is a constant neighborIterator or not.
* Note that you are supposed to obtain an instance of this class from UnitigMap::getPredecessors() so you shouldn't
* have to instantiate an object BackwardCDBG and its template parameters yourself. The unitig data and graph
* data types should be the same as the ones used for the CompactedDBG the iterator is from. No specific order
* (such as a lexicographic one) is assumed during iteration.
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
* @brief Wrapper for class neighborIterator to iterate over the predecessors of a reference unitig used in a UnitigMap object.
* A ForwardCDBG object has 3 template parameters: the type of data associated with the unitigs of the graph,
* the type of data associated with the graph and a boolean indicating if this is a constant neighborIterator or not.
* Note that you are supposed to obtain an instance of this class from UnitigMap::getSuccessors() so you shouldn't
* have to instantiate an object ForwardCDBG and its template parameters yourself. The unitig data and graph
* data types should be the same as the ones used for the CompactedDBG the iterator is from. No specific order
* (such as a lexicographic one) is assumed during iteration.
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
