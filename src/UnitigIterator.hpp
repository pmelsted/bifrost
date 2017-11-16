#ifndef UNITIG_ITERATOR_HPP
#define UNITIG_ITERATOR_HPP

#include "UnitigMap.hpp"
#include "KmerHashTable.h"
#include "CompressedCoverage.hpp"

/** @file src/UnitigIterator.hpp
* The unitigIterator type interface.
* Code snippets using this interface are provided in snippets.hpp.
*/

template<typename T> class CompactedDBG;

/** @class unitigIterator
* @brief Iterator for the unitigs of a Compacted de Bruijn graph.
* The first template argument, type T, is the type of data associated with the unitigs in the
* Compacted de Bruijn graph. Second template argument indicates whether the iterator is a constant.
* Note that no specific order (such as a lexicographic one) is assumed during iteration.
* An example of using such a class is shown in src/snippets.hpp.
*/
template<typename T = void, bool is_const = true>
class unitigIterator : public std::iterator<std::input_iterator_tag, UnitigMap<T>, int> {

    public:

        typedef typename std::conditional<is_const, const UnitigMap<T>&, UnitigMap<T>&>::type UnitigMap_ref_t;
        typedef typename std::conditional<is_const, const UnitigMap<T>*, UnitigMap<T>*>::type UnitigMap_ptr_t;

        unitigIterator();
        unitigIterator(CompactedDBG<T>* cdbg_);
        unitigIterator(const unitigIterator& o);

        unitigIterator& operator++();
        unitigIterator operator++(int);

        bool operator==(const unitigIterator& o);
        bool operator!=(const unitigIterator& o);

        UnitigMap_ref_t operator*();
        UnitigMap_ptr_t operator->();

    private:

        size_t i;

        size_t v_unitigs_sz;
        size_t v_kmers_sz;
        size_t h_kmers_ccov_sz;
        size_t sz;

        bool invalid;

        typename KmerHashTable<CompressedCoverage_t<T>>::const_iterator it_h_kmers_ccov;

        UnitigMap<T> um;

        CompactedDBG<T>* cdbg;
};

#include "UnitigIterator.tpp"

#endif
