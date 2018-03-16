#ifndef UNITIG_ITERATOR_HPP
#define UNITIG_ITERATOR_HPP

#include "UnitigMap.hpp"
#include "KmerHashTable.hpp"
#include "CompressedCoverage.hpp"

/** @file src/UnitigIterator.hpp
* The unitigIterator type interface.
* Code snippets using this interface are provided in snippets.hpp.
*/

template<typename U, typename G> class CompactedDBG;

/** @class unitigIterator
* @brief Iterator for the unitigs of a Compacted de Bruijn graph.
* The first template argument, type T, is the type of data associated with the unitigs in the
* Compacted de Bruijn graph. Second template argument indicates whether the iterator is a constant.
* Note that no specific order (such as a lexicographic one) is assumed during iteration.
* An example of using such a class is shown in src/snippets.hpp.
*/
template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
class unitigIterator : public std::iterator<std::input_iterator_tag, UnitigMap<Unitig_data_t, Graph_data_t, is_const>, int> {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    public:

        typedef typename std::conditional<is_const, const CompactedDBG<U, G>*, CompactedDBG<U, G>*>::type CompactedDBG_ptr_t;

        unitigIterator();
        unitigIterator(CompactedDBG_ptr_t cdbg_);
        unitigIterator(const unitigIterator& o);

        unitigIterator& operator++();
        unitigIterator operator++(int);

        bool operator==(const unitigIterator& o) const;
        bool operator!=(const unitigIterator& o) const;

        const UnitigMap<U, G, is_const>& operator*() const;
        const UnitigMap<U, G, is_const>* operator->() const;

    private:

        size_t i;

        size_t v_unitigs_sz;
        size_t v_kmers_sz;
        size_t h_kmers_ccov_sz;
        size_t sz;

        bool invalid;

        typename KmerHashTable<CompressedCoverage_t<U>>::const_iterator it_h_kmers_ccov;

        UnitigMap<U, G, is_const> um;

        CompactedDBG_ptr_t cdbg;
};

#include "UnitigIterator.tcc"

#endif
