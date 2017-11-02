#ifndef SUPER_UNITIG_ITERATOR_HPP
#define SUPER_UNITIG_ITERATOR_HPP

#include "UnitigMap.hpp"
#include "KmerHashTable.h"
#include "CompressedCoverage.hpp"

template<typename T> class superCompactedDBG;

template<typename T = void, bool is_const = true>
class unitigIterator : public std::iterator<std::input_iterator_tag, UnitigMap, int> {

    public:

        typedef typename std::conditional<is_const, const UnitigMap&, UnitigMap&>::type UnitigMap_ref_t;
        typedef typename std::conditional<is_const, const UnitigMap*, UnitigMap*>::type UnitigMap_ptr_t;

        unitigIterator();
        unitigIterator(const CompactedDBG<T>* cdbg_);
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

        UnitigMap um;

        const CompactedDBG<T>* cdbg;
};

#include "UnitigIterator.tpp"

#endif
