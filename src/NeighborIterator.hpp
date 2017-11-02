#ifndef SUPER_NEIGHBOR_ITERATOR_HPP
#define SUPER_NEIGHBOR_ITERATOR_HPP

#include "UnitigMap.hpp"
#include "Kmer.hpp"

template<typename T> class CompactedDBG;

template<typename T = void, bool is_const = true>
class neighborIterator : public std::iterator<std::input_iterator_tag, UnitigMap, int> {

    public:

        typedef typename std::conditional<is_const, const UnitigMap&, UnitigMap&>::type UnitigMap_ref_t;
        typedef typename std::conditional<is_const, const UnitigMap*, UnitigMap*>::type UnitigMap_ptr_t;

        neighborIterator();
        neighborIterator(const Kmer km_, const CompactedDBG<T>* cdbg_, const bool is_forward_);
        neighborIterator(const UnitigMap& um_, const CompactedDBG<T>* cdbg_, const bool is_forward_);
        neighborIterator(const neighborIterator& o);

        neighborIterator& operator++();
        neighborIterator operator++(int);

        bool operator==(const neighborIterator& o);
        bool operator!=(const neighborIterator& o);

        UnitigMap_ref_t operator*();
        UnitigMap_ptr_t operator->();

    protected:

        int i;

        const bool is_fw;

        UnitigMap um;

        Kmer km;

        const CompactedDBG<T>* cdbg;
};

template<typename T>
class BackwardCDBG{

    public:

        explicit BackwardCDBG(const UnitigMap& um_, const CompactedDBG<T>& cdbg_);

        typename CompactedDBG<T>::neighbor_iterator begin();
        typename CompactedDBG<T>::neighbor_iterator end();

        typename CompactedDBG<T>::const_neighbor_iterator begin() const;
        typename CompactedDBG<T>::const_neighbor_iterator end() const;

    private:

        const UnitigMap& um;
        const CompactedDBG<T>& cdbg;
};

template<typename T>
class ForwardCDBG{

    public:

        explicit ForwardCDBG(const UnitigMap& um_, const CompactedDBG<T>& cdbg_);

        typename CompactedDBG<T>::neighbor_iterator begin();
        typename CompactedDBG<T>::neighbor_iterator end();

        typename CompactedDBG<T>::const_neighbor_iterator begin() const;
        typename CompactedDBG<T>::const_neighbor_iterator end() const;

    private:

        const UnitigMap& um;
        const CompactedDBG<T>& cdbg;
};

#include "NeighborIterator.tpp"

#endif
