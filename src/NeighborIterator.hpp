#ifndef SUPER_NEIGHBOR_ITERATOR_HPP
#define SUPER_NEIGHBOR_ITERATOR_HPP

#include "Kmer.hpp"

template<typename T> class CompactedDBG;
template<typename T> class UnitigMap;

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

    protected:

        int i;

        const bool is_fw;
        const bool is_km;

        UnitigMap<T> um;

        Kmer km_head;
        Kmer km_tail;

        CompactedDBG<T>* cdbg;
};

template<typename T = void, bool is_const = true>
class BackwardCDBG{

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

template<typename T = void, bool is_const = true>
class ForwardCDBG{

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
