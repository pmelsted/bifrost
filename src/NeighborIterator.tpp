#include "CompactedDBG.hpp"

template<typename T, bool is_const>
neighborIterator<T, is_const>::neighborIterator() : i(3), is_fw(true), cdbg(NULL) {}

template<typename T, bool is_const>
neighborIterator<T, is_const>::neighborIterator(const Kmer km_, const CompactedDBG<T>* cdbg_, const bool is_forward_) : i(-1), is_fw(is_forward_), km(km_), cdbg(cdbg_) {

    if ((cdbg == NULL) || cdbg->invalid) i = 3;
}

template<typename T, bool is_const>
neighborIterator<T, is_const>::neighborIterator(const UnitigMap& um_, const CompactedDBG<T>* cdbg_, const bool is_forward_) : i(-1), is_fw(is_forward_), cdbg(cdbg_) {

    if (um_.isEmpty || (cdbg == NULL) || cdbg->invalid) i = 3;
    else if (is_fw) km = cdbg->getTail(um_);
    else km = cdbg->getHead(um_);
}

template<typename T, bool is_const>
neighborIterator<T, is_const>::neighborIterator(const neighborIterator& o) : i(o.i), is_fw(o.is_fw), um(o.um), km(o.km), cdbg(o.cdbg) {}

template<typename T, bool is_const>
neighborIterator<T, is_const>& neighborIterator<T, is_const>::operator++() {

    if ((cdbg == NULL) || cdbg->invalid || (i >= 3)) return *this;

    while (i < 3){

        i++;

        um = cdbg->find(is_fw ? km.forwardBase(alpha[i]) : km.backwardBase(alpha[i]));

        if (!um.isEmpty) break;
    }

    return *this;
}

template<typename T, bool is_const>
neighborIterator<T, is_const> neighborIterator<T, is_const>::operator++(int) {

    neighborIterator tmp(*this);
    operator++();

    return tmp;
}

template<typename T, bool is_const>
bool neighborIterator<T, is_const>::operator==(const neighborIterator& o) {

    if ((i >= 3) || (o.i >= 3)) return (i >= 3) && (o.i >= 3);
    return (is_fw == o.is_fw) && (km == o.km) && (cdbg == o.cdbg) && (um == o.um);
}

template<typename T, bool is_const>
bool neighborIterator<T, is_const>::operator!=(const neighborIterator& o) { return !operator==(o); }

template<typename T, bool is_const>
typename neighborIterator<T, is_const>::UnitigMap_ref_t neighborIterator<T, is_const>::operator*() { return um; }

template<typename T, bool is_const>
typename neighborIterator<T, is_const>::UnitigMap_ptr_t neighborIterator<T, is_const>::operator->() { return &um; }





template<typename T>
BackwardCDBG<T>::BackwardCDBG(const UnitigMap& um_, const CompactedDBG<T>& cdbg_) : um(um_), cdbg(cdbg_) {}

template<typename T>
typename CompactedDBG<T>::neighbor_iterator BackwardCDBG<T>::begin() { return cdbg.bw_begin(um); }

template<typename T>
typename CompactedDBG<T>::neighbor_iterator BackwardCDBG<T>::end() { return cdbg.bw_end(); }

template<typename T>
typename CompactedDBG<T>::const_neighbor_iterator BackwardCDBG<T>::begin() const { return cdbg.bw_begin(um); }

template<typename T>
typename CompactedDBG<T>::const_neighbor_iterator BackwardCDBG<T>::end() const { return cdbg.bw_end(); }





template<typename T>
ForwardCDBG<T>::ForwardCDBG(const UnitigMap& um_, const CompactedDBG<T>& cdbg_) : um(um_), cdbg(cdbg_) {}

template<typename T>
typename CompactedDBG<T>::neighbor_iterator ForwardCDBG<T>::begin() { return cdbg.fw_begin(um); }

template<typename T>
typename CompactedDBG<T>::neighbor_iterator ForwardCDBG<T>::end() { return cdbg.fw_end(); }

template<typename T>
typename CompactedDBG<T>::const_neighbor_iterator ForwardCDBG<T>::begin() const { return cdbg.fw_begin(um); }

template<typename T>
typename CompactedDBG<T>::const_neighbor_iterator ForwardCDBG<T>::end() const { return cdbg.fw_end(); }
