#include "CompactedDBG.hpp"
#include "UnitigMap.hpp"

template<typename T, bool is_const>
neighborIterator<T, is_const>::neighborIterator() : i(7), is_fw(true), is_km(false), cdbg(nullptr) {}

template<typename T, bool is_const>
neighborIterator<T, is_const>::neighborIterator(const UnitigMap<T>& um_, const bool is_forward_) : i(-1), is_fw(is_forward_), is_km(um_.isShort || um_.isAbundant), cdbg(um_.cdbg) {

    if (um_.isEmpty || (um_.cdbg == nullptr) || cdbg->invalid) i = 7;
    else if (is_fw){

        km_tail = um_.getTail();
        km_head = um_.getHead().twin();
    }
    else {

        km_head = um_.getHead();
        km_tail = um_.getTail().twin();
    }
}

template<typename T, bool is_const>
neighborIterator<T, is_const>::neighborIterator(const neighborIterator& o) : i(o.i), is_fw(o.is_fw), is_km(o.is_km), um(o.um), km_head(o.km_head), km_tail(o.km_tail), cdbg(o.cdbg) {}

template<typename T, bool is_const>
neighborIterator<T, is_const>& neighborIterator<T, is_const>::operator++() {

    if ((cdbg == NULL) || cdbg->invalid || (i >= 7)) return *this;

    while (i < 7){

        i++;

        if (i <= 3) um = cdbg->find(is_fw ? km_tail.forwardBase(alpha[i]) : km_head.backwardBase(alpha[i]));
        else if (is_km){

            i = 7;
            um = UnitigMap<T>();
        }
        else um = cdbg->find(is_fw ? km_head.forwardBase(alpha[i-4]) : km_tail.backwardBase(alpha[i-4]));

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

    if ((i >= 7) || (o.i >= 7)) return (i >= 7) && (o.i >= 7);
    return (is_fw == o.is_fw) && (km_head == o.km_head) && (km_tail == o.km_tail) && (cdbg == o.cdbg) && (um == o.um);
}

template<typename T, bool is_const>
bool neighborIterator<T, is_const>::operator!=(const neighborIterator& o) { return !operator==(o); }

template<typename T, bool is_const>
typename neighborIterator<T, is_const>::UnitigMap_ref_t neighborIterator<T, is_const>::operator*() { return um; }

template<typename T, bool is_const>
typename neighborIterator<T, is_const>::UnitigMap_ptr_t neighborIterator<T, is_const>::operator->() { return &um; }





template<typename T, bool is_const>
BackwardCDBG<T, is_const>::BackwardCDBG(BackwardCDBG<T, is_const>::UnitigMap_ref_t um_) : um(um_) {}

template<typename T, bool is_const>
typename UnitigMap<T>::neighbor_iterator BackwardCDBG<T, is_const>::begin() { return um.bw_begin(); }

template<typename T, bool is_const>
typename UnitigMap<T>::neighbor_iterator BackwardCDBG<T, is_const>::end() { return um.bw_end(); }

template<typename T, bool is_const>
typename UnitigMap<T>::const_neighbor_iterator BackwardCDBG<T, is_const>::begin() const { return um.bw_begin(); }

template<typename T, bool is_const>
typename UnitigMap<T>::const_neighbor_iterator BackwardCDBG<T, is_const>::end() const { return um.bw_end(); }





template<typename T, bool is_const>
ForwardCDBG<T, is_const>::ForwardCDBG(ForwardCDBG<T, is_const>::UnitigMap_ref_t um_) : um(um_) {}

template<typename T, bool is_const>
typename UnitigMap<T>::neighbor_iterator ForwardCDBG<T, is_const>::begin() { return um.fw_begin(); }

template<typename T, bool is_const>
typename UnitigMap<T>::neighbor_iterator ForwardCDBG<T, is_const>::end() { return um.fw_end(); }

template<typename T, bool is_const>
typename UnitigMap<T>::const_neighbor_iterator ForwardCDBG<T, is_const>::begin() const { return um.fw_begin(); }

template<typename T, bool is_const>
typename UnitigMap<T>::const_neighbor_iterator ForwardCDBG<T, is_const>::end() const { return um.fw_end(); }
