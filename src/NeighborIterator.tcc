#include "CompactedDBG.hpp"
#include "UnitigMap.hpp"

/** Constructor.
* @return an empty neighborIterator.
*/
template<typename T, bool is_const>
neighborIterator<T, is_const>::neighborIterator() : i(4), is_fw(true), cdbg(nullptr) {}

/** Constructor.
* @param um_ is a mapped unitig from which the neighbors are to be iterated
* @param is_forward_ indicates if the iterator must iterate over the successors (true) or the predecessors (false)
* @return a neighborIterator.
*/
template<typename T, bool is_const>
neighborIterator<T, is_const>::neighborIterator(const UnitigMap<T>& um_, const bool is_forward_) : i(-1), is_fw(is_forward_), cdbg(um_.cdbg) {

    if (um_.isEmpty || (um_.cdbg == nullptr) || cdbg->invalid) i = 4;
    else {

        km_head = um_.getHead();
        km_tail = um_.getTail();
    }
}

/** Copy constructor.
* @return a copy of a neighborIterator.
*/
template<typename T, bool is_const>
neighborIterator<T, is_const>::neighborIterator(const neighborIterator& o) : i(o.i), is_fw(o.is_fw), um(o.um), km_head(o.km_head), km_tail(o.km_tail), cdbg(o.cdbg) {}

/** Prefix increment, iterate over the next neighbor (predecessor or successor).
* This iterator considers the forward and reverse-complement strand of the mapped unitig so the mapped
* unitig can have more than 4 predecessors/successors (up to 4 in one direction and 4 in the other
* direction). Also, note that no specific order (such as a lexicographic one) is assumed during iteration.
*/
template<typename T, bool is_const>
neighborIterator<T, is_const>& neighborIterator<T, is_const>::operator++() {

    if ((cdbg == NULL) || cdbg->invalid || (i >= 4)) return *this;

    ++i;

    while (i < 4){

        um = cdbg->find(is_fw ? km_tail.forwardBase(alpha[i]) : km_head.backwardBase(alpha[i]), true);

        if (!um.isEmpty) break;

        ++i;
    }

    return *this;
}

/** Postfix increment, iterate over the next neighbor (predecessor or successor).
* This iterator considers the forward and reverse-complement strand of the mapped unitig so the mapped
* unitig can have more than 4 predecessors/successors (up to 4 in one direction and 4 in the other
* direction). Also, note that no specific order (such as a lexicographic one) is assumed during iteration.
*/
template<typename T, bool is_const>
neighborIterator<T, is_const> neighborIterator<T, is_const>::operator++(int) {

    neighborIterator tmp(*this);
    operator++();

    return tmp;
}

/** Check if two neighborIterators are the same.
* @param o is another neighborIterator.
* @return a boolean indicating whether the two neighborIterators are the same.
*/
template<typename T, bool is_const>
bool neighborIterator<T, is_const>::operator==(const neighborIterator& o) {

    if ((i >= 4) || (o.i >= 4)) return (i >= 4) && (o.i >= 4);
    return (is_fw == o.is_fw) && (km_head == o.km_head) && (km_tail == o.km_tail) && (cdbg == o.cdbg) && (um == o.um);
}

/** Check if two neighborIterators are different.
* @param o is another neighborIterator.
* @return a boolean indicating whether the two neighborIterators are different.
*/
template<typename T, bool is_const>
bool neighborIterator<T, is_const>::operator!=(const neighborIterator& o) { return !operator==(o); }

/** Return a UnitigMap object reference which contains information about the mapped neighbor.
* @return a UnitigMap object reference which contains information about the mapped neighbor.
*/
template<typename T, bool is_const>
typename neighborIterator<T, is_const>::UnitigMap_ref_t neighborIterator<T, is_const>::operator*() { return um; }

/** Return a UnitigMap object pointer which contains information about the mapped neighbor.
* @return a UnitigMap object pointer which contains information about the mapped neighbor.
*/
template<typename T, bool is_const>
typename neighborIterator<T, is_const>::UnitigMap_ptr_t neighborIterator<T, is_const>::operator->() { return &um; }




/** Constructor.
* @param um_ is a reference to a mapped unitig object UnitigMap.
*/
template<typename T, bool is_const>
BackwardCDBG<T, is_const>::BackwardCDBG(BackwardCDBG<T, is_const>::UnitigMap_ref_t um_) : um(um_) {}

/** Return an iterator over the predecessors of a mapped unitig. The returned iterator is initiated
* over the first such predecessor, if there is one.
* @return an iterator over the predecessors of a mapped unitig. It is initiated over the first such
* predecessor, if there is one.
*/
template<typename T, bool is_const>
typename UnitigMap<T>::neighbor_iterator BackwardCDBG<T, is_const>::begin() { return um.bw_begin(); }

/** Return an iterator over the past-the-last predecessor of a mapped unitig.
* @return an iterator over the past-the-last predecessor of a mapped unitig.
*/
template<typename T, bool is_const>
typename UnitigMap<T>::neighbor_iterator BackwardCDBG<T, is_const>::end() { return um.bw_end(); }

/** Return a constant iterator over the predecessors of a mapped unitig. The returned iterator is initiated
* over the first such predecessor, if there is one.
* @return a constant  iterator over the predecessors of a mapped unitig. It is initiated over the first such
* predecessor, if there is one.
*/
template<typename T, bool is_const>
typename UnitigMap<T>::const_neighbor_iterator BackwardCDBG<T, is_const>::begin() const { return um.bw_begin(); }

/** Return a constant iterator over the past-the-last predecessor of a mapped unitig.
* @return a constant iterator over the past-the-last predecessor of a mapped unitig.
*/
template<typename T, bool is_const>
typename UnitigMap<T>::const_neighbor_iterator BackwardCDBG<T, is_const>::end() const { return um.bw_end(); }




/** Constructor.
* @param um_ is a reference to a mapped unitig object UnitigMap.
*/
template<typename T, bool is_const>
ForwardCDBG<T, is_const>::ForwardCDBG(ForwardCDBG<T, is_const>::UnitigMap_ref_t um_) : um(um_) {}

/** Return an iterator over the successors of a mapped unitig. The returned iterator is initiated
* over the first such successor, if there is one.
* @return an iterator over the successors of a mapped unitig. It is initiated over the first such
* successor, if there is one.
*/
template<typename T, bool is_const>
typename UnitigMap<T>::neighbor_iterator ForwardCDBG<T, is_const>::begin() { return um.fw_begin(); }

/** Return an iterator over the past-the-last successor of a mapped unitig.
* @return an iterator over the past-the-last successor of a mapped unitig.
*/
template<typename T, bool is_const>
typename UnitigMap<T>::neighbor_iterator ForwardCDBG<T, is_const>::end() { return um.fw_end(); }

/** Return a constant iterator over the successors of a mapped unitig. The returned iterator is initiated
* over the first such successor, if there is one.
* @return a constant  iterator over the successors of a mapped unitig. It is initiated over the first such
* successor, if there is one.
*/
template<typename T, bool is_const>
typename UnitigMap<T>::const_neighbor_iterator ForwardCDBG<T, is_const>::begin() const { return um.fw_begin(); }

/** Return a constant iterator over the past-the-last successor of a mapped unitig.
* @return a constant iterator over the past-the-last successor of a mapped unitig.
*/
template<typename T, bool is_const>
typename UnitigMap<T>::const_neighbor_iterator ForwardCDBG<T, is_const>::end() const { return um.fw_end(); }
