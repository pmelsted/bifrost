#include "CompactedDBG.hpp"
#include "NeighborIterator.hpp"

template<typename T>
UnitigMap<T>::UnitigMap(size_t p_unitig, size_t i, size_t l, size_t sz, bool short_, bool abundance, bool strd, CompactedDBG<T>& cdbg_) :
                        pos_unitig(p_unitig), dist(i), len(l), size(sz), cdbg(&cdbg_), strand(strd), isShort(short_), isAbundant(abundance),
                        selfLoop(false), isEmpty(false), isTip(false), isIsolated(false) {}

/** UnitigMap constructor.
* @return an empty UnitigMap.
*/
template<typename T>
UnitigMap<T>::UnitigMap(size_t l) : len(l), isTip(false), isIsolated(false), isShort(false), isAbundant(false), isEmpty(true),
                                    selfLoop(false), strand(true), pos_unitig(0), dist(0), size(0), cdbg(nullptr) {}

/** Check if two UnitigMaps are the same.
* @return a boolean indicating if the two UnitigMaps are the same.
*/
template<typename T>
bool UnitigMap<T>::operator==(const UnitigMap& o) {

    return  (pos_unitig == o.pos_unitig) && (dist == o.dist) && (len == o.len) && (size == o.size) && (strand == o.strand) &&
            (selfLoop == o.selfLoop) && (isEmpty == o.isEmpty) && (isShort == o.isShort) && (isAbundant == o.isAbundant) &&
            (isIsolated == o.isIsolated) && (isTip == o.isTip) && (cdbg == o.cdbg);
}

/** Check if two UnitigMaps are different.
* @return a boolean indicating if the two UnitigMaps are different.
*/
template<typename T> bool UnitigMap<T>::operator!=(const UnitigMap& o) { return !operator==(o); }

/** Return a string containing the sequence of the mapped unitig.
* @return a string containing the sequence of the mapped unitig or
* an empty string if there is no mapping (UnitigMap<T>::isEmpty = true).
*/
template<typename T>
string UnitigMap<T>::toString() const {

    if (isEmpty) return string();
    if (isShort) return cdbg->v_kmers[pos_unitig].first.toString();
    if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().toString();

    return cdbg->v_unitigs[pos_unitig]->seq.toString();
}

/** Return a string containing the reverse-complement sequence of the mapped unitig.
* @return a string containing the reverse-complement sequence of the mapped unitig or
* an empty string if there is no mapping (UnitigMap<T>::isEmpty = true).
*/
template<typename T>
string UnitigMap<T>::reverseToString() const {

    if (isEmpty) return string();
    if (isShort) return cdbg->v_kmers[pos_unitig].first.twin().toString();
    if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().twin().toString();

    return cdbg->v_unitigs[pos_unitig]->seq.rev().toString();
}

template<typename T>
size_t UnitigMap<T>::lcp(const char* s, const size_t pos_s, const size_t pos_um_seq, const bool um_reversed) const {

    if (isEmpty || (pos_s >= strlen(s))) return 0;

    if (isShort || isAbundant){

        if (pos_um_seq >= Kmer::k) return 0;

        char km_str[Kmer::k + 1];

        const Kmer km = isShort ? cdbg->v_kmers[pos_unitig].first : cdbg->h_kmers_ccov.find(pos_unitig).getKey();

        um_reversed ? km.twin().toString(km_str) : km.toString(km_str);

        return cstrMatch(&s[pos_s], &km_str[pos_um_seq]);
    }

    if (pos_um_seq >= cdbg->v_unitigs[pos_unitig]->length()) return 0;

    return cdbg->v_unitigs[pos_unitig]->seq.jump(s, pos_s, pos_um_seq, um_reversed);
}

/** Return the head k-mer of the mapped unitig.
* @return a Kmer object which is either the head k-mer of the mapped unitig or
* an empty k-mer if there is no mapping (UnitigMap<T>::isEmpty = true).
*/
template<typename T>
Kmer UnitigMap<T>::getHead() const {

    if (!isEmpty){

        if (isShort) return cdbg->v_kmers[pos_unitig].first;
        if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

        return cdbg->v_unitigs[pos_unitig]->seq.getKmer(0);
    }

    Kmer km;

    km.set_empty();

    return km;
}

/** Return the tail k-mer of the mapped unitig.
* @return a Kmer object which is either the tail k-mer of the mapped unitig or
* an empty k-mer if there is no mapping (UnitigMap<T>::isEmpty = true).
*/
template<typename T>
Kmer UnitigMap<T>::getTail() const {

    if (!isEmpty){

        if (isShort) return cdbg->v_kmers[pos_unitig].first;
        if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

        return cdbg->v_unitigs[pos_unitig]->seq.getKmer(cdbg->v_unitigs[pos_unitig]->numKmers() - 1);
    }

    Kmer km;

    km.set_empty();

    return km;
}

/** Return the k-mer starting at position pos in the mapped unitig.
* @param pos is the start position of the k-mer to extract.
* @return a Kmer object which is either the k-mer starting at position pos in the mapped unitig or
* an empty k-mer if there is either no mapping (UnitigMap<T>::isEmpty = true) or the position is
* greater than length(unitig)-k.
*/
template<typename T>
Kmer UnitigMap<T>::getKmer(const size_t pos) const {

    if (!isEmpty){

        if (isShort && (pos == 0)) return cdbg->v_kmers[pos_unitig].first;
        if (isAbundant && (pos == 0)) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

        if (!isShort && !isAbundant && (pos < cdbg->v_unitigs[pos_unitig]->numKmers())) {

            return cdbg->v_unitigs[pos_unitig]->seq.getKmer(pos);
        }
    }

    Kmer km;

    km.set_empty();

    return km;
}

/** Return a pointer to the data associated with the mapped unitig.
* @return a constant pointer to the data associated with the mapped unitig or
* a constant null pointer if there is no mapping (UnitigMap<T>::isEmpty = true) or no data
* associated with the unitigs (T = void).
*/
template<typename T>
const T* UnitigMap<T>::getData() const {

    if (isEmpty || !cdbg->has_data) return nullptr;
    if (isShort) return cdbg->v_kmers[pos_unitig].second.getData();
    if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig)->getData();

    return cdbg->v_unitigs[pos_unitig]->getData();
}

template<typename T>
T* UnitigMap<T>::getData() {

    if (isEmpty || !cdbg->has_data) return nullptr;
    if (isShort) return cdbg->v_kmers[pos_unitig].second.getData();
    if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig)->getData();

    return cdbg->v_unitigs[pos_unitig]->getData();
}

/** Set the data associated with a mapped unitig. The function does not set anything if
* there is no mapping (UnitigMap<T>::isEmpty = true) or no data associated with the unitigs (T = void).
* @param data is a pointer to the data that will be copied to the data associated with the mapped unitig.
*/
template<typename T>
void UnitigMap<T>::setData(const T* const data) const {

    if (isEmpty || !cdbg->has_data){

        // Raise error message with cerr
    }
    else if (isShort) cdbg->v_kmers[pos_unitig].second.setData(data);
    else if (isAbundant) cdbg->h_kmers_ccov.find(pos_unitig)->setData(data);
    else cdbg->v_unitigs[pos_unitig]->setData(data);
}

/** Merge the data of the mapped unitig with the data of another mapped unitig.
* This function calls the function CDBG_Data_t<T>::join but it does not set anything
* if there is no mapping (UnitigMap<T>::isEmpty = true) or no data associated with the unitigs (T = void).
* @param um is a reference to a mapped unitig that will be merged with the current mapped unitig.
*/
template<typename T>
void UnitigMap<T>::mergeData(const UnitigMap<T>& um) {

    getData()->join(*this, um);
}

template<> inline void UnitigMap<void>::mergeData(const UnitigMap<void>& um) { /* Does nothing */ }

/** Create new data to associate with a new unitig which is a subsequence of the current mapped unitig.
* This function calls the function CDBG_Data_t<T>::split but it does not create anything
* if there is no mapping (UnitigMap<T>::isEmpty = true) or no data associated with the unitigs (T = void).
* @param pos is the start position of the new unitig within the current mapped unitig.
* @param len is the length of the new unitig.
* @return a Unitig object containing the new data. The unitig sequence and coverage are not set, only
* the data (if there are some).
*/
template<typename T>
Unitig<T> UnitigMap<T>::splitData(const bool last_split) {

    Unitig<T> unitig;

    getData()->sub(*this, unitig.data, last_split);

    return unitig;
}

template<> inline Unitig<void> UnitigMap<void>::splitData(const bool last_split) { return Unitig<void>(); }

/** Return a BackwardCDBG object that can create iterators (through BackwardCDBG::begin() and
* BackwardCDBG::end()) over the predecessors of the current mapped unitig.
* @return a BackwardCDBG object that can create iterators (through BackwardCDBG::begin() and
* BackwardCDBG::end()) over the predecessors of the current mapped unitig.
*/
template<typename T>
BackwardCDBG<T, false> UnitigMap<T>::getPredecessors() { return BackwardCDBG<T, false>(*this); }

/** Return a BackwardCDBG object that can create iterators (through BackwardCDBG::begin() and
* BackwardCDBG::end()) over the successors of the current mapped unitig.
* @return a BackwardCDBG object that can create iterators (through BackwardCDBG::begin() and
* BackwardCDBG::end()) over the successors of the current mapped unitig.
*/
template<typename T>
ForwardCDBG<T, false> UnitigMap<T>::getSuccessors() { return ForwardCDBG<T, false>(*this); }

/** Return a constant BackwardCDBG object that can create iterators (through BackwardCDBG::begin() and
* BackwardCDBG::end()) over the predecessors of the current mapped unitig.
* @return a constant BackwardCDBG object that can create iterators (through BackwardCDBG::begin() and
* BackwardCDBG::end()) over the predecessors of the current mapped unitig.
*/
template<typename T>
BackwardCDBG<T, true> UnitigMap<T>::getPredecessors() const { return BackwardCDBG<T, true>(*this); }

/** Return a constant BackwardCDBG object that can create iterators (through BackwardCDBG::begin() and
* BackwardCDBG::end()) over the successors of the current mapped unitig.
* @return a constant BackwardCDBG object that can create iterators (through BackwardCDBG::begin() and
* BackwardCDBG::end()) over the successors of the current mapped unitig.
*/
template<typename T>
ForwardCDBG<T, true> UnitigMap<T>::getSuccessors() const { return ForwardCDBG<T, true>(*this); }

template<typename T>
typename UnitigMap<T>::neighbor_iterator UnitigMap<T>::bw_begin() {

    neighbor_iterator it(*this, false);
    it++;
    return it;
}

template<typename T>
typename UnitigMap<T>::const_neighbor_iterator UnitigMap<T>::bw_begin() const {

    const_neighbor_iterator it(*this, false);
    it++;
    return it;
}

template<typename T>
typename UnitigMap<T>::neighbor_iterator UnitigMap<T>::bw_end() { return neighbor_iterator(); }

template<typename T>
typename UnitigMap<T>::const_neighbor_iterator UnitigMap<T>::bw_end() const { return const_neighbor_iterator(); }

template<typename T>
typename UnitigMap<T>::neighbor_iterator UnitigMap<T>::fw_begin() {

    neighbor_iterator it(*this, true);
    it++;
    return it;
}

template<typename T>
typename UnitigMap<T>::const_neighbor_iterator UnitigMap<T>::fw_begin() const {

    const_neighbor_iterator it(*this, true);
    it++;
    return it;
}

template<typename T>
typename UnitigMap<T>::neighbor_iterator UnitigMap<T>::fw_end() { return neighbor_iterator(); }

template<typename T>
typename UnitigMap<T>::const_neighbor_iterator UnitigMap<T>::fw_end() const { return const_neighbor_iterator(); }
