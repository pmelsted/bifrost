#include "CompactedDBG.hpp"
#include "NeighborIterator.hpp"

template<typename T>
UnitigMap<T>::UnitigMap(size_t p_unitig, size_t i, size_t l, size_t sz, bool short_, bool abundance, bool strd, CompactedDBG<T>& cdbg_) :
                        pos_unitig(p_unitig), dist(i), len(l), size(sz), cdbg(&cdbg_), strand(strd), isShort(short_), isAbundant(abundance),
                        selfLoop(false), isEmpty(false), isTip(false), isIsolated(false) {}

template<typename T>
UnitigMap<T>::UnitigMap(size_t l) : len(l), isTip(false), isIsolated(false), isShort(false), isAbundant(false), isEmpty(true),
                                    selfLoop(false), strand(true), pos_unitig(0), dist(0), size(0), cdbg(nullptr) {}
template<typename T>
bool UnitigMap<T>::operator==(const UnitigMap& o) {

    return  (pos_unitig == o.pos_unitig) && (dist == o.dist) && (len == o.len) && (size == o.size) && (strand == o.strand) &&
            (selfLoop == o.selfLoop) && (isEmpty == o.isEmpty) && (isShort == o.isShort) && (isAbundant == o.isAbundant) &&
            (isIsolated == o.isIsolated) && (isTip == o.isTip) && (cdbg == o.cdbg);
}

template<typename T> bool UnitigMap<T>::operator!=(const UnitigMap& o) { return !operator==(o); }

template<typename T>
string UnitigMap<T>::toString() const {

    if (isEmpty) return string();
    else if (isShort) return cdbg->v_kmers[pos_unitig].first.toString();
    else if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig)->first.toString();
    else return cdbg->v_unitigs[pos_unitig]->seq.toString();
}

template<typename T>
Kmer UnitigMap<T>::getHead() const {

    if (!isEmpty){

        if (isShort) return cdbg->v_kmers[pos_unitig].first;
        else if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig)->first;
        else return cdbg->v_unitigs[pos_unitig]->seq.getKmer(0);
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename T>
Kmer UnitigMap<T>::getTail() const {

    if (!isEmpty){

        if (isShort) return cdbg->v_kmers[pos_unitig].first;
        else if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig)->first;
        else return cdbg->v_unitigs[pos_unitig]->seq.getKmer(cdbg->v_unitigs[pos_unitig]->numKmers() - 1);
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename T>
const T* UnitigMap<T>::getData() const {

    if (isEmpty || !cdbg->has_data) return nullptr;
    else if (isShort) return cdbg->v_kmers[pos_unitig].second.getData();
    else if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig)->second.getData();
    else return cdbg->v_unitigs[pos_unitig]->getData();
}

template<typename T>
T* UnitigMap<T>::getData() {

    if (isEmpty || !cdbg->has_data) return nullptr;
    else if (isShort) return cdbg->v_kmers[pos_unitig].second.getData();
    else if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig)->second.getData();
    else return cdbg->v_unitigs[pos_unitig]->getData();
}

template<typename T>
void UnitigMap<T>::setData(const T* const data) {

    if (isEmpty || !cdbg->has_data){

        // Raise error message with cerr
    }
    else if (isShort) cdbg->v_kmers[pos_unitig].second.setData(data);
    else if (isAbundant) cdbg->h_kmers_ccov.find(pos_unitig)->second.setData(data);
    else cdbg->v_unitigs[pos_unitig]->setData(data);
}

template<typename T>
void UnitigMap<T>::mergeData(const UnitigMap& um) {

    getData()->join(*um.getData(), *cdbg);
}

template<> void UnitigMap<void>::mergeData(const UnitigMap& um) { /* Does nothing */ }

template<typename T>
Unitig<T> UnitigMap<T>::splitData(const size_t pos, const size_t len) {

    Unitig<T> unitig;

    getData()->split(pos, len, unitig.data, *cdbg);

    return unitig;
}

template<> Unitig<void> UnitigMap<void>::splitData(const size_t pos, const size_t len) { return Unitig<void>(); }

template<typename T>
BackwardCDBG<T, false> UnitigMap<T>::getPredecessors() { return BackwardCDBG<T, false>(*this); }

template<typename T>
ForwardCDBG<T, false> UnitigMap<T>::getSuccessors() { return ForwardCDBG<T, false>(*this); }

template<typename T>
BackwardCDBG<T, true> UnitigMap<T>::getPredecessors() const { return BackwardCDBG<T, true>(*this); }

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
