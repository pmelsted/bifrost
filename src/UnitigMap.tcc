#ifndef BFG_UNITIGMAP_TCC
#define BFG_UNITIGMAP_TCC

#include "CompactedDBG.hpp"
#include "NeighborIterator.hpp"

template<typename U, typename G, bool is_const>
UnitigMap<U, G, is_const>::UnitigMap(const size_t length, CompactedDBG_ptr_t cdbg_) :   len(length), pos_unitig(0), dist(0), size(0), isShort(false),
                                                                                        isAbundant(false), isEmpty(true), strand(true), cdbg(cdbg_) {}

template<typename U, typename G, bool is_const>
UnitigMap<U, G, is_const>::UnitigMap(const size_t start, const size_t length, const size_t unitig_sz, const bool strand) :
                                    pos_unitig(0), dist(start), len(length), size(unitig_sz), cdbg(nullptr), strand(strand), isShort(false),
                                    isAbundant(false), isEmpty(false) {}

template<typename U, typename G, bool is_const>
bool UnitigMap<U, G, is_const>::operator==(const UnitigMap& o) const {

    return  (pos_unitig == o.pos_unitig) && (dist == o.dist) && (len == o.len) && (size == o.size) && (strand == o.strand) &&
            (isEmpty == o.isEmpty) && (isShort == o.isShort) && (isAbundant == o.isAbundant) && (cdbg == o.cdbg);
}

template<typename U, typename G, bool is_const>
bool UnitigMap<U, G, is_const>::operator!=(const UnitigMap& o) const {

    return  (pos_unitig != o.pos_unitig) || (dist != o.dist) || (len != o.len) || (size != o.size) || (strand != o.strand) ||
            (isEmpty != o.isEmpty) || (isShort != o.isShort) || (isAbundant != o.isAbundant) || (cdbg != o.cdbg);
}

template<typename U, typename G, bool is_const>
string UnitigMap<U, G, is_const>::toString() const {

    if (isEmpty) return string();

    if (strand){

        if (isShort) return cdbg->v_kmers[pos_unitig].first.toString();
        if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().toString();

        return cdbg->v_unitigs[pos_unitig]->seq.toString(dist, len + cdbg->k_ - 1);
    }
    else {

        if (isShort) return cdbg->v_kmers[pos_unitig].first.twin().toString();
        if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().twin().toString();

        return reverse_complement(cdbg->v_unitigs[pos_unitig]->seq.toString(dist, len + cdbg->k_ - 1));
    }
}

template<typename U, typename G, bool is_const>
size_t UnitigMap<U, G, is_const>::lcp(const char* s, const size_t pos_s, const size_t pos_um_seq, const bool um_reversed) const {

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

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getUnitigHead() const {

    if (!isEmpty){

        if (isShort) return cdbg->v_kmers[pos_unitig].first;
        if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

        return cdbg->v_unitigs[pos_unitig]->seq.getKmer(0);
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getUnitigTail() const {

    if (!isEmpty){

        if (isShort) return cdbg->v_kmers[pos_unitig].first;
        if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

        return cdbg->v_unitigs[pos_unitig]->seq.getKmer(cdbg->v_unitigs[pos_unitig]->numKmers() - 1);
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getUnitigKmer(const size_t pos) const {

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

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getMappedHead() const {

    if (!isEmpty){

        if (strand){

            if (isShort) return cdbg->v_kmers[pos_unitig].first;
            if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

            return cdbg->v_unitigs[pos_unitig]->seq.getKmer(dist);
        }
        else {

            if (isShort) return cdbg->v_kmers[pos_unitig].first.twin();
            if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().twin();

            return cdbg->v_unitigs[pos_unitig]->seq.getKmer(dist + len - 1).twin();
        }
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getMappedTail() const {

    if (!isEmpty){

        if (strand){

            if (isShort) return cdbg->v_kmers[pos_unitig].first;
            if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

            return cdbg->v_unitigs[pos_unitig]->seq.getKmer(dist + len - 1);
        }
        else {

            if (isShort) return cdbg->v_kmers[pos_unitig].first.twin();
            if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().twin();

            return cdbg->v_unitigs[pos_unitig]->seq.getKmer(dist).twin();
        }
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getMappedKmer(const size_t pos) const {

    if (!isEmpty && (pos < len)){

        if (strand){

            if (isShort && (pos + dist == 0)) return cdbg->v_kmers[pos_unitig].first;
            if (isAbundant && (pos + dist == 0)) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

            if (!isShort && !isAbundant && (pos < cdbg->v_unitigs[pos_unitig]->numKmers())) {

                return cdbg->v_unitigs[pos_unitig]->seq.getKmer(pos);
            }
        }
        else {

            if (isShort && (pos + dist == 0)) return cdbg->v_kmers[pos_unitig].first.twin();
            if (isAbundant && (pos + dist == 0)) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().twin();

            if (!isShort && !isAbundant && (pos < cdbg->v_unitigs[pos_unitig]->numKmers())) {

                return cdbg->v_unitigs[pos_unitig]->seq.getKmer(pos).twin();
            }
        }
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename U, typename G, bool is_const>
UnitigMap<U, G, is_const> UnitigMap<U, G, is_const>::getKmerMapping(const size_t pos) const {

    UnitigMap<U, G, is_const> um(*this);

    um.dist = pos;
    um.len = 1;

    if (pos >= um.size) um.isEmpty = true;

    return um;
}

template<typename U, typename G, bool is_const>
typename UnitigMap<U, G, is_const>::Unitig_data_ptr_t UnitigMap<U, G, is_const>::getData() const {

    return getData_<is_void<U>::value>();
}

template<typename U, typename G, bool is_const>
typename UnitigMap<U, G, is_const>::UnitigMap_BW UnitigMap<U, G, is_const>::getPredecessors() const {

    return BackwardCDBG<U, G, is_const>(*this);
}

template<typename U, typename G, bool is_const>
typename UnitigMap<U, G, is_const>::UnitigMap_FW UnitigMap<U, G, is_const>::getSuccessors() const {

    return ForwardCDBG<U, G, is_const>(*this);
}

template<typename U, typename G, bool is_const>
template<bool is_void>
typename std::enable_if<!is_void, typename UnitigMap<U, G, is_const>::Unitig_data_ptr_t>::type UnitigMap<U, G, is_const>::getData_() const {

    if (isEmpty) return nullptr;
    if (isShort) return cdbg->v_kmers[pos_unitig].second.getData();
    if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig)->getData();

    return cdbg->v_unitigs[pos_unitig]->getData();
}

template<typename U, typename G, bool is_const>
template<bool is_void>
typename std::enable_if<is_void, typename UnitigMap<U, G, is_const>::Unitig_data_ptr_t>::type UnitigMap<U, G, is_const>::getData_() const {

    return nullptr;
}

template<typename U, typename G, bool is_const>
template<bool is_void>
typename std::enable_if<!is_void, void>::type UnitigMap<U, G, is_const>::mergeData_(const UnitigMap<U, G, is_const>& um) const {

    U::join(*this, um);
}

template<typename U, typename G, bool is_const>
template<bool is_void>
typename std::enable_if<is_void, void>::type UnitigMap<U, G, is_const>::mergeData_(const UnitigMap<U, G, is_const>& um) const {}

template<typename U, typename G, bool is_const>
void UnitigMap<U, G, is_const>::mergeData(const UnitigMap<U, G, is_const>& um) const {

    mergeData_<is_void<U>::value>(um);
}


template<typename U, typename G, bool is_const>
template<bool is_void>
typename std::enable_if<!is_void, Unitig<U>>::type UnitigMap<U, G, is_const>::splitData_(const bool last_split) const {

    Unitig<U> unitig;

    U::sub(*this, &(unitig.data), last_split);

    return unitig;
}

template<typename U, typename G, bool is_const>
template<bool is_void>
typename std::enable_if<is_void, Unitig<U>>::type UnitigMap<U, G, is_const>::splitData_(const bool last_split) const {

    return Unitig<U>();
}

template<typename U, typename G, bool is_const>
Unitig<U> UnitigMap<U, G, is_const>::splitData(const bool last_split) const {

    return splitData_<is_void<U>::value>(last_split);
}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> UnitigMap<U, G, is_const>::bw_begin() const {

    neighborIterator<U, G, is_const> it(*this, false);
    it++;
    return it;
}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> UnitigMap<U, G, is_const>::bw_end() const { return neighborIterator<U, G, is_const>(); }

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> UnitigMap<U, G, is_const>::fw_begin() const {

    neighborIterator<U, G, is_const> it(*this, true);
    it++;
    return it;
}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> UnitigMap<U, G, is_const>::fw_end() const { return neighborIterator<U, G, is_const>(); }

template<typename U, typename G, bool is_const>
void UnitigMap<U, G, is_const>::partialCopy(const UnitigMap<U, G, is_const>& o) {

    pos_unitig = o.pos_unitig;
    dist = o.dist;
    len = o.len;
    size = o.size;

    strand = o.strand;
    isEmpty = o.isEmpty;

    isShort = o.isShort;
    isAbundant = o.isAbundant;
}

template<typename U, typename G, bool is_const>
UnitigMap<U, G, is_const>::UnitigMap(size_t p_unitig, size_t i, size_t l, size_t sz, bool short_, bool abundance, bool strd, CompactedDBG_ptr_t cdbg_) :
                                    pos_unitig(p_unitig), dist(i), len(l), size(sz), cdbg(cdbg_), strand(strd), isShort(short_), isAbundant(abundance),
                                    isEmpty(false) {}

#endif
