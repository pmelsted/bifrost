#include "CompactedDBG.hpp"

template<typename T, bool is_const>
unitigIterator<T,is_const>::unitigIterator() :  i(0), v_unitigs_sz(0), v_kmers_sz(0), h_kmers_ccov_sz(0), sz(0), invalid(true), cdbg(NULL) {}

template<typename T, bool is_const>
unitigIterator<T,is_const>::unitigIterator(CompactedDBG<T>* cdbg_) :
                i(0), v_unitigs_sz(0), v_kmers_sz(0), h_kmers_ccov_sz(0), sz(0), invalid(true), cdbg(cdbg_),
                it_h_kmers_ccov((cdbg_ == NULL) || cdbg_->invalid ? typename KmerHashTable<CompressedCoverage_t<T>>::const_iterator() : cdbg_->h_kmers_ccov.begin()){

    if ((cdbg != NULL) && !cdbg->invalid && (cdbg->size() != 0)){

        invalid = false;

        v_unitigs_sz = cdbg->v_unitigs.size();
        v_kmers_sz = cdbg->v_kmers.size();
        h_kmers_ccov_sz = cdbg->h_kmers_ccov.size();

        sz = v_unitigs_sz + v_kmers_sz + h_kmers_ccov_sz;
    }
}

template<typename T, bool is_const>
unitigIterator<T,is_const>::unitigIterator(const unitigIterator& o) :   i(o.i), v_unitigs_sz(o.v_unitigs_sz), v_kmers_sz(o.v_kmers_sz),
                                            h_kmers_ccov_sz(o.h_kmers_ccov_sz), sz(o.sz), invalid(o.invalid),
                                            um(o.um), cdbg(o.cdbg) {}
template<typename T, bool is_const>
unitigIterator<T,is_const>& unitigIterator<T,is_const>::operator++() {

    if (invalid) return *this;

    if ((cdbg == NULL) || cdbg->invalid || (i >= sz)){

        invalid = true;
        return *this;
    }

    if (i < v_unitigs_sz) um = UnitigMap<T>(i, 0, 1, cdbg->v_unitigs[i]->seq.size(), false, false, true, *cdbg);
    else if (i < (v_unitigs_sz + v_kmers_sz)){

        um = UnitigMap<T>(i - v_unitigs_sz, 0, 1, cdbg->getK(), true, false, true, *cdbg);
    }
    else {

        um = UnitigMap<T>(it_h_kmers_ccov.getHash(), 0, 1, cdbg->getK(), false, true, true, *cdbg);

        it_h_kmers_ccov++;
    }

    i++;

    return *this;
}

template<typename T, bool is_const>
unitigIterator<T,is_const> unitigIterator<T,is_const>::operator++(int) {

    unitigIterator<T,is_const> tmp(*this);
    operator++();

    return tmp;
}

template<typename T, bool is_const>
bool unitigIterator<T,is_const>::operator==(const unitigIterator& o) {

    if (invalid || o.invalid) return invalid && o.invalid;
    return  (i == o.i) && (v_unitigs_sz == o.v_unitigs_sz) && (v_kmers_sz == o.v_kmers_sz) &&
            (h_kmers_ccov_sz == o.h_kmers_ccov_sz) && (sz == o.sz) && (it_h_kmers_ccov == o.it_h_kmers_ccov) &&
            (cdbg == o.cdbg) && (um == o.um);
}

template<typename T, bool is_const>
bool unitigIterator<T,is_const>::operator!=(const unitigIterator& o) { return !operator==(o); }

template<typename T, bool is_const>
typename unitigIterator<T,is_const>::UnitigMap_ref_t unitigIterator<T,is_const>::operator*() { return um; }

template<typename T, bool is_const>
typename unitigIterator<T,is_const>::UnitigMap_ptr_t unitigIterator<T,is_const>::operator->() { return &um; }
