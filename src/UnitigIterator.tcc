#include "CompactedDBG.hpp"

/** Constructor.
* @return an empty unitigIterator.
*/
template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>::unitigIterator() :  i(0), v_unitigs_sz(0), v_kmers_sz(0), h_kmers_ccov_sz(0), sz(0), invalid(true), cdbg(nullptr) {}

/** Constructor.
* @param cdbg_ is a pointer to a Compacted de Bruijn graph. The unitigIterator will iterate over the unitigs of this graph.
* If the unitigIterator is a constant, the CompactedDBG pointer is a constant pointer and the accessed unitigs can't be modified.
* @return a unitigIterator.
*/
template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>::unitigIterator(CompactedDBG_ptr_t cdbg_) :
                i(0), v_unitigs_sz(0), v_kmers_sz(0), h_kmers_ccov_sz(0), sz(0), invalid(true), cdbg(cdbg_),
                it_h_kmers_ccov((cdbg_ == nullptr) || cdbg_->invalid ? typename KmerHashTable<CompressedCoverage_t<U>>::const_iterator() : cdbg_->h_kmers_ccov.begin()){

    if ((cdbg != nullptr) && !cdbg->invalid && (cdbg->size() != 0)){

        invalid = false;

        v_unitigs_sz = cdbg->v_unitigs.size();
        v_kmers_sz = cdbg->v_kmers.size();
        h_kmers_ccov_sz = cdbg->h_kmers_ccov.size();

        sz = v_unitigs_sz + v_kmers_sz + h_kmers_ccov_sz;
    }
}

/** Copy constructor.
* @return a copy of a unitigIterator.
*/
template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>::unitigIterator(const unitigIterator& o) :   i(o.i), v_unitigs_sz(o.v_unitigs_sz), v_kmers_sz(o.v_kmers_sz),
                                                                            it_h_kmers_ccov(o.it_h_kmers_ccov), h_kmers_ccov_sz(o.h_kmers_ccov_sz),
                                                                            sz(o.sz), invalid(o.invalid), um(o.um), cdbg(o.cdbg) {}

/** Prefix increment, iterate over the next unitig.
* Note that no specific order (such as a lexicographic one) is assumed during iteration.
*/
template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>& unitigIterator<U, G, is_const>::operator++() {

    if (invalid) return *this;

    if ((cdbg == nullptr) || cdbg->invalid || (i >= sz)){

        invalid = true;
        return *this;
    }

    if (i < v_unitigs_sz){

        um = UnitigMap<U, G, is_const>(i, 0, cdbg->v_unitigs[i]->seq.size() - cdbg->getK() + 1,
                                       cdbg->v_unitigs[i]->seq.size(), false, false, true, cdbg);
    }
    else if (i < (v_unitigs_sz + v_kmers_sz)){

        um = UnitigMap<U, G, is_const>(i - v_unitigs_sz, 0, 1, cdbg->getK(), true, false, true, cdbg);
    }
    else {

        um = UnitigMap<U, G, is_const>(it_h_kmers_ccov.getHash(), 0, 1, cdbg->getK(), false, true, true, cdbg);

        ++it_h_kmers_ccov;
    }

    ++i;

    return *this;
}

/** Postfix increment, iterate over the next unitig.
* Note that no specific order (such as a lexicographic one) is assumed during iteration.
*/
template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const> unitigIterator<U, G, is_const>::operator++(int) {

    unitigIterator<U, G, is_const> tmp(*this);
    operator++();

    return tmp;
}

/** Check if two unitigIterator are the same.
* @param o is another unitigIterator.
* @return a boolean indicating whether the two unitigIterator are the same.
*/
template<typename U, typename G, bool is_const>
bool unitigIterator<U, G, is_const>::operator==(const unitigIterator& o) const {

    if (invalid || o.invalid) return invalid && o.invalid;
    return  (i == o.i) && (v_unitigs_sz == o.v_unitigs_sz) && (v_kmers_sz == o.v_kmers_sz) &&
            (h_kmers_ccov_sz == o.h_kmers_ccov_sz) && (sz == o.sz) && (it_h_kmers_ccov == o.it_h_kmers_ccov) &&
            (cdbg == o.cdbg) && (um == o.um);
}

/** Check if two unitigIterator are different.
* @param o is another unitigIterator.
* @return a boolean indicating whether the two unitigIterator are different.
*/
template<typename U, typename G, bool is_const>
bool unitigIterator<U, G, is_const>::operator!=(const unitigIterator& o) const { return !operator==(o); }

/** Return a UnitigMap object reference which contains information about the mapped unitig.
* @return a UnitigMap object reference which contains information about the mapped unitig.
*/
template<typename U, typename G, bool is_const>
const UnitigMap<U, G, is_const>& unitigIterator<U, G, is_const>::operator*() const { return um; }

/** Return a UnitigMap object pointer which contains information about the mapped unitig.
* @return a UnitigMap object pointer which contains information about the mapped unitig.
*/
template<typename U, typename G, bool is_const>
const UnitigMap<U, G, is_const>* unitigIterator<U, G, is_const>::operator->() const { return &um; }
