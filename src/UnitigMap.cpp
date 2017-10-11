#include "UnitigMap.hpp"
#include "CompactedDBG.hpp"

UnitigMap::iterator UnitigMap::begin() {

    UnitigMap::iterator it(*this, dbg, true);
    it++;
    return it;
}

UnitigMap::const_iterator UnitigMap::begin() const {

    UnitigMap::const_iterator it(*this, dbg, true);
    it++;
    return it;
}

UnitigMap::iterator UnitigMap::rbegin() {

    UnitigMap::iterator it(*this, dbg, false);
    it++;
    return it;
}

UnitigMap::const_iterator UnitigMap::rbegin() const {

    UnitigMap::const_iterator it(*this, dbg, false);
    it++;
    return it;
}

UnitigMap::iterator UnitigMap::end() { return UnitigMap::iterator(); }
UnitigMap::const_iterator UnitigMap::end() const { return UnitigMap::const_iterator(); }
UnitigMap::iterator UnitigMap::rend() { return UnitigMap::iterator(); }
UnitigMap::const_iterator UnitigMap::rend() const { return UnitigMap::const_iterator(); }

string UnitigMap::toString() const {

    if (isEmpty) return string();
    else if (isShort) return dbg->v_kmers[pos_unitig].first.toString();
    else if (isAbundant) return dbg->h_kmers_ccov.find(pos_unitig)->first.toString();
    else return dbg->v_unitigs[pos_unitig]->seq.toString();
}

Kmer UnitigMap::getHead() const {

    if (!isEmpty){

        if (isShort) return dbg->v_kmers[pos_unitig].first;
        else if (isAbundant) return dbg->h_kmers_ccov.find(pos_unitig)->first;
        else return dbg->v_unitigs[pos_unitig]->seq.getKmer(dbg->v_unitigs[pos_unitig]->numKmers() - 1);
    }

    Kmer km;

    km.set_empty();

    return km;
}

Kmer UnitigMap::getTail() const {

    if (!isEmpty){

        if (isShort) return dbg->v_kmers[pos_unitig].first;
        else if (isAbundant) return dbg->h_kmers_ccov.find(pos_unitig)->first;
        else return dbg->v_unitigs[pos_unitig]->seq.getKmer(0);
    }

    Kmer km;

    km.set_empty();

    return km;
}

rUnitigMap UnitigMap::backward() const { return rUnitigMap(*this); }

UnitigMap::const_iterator rUnitigMap::begin() const { return um_.rbegin(); }
UnitigMap::const_iterator rUnitigMap::end() const { return um_.rend(); }

template <bool is_const>
neighborIterator<is_const>::neighborIterator() : i_(3), is_fw_(true), dbg_(NULL) {}

template <bool is_const>
neighborIterator<is_const>::neighborIterator(const Kmer km, const CompactedDBG* dbg, const bool is_forward) : i_(-1), is_fw_(is_forward), km_(km), dbg_(dbg) {

    if ((dbg_ == NULL) || dbg_->invalid) i_ = 3;
}

template <bool is_const>
neighborIterator<is_const>::neighborIterator(const UnitigMap& um, const CompactedDBG* dbg, const bool is_forward) : i_(-1), is_fw_(is_forward), dbg_(dbg) {

    if (um.isEmpty || (dbg_ == NULL) || dbg_->invalid) i_ = 3;
    else if (is_fw_) km_ = um.getTail();
    else km_ = um.getHead();
}

template <bool is_const>
neighborIterator<is_const>::neighborIterator(const neighborIterator& o) : i_(o.i_), is_fw_(o.is_fw_), um_(o.um_), km_(o.km_), dbg_(o.dbg_) {}

template <bool is_const>
neighborIterator<is_const>& neighborIterator<is_const>::operator++() {

    if ((dbg_ == NULL) || dbg_->invalid || (i_ >= 3)) return *this;

    while (i_ < 3){

        i_++;

        um_ = dbg_->find(is_fw_ ? km_.forwardBase(alpha[i_]) : km_.backwardBase(alpha[i_]));

        if (!um_.isEmpty) break;
    }

    return *this;
}

template <bool is_const>
neighborIterator<is_const> neighborIterator<is_const>::operator++(int) {

    neighborIterator tmp(*this);
    operator++();

    return tmp;
}

template <bool is_const>
bool neighborIterator<is_const>::operator==(const neighborIterator& o) {

    if ((i_ >= 3) || (o.i_ >= 3)) return (i_ >= 3) && (o.i_ >= 3);
    return (is_fw_ == o.is_fw_) && (km_ == o.km_) && (dbg_ == o.dbg_) && (um_ == o.um_);
}

template <bool is_const>
bool neighborIterator<is_const>::operator!=(const neighborIterator& o) { return !operator==(o); }

template <bool is_const>
typename neighborIterator<is_const>::UnitigMap_ref_t neighborIterator<is_const>::operator*() { return um_; }

template <bool is_const>
typename neighborIterator<is_const>::UnitigMap_ptr_t neighborIterator<is_const>::operator->() { return &um_; }
