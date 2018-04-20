#ifndef BFG_HASHID_TCC
#define BFG_HASHID_TCC

template<typename U>
DataAccessor<U>::DataAccessor(const uint8_t id) : da_id(id) {}

template<typename U>
const U* DataAccessor<U>::getData(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getCompactedDBG() != nullptr)){

        const DataStorage<U>* ds = um.getCompactedDBG()->getData();

        if (ds != nullptr) return ds->getData(um);
    }

    return nullptr;
}

template<> inline const void* DataAccessor<void>::getData(const const_UnitigColorMap<void>& um) const { return nullptr; }

template<typename U>
U* DataAccessor<U>::getData(const UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getCompactedDBG() != nullptr)){

        DataStorage<U>* ds = um.getCompactedDBG()->getData();

        if (ds != nullptr) return ds->getData(um);
    }

    return nullptr;
}

template<> inline void* DataAccessor<void>::getData(const UnitigColorMap<void>& um) const { return nullptr; }

template<typename U>
const UnitigColors* DataAccessor<U>::getUnitigColors(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getCompactedDBG() != nullptr)){

        const DataStorage<U>* ds = um.getCompactedDBG()->getData();

        if (ds != nullptr) return ds->getUnitigColors(um);
    }

    return nullptr;
}

template<typename U>
UnitigColors* DataAccessor<U>::getUnitigColors(const UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getCompactedDBG() != nullptr)){

        DataStorage<U>* ds = um.getCompactedDBG()->getData();

        if (ds != nullptr) return ds->getUnitigColors(um);
    }

    return nullptr;
}

template<typename U>
UnitigColors DataAccessor<U>::getSubUnitigColors(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getCompactedDBG() != nullptr)){

        const DataStorage<U>* ds = um.getCompactedDBG()->getData();

        if (ds != nullptr) return ds->getSubUnitigColors(um);
    }

    return UnitigColors();
}

template<typename U>
vector<string> DataAccessor<U>::getSubUnitigColorNames(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getCompactedDBG() != nullptr)){

        const DataStorage<U>* ds = um.getCompactedDBG()->getData();

        if (ds != nullptr) return ds->getSubUnitigColorNames(um);
    }

    return vector<string>();
}

template<typename U>
void DataAccessor<U>::join(const UnitigColorMap<U>& um_dest, const UnitigColorMap<U>& um_src){

    DataStorage<U>* ds = um_dest.getCompactedDBG()->getData();

    DataAccessor<U>* da_dest = um_dest.getData();
    UnitigColors* cs_dest = da_dest->getUnitigColors(um_dest);

    const DataAccessor<U>* da_src = um_src.getData();
    UnitigColors* cs_src = da_src->getUnitigColors(um_src);

    U* data_unitig_dest = da_dest->getData(um_dest);
    U* data_unitig_src = da_src->getData(um_src);

    if (cs_dest != nullptr){ // If a colorset exists for um_dest

        const Kmer head = um_dest.getUnitigHead();
        const Kmer new_head = um_dest.strand ? head : um_dest.getUnitigTail().twin();

        ds->joinUnitigColors(um_dest, um_src); // Join the color sets

        // TODO: Insert in tombstone if available

        if (new_head != head){

            // Insert new colorset with corresponding head into overflow of k-mers
            ds->overflow.insert(new_head, cs_dest - ds->color_sets);

            *da_dest = DataAccessor<U>();
        }
    }

    if (cs_src != nullptr){

        const uint64_t h_src = ds->getHash(um_src) % ds->nb_color_sets;

        ds->unitig_cs_link[h_src >> 6] &= ~(1ULL << (h_src & 0x3F));
    }

    if ((data_unitig_dest != nullptr) && (data_unitig_src != nullptr)){

        U::join(um_dest, um_src);
    }
}

template<>
inline void DataAccessor<void>::join(const UnitigColorMap<void>& um_dest, const UnitigColorMap<void>& um_src){

    DataStorage<void>* ds = um_dest.getCompactedDBG()->getData();

    DataAccessor<void>* da_dest = um_dest.getData();
    UnitigColors* cs_dest = da_dest->getUnitigColors(um_dest);

    const DataAccessor<void>* da_src = um_src.getData();
    UnitigColors* cs_src = da_src->getUnitigColors(um_src);

    if (cs_dest != nullptr){ // If a colorset exists for um_dest

        const Kmer head = um_dest.getUnitigHead();
        const Kmer new_head = um_dest.strand ? head : um_dest.getUnitigTail().twin();

        ds->joinUnitigColors(um_dest, um_src); // Join the color sets

        // TODO: Insert in tombstone if available

        if (new_head != head){

            // Insert new colorset with corresponding head into overflow of k-mers
            ds->overflow.insert(new_head, cs_dest - ds->color_sets);

            *da_dest = DataAccessor<void>();
        }
    }

    if (cs_src != nullptr){

        const uint64_t h_src = ds->getHash(um_src) % ds->nb_color_sets;

        ds->unitig_cs_link[h_src >> 6] &= ~(1ULL << (h_src & 0x3F));
    }
}

template<typename U>
void DataAccessor<U>::sub(DataAccessor<U>* data_dest, const UnitigColorMap<U>& um_src, const bool last_extraction) {

    DataStorage<U>* ds = um_src.getCompactedDBG()->getData();
    UnitigColors cs = ds->getSubUnitigColors(um_src);

    if (cs.size() != 0){

        const Kmer km = um_src.getMappedHead();

        UnitigColors* cs_um_src = (last_extraction ? ds->getUnitigColors(um_src) : ds->insert());

        *cs_um_src = move(cs);

        ds->overflow.insert(km, cs_um_src - ds->color_sets);
    }

    U* data_unitig = ds->getData(um_src);

    if (data_unitig != nullptr){

        const UnitigMapBase umb(0, um_src.len, um_src.len + Kmer::k - 1, true);

        U::sub(data_unitig, cs, umb, um_src, last_extraction);
    }
}

template<>
inline void DataAccessor<void>::sub(DataAccessor<void>* data_dest, const UnitigColorMap<void>& um_src, const bool last_extraction) {

    DataStorage<void>* ds = um_src.getCompactedDBG()->getData();
    UnitigColors cs = ds->getSubUnitigColors(um_src);

    if (cs.size() != 0){

        const Kmer km = um_src.getMappedHead();

        UnitigColors* cs_um_src = (last_extraction ? ds->getUnitigColors(um_src) : ds->insert());

        *cs_um_src = move(cs);

        ds->overflow.insert(km, cs_um_src - ds->color_sets);
    }
}

template<typename U>
string DataAccessor<U>::serialize() const { return std::to_string(da_id); }

#endif
