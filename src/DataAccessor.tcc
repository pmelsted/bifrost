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

template<> const void* DataAccessor<void>::getData(const const_UnitigColorMap<void>& um) const { return nullptr; }

template<typename U>
U* DataAccessor<U>::getData(const UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getCompactedDBG() != nullptr)){

        DataStorage<U>* ds = um.getCompactedDBG()->getData();

        if (ds != nullptr) return ds->getData(um);
    }

    return nullptr;
}

template<> void* DataAccessor<void>::getData(const UnitigColorMap<void>& um) const { return nullptr; }

/** Get the color set of a unitig.
* @param um is a UnitigMap representing a mapping to a unitig from the graph.
* @return a constant pointer to the color set matching um. If no such color set
* is found, the pointer is nullptr.
*/
template<typename U>
const UnitigColors<U>* DataAccessor<U>::getUnitigColors(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getCompactedDBG() != nullptr)){

        const DataStorage<U>* ds = um.getCompactedDBG()->getData();

        if (ds != nullptr) return ds->getUnitigColors(um);
    }

    return nullptr;
}

/** Get the color set of a unitig.
* @param um is a UnitigMap representing a mapping to a unitig from the graph.
* @return a constant pointer to the color set matching um. If no such color set
* is found, the pointer is nullptr.
*/
template<typename U>
UnitigColors<U>* DataAccessor<U>::getUnitigColors(const UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getCompactedDBG() != nullptr)){

        DataStorage<U>* ds = um.getCompactedDBG()->getData();

        if (ds != nullptr) return ds->getUnitigColors(um);
    }

    return nullptr;
}

/** Extract the color set matching a sub-unitig (see UnitigMap).
* @param um is a UnitigMap representing a mapping to a unitig from the graph.
* @return a new color set
*/
template<typename U>
UnitigColors<U> DataAccessor<U>::getSubUnitigColors(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getCompactedDBG() != nullptr)){

        const DataStorage<U>* ds = um.getCompactedDBG()->getData();

        if (ds != nullptr) return ds->getSubUnitigColors(um);
    }

    return UnitigColors<U>();
}

/** Same as ColoredCDBG::extractColors but extract the color names instead.
* @param um is a UnitigMap representing a mapping to a unitig from the graph.
* @return a vector a string. Each string is the name of a color.
*/
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
    UnitigColors<U>* cs_dest = da_dest->getUnitigColors(um_dest);

    const DataAccessor<U>* da_src = um_src.getData();
    UnitigColors<U>* cs_src = da_src->getUnitigColors(um_src);

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

    if (cs_src != nullptr) cs_src->setUnoccupied();

    if ((data_unitig_dest != nullptr) && (data_unitig_src != nullptr)){

        U::join(um_dest, um_src);
    }
}

template<>
void DataAccessor<void>::join(const UnitigColorMap<void>& um_dest, const UnitigColorMap<void>& um_src){

    DataAccessor<void>* da_dest = um_dest.getData();
    UnitigColors<void>* cs_dest = da_dest->getUnitigColors(um_dest);

    const DataAccessor<void>* da_src = um_src.getData();
    UnitigColors<void>* cs_src = da_src->getUnitigColors(um_src);

    if (cs_dest != nullptr){ // If a colorset exists for um_dest

        const Kmer head = um_dest.getUnitigHead();
        const Kmer new_head = um_dest.strand ? head : um_dest.getUnitigTail().twin();

        DataStorage<void>* ds = um_dest.getCompactedDBG()->getData();

        ds->joinUnitigColors(um_dest, um_src); // Join the color sets

        // TODO: Insert in tombstone if available

        if (new_head != head){

            // Insert new colorset with corresponding head into overflow of k-mers
            ds->overflow.insert(new_head, cs_dest - ds->color_sets);

            *da_dest = DataAccessor<void>();
        }
    }

    if (cs_src != nullptr) cs_src->setUnoccupied();
}

template<typename U>
void DataAccessor<U>::sub(const UnitigColorMap<U>& um, DataAccessor<U>* data_dest, const bool last_extraction) {

    DataStorage<U>* ds = um.getCompactedDBG()->getData();
    UnitigColors<U> cs = ds->getSubUnitigColors(um);

    if (cs.size() != 0){

        //const Kmer km = um.getKmer(um.dist);
        const Kmer km = um.getMappedHead();

        UnitigColors<U>* cs_um = (last_extraction ? ds->getUnitigColors(um) : ds->insert());

        *cs_um = move(cs);

        ds->overflow.insert(km, cs_um - ds->color_sets);
    }

    U* data_unitig = ds->getData(um);

    if (data_unitig != nullptr) U::sub(um, data_unitig, last_extraction);
}

template<>
void DataAccessor<void>::sub(const UnitigColorMap<void>& um, DataAccessor<void>* data_dest, const bool last_extraction) {

    DataStorage<void>* ds = um.getCompactedDBG()->getData();
    UnitigColors<void> cs = ds->getSubUnitigColors(um);

    if (cs.size() != 0){

        //const Kmer km = um.getKmer(um.dist);
        const Kmer km = um.getMappedHead();

        UnitigColors<void>* cs_um = (last_extraction ? ds->getUnitigColors(um) : ds->insert());

        *cs_um = move(cs);

        ds->overflow.insert(km, cs_um - ds->color_sets);
    }
}

template<typename U>
string DataAccessor<U>::serialize() const { return std::to_string(da_id); }
