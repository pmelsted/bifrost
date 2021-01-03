#ifndef BIFROST_DATA_ACCESSOR_TCC
#define BIFROST_DATA_ACCESSOR_TCC

template<typename U>
void DataAccessor<U>::clear(const UnitigColorMap<U>& um) {

    if (!um.isEmpty){

        dac.cs.clear();
        dac.data.clear(um);
    }
}

template<>
inline void DataAccessor<void>::clear(const UnitigColorMap<void>& um) {

    if (!um.isEmpty) dac.cs.clear();
}

template<typename U>
const U* DataAccessor<U>::getData(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty) return &data;

    return nullptr;
}

template<>
inline const void* DataAccessor<void>::getData(const const_UnitigColorMap<void>& um) const {

    return nullptr;
}

template<typename U>
U* DataAccessor<U>::getData(const UnitigColorMap<U>& um) {

    if (!um.isEmpty) return &data;

    return nullptr;
}

template<>
inline void* DataAccessor<void>::getData(const UnitigColorMap<void>& um) {

    return nullptr;
}

template<typename U>
const UnitigColors& DataAccessor<U>::getUnitigColors(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty) return cs;

    return nullptr;
}

template<typename U>
UnitigColors& DataAccessor<U>::getUnitigColors(const UnitigColorMap<U>& um) {

    if (!um.isEmpty) return cs;

    return nullptr;
}

template<typename U>
UnitigColors DataAccessor<U>::getSubUnitigColors(const const_UnitigColorMap<U>& um) const {

    UnitigColors new_cs;

    if (!um.isEmpty && (um.getGraph() != nullptr)){

        const UnitigColors& cs = getUnitigColors(um);
        const size_t nb_colors = um.getGraph()->getData()->color_names.size();

        UnitigColorMap<U> um_tmp(0, 1, um.len + um.getGraph()->getK() - 1, um.strand);

        if ((um.len * nb_colors * 16) < cs->size(um)){

            const size_t end = um.dist + um.len;
            const size_t um_km_sz = um.size - um.getGraph()->getK() + 1;

            for (size_t colorID = 0; colorID < nb_colors; ++colorID){

                for (size_t km_dist = um.dist; km_dist < end; ++km_dist){

                    if (cs->contains(colorID * um_km_sz + km_dist)){

                        um_tmp.dist = um.strand ? (km_dist - um.dist) : (um.dist + um.len - km_dist - 1);
                        new_cs.add(um_tmp, colorID);
                    }
                }
            }
        }
        else {

            UnitigColors::const_iterator it(cs->begin(um));
            const UnitigColors::const_iterator it_end(cs->end());

            while (it != it_end){

                const size_t km_dist = it.getKmerPosition();

                um_tmp.dist = um.strand ? (km_dist - um.dist) : (um.dist + um.len - km_dist - 1);
                new_cs.add(um_tmp, it.getColorID());

                ++it;
            }
        }
    }

    return new_cs;
}

template<typename U>
vector<string> DataAccessor<U>::getSubUnitigColorNames(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getGraph() != nullptr)){

        const DataStorage<U>* ds = um.getGraph()->getData();

        if (ds != nullptr) return ds->getSubUnitigColorNames(um);
    }

    return vector<string>();
}

template<typename U>
void DataAccessor<U>::concat(const UnitigColorMap<U>& um_a, const UnitigColorMap<U>& um_b){

    clear();

    dac.uc = concatUnitigColors(um_a, um_b);

    if ((um_a.getData()->getData(um_a) != nullptr) || (um_b.getData()->getData(um_b) != nullptr)) dac.data.concat(um_a, um_b);
}

template<>
inline void DataAccessor<void>::concat(const UnitigColorMap<void>& um_a, const UnitigColorMap<void>& um_b){

    clear();

    dac.uc = concatUnitigColors(um_a, um_b);
}

template<typename U>
void DataAccessor<U>::merge(const UnitigColorMap<U>& um_a, const const_UnitigColorMap<U>& um_b){

    mergeUnitigColors(um_a, um_b);

    um_a.getData()->getData(um_a)->merge(um_a, um_b);
}

template<>
inline void DataAccessor<void>::merge(const UnitigColorMap<void>& um_dest, const const_UnitigColorMap<void>& um_src){

    mergeUnitigColors(um_a, um_b);
}

template<typename U>
void DataAccessor<U>::extract(const UnitigColorMap<U>& um_a, const bool last_extraction) {

    clear();

    dac.uc = um_a.getData()->getSubUnitigColors(um_a);
    dac.data = um_a.getData()->getData()->extract(um_a, last_extraction);
}

template<>
inline void DataAccessor<void>::extract(const UnitigColorMap<void>& um_a, const bool last_extraction) {

    clear();

    dac.uc = um_a.getData()->getSubUnitigColors(um_a);
}

template<typename U>
UnitigColors DataAccessor<U>::concatUnitigColors(const const_UnitigColorMap<U>& um_a, const const_UnitigColorMap<U>& um_b) const {

    UnitigColors new_cs;

    if (!um_a.isEmpty && !um_b.isEmpty && (um_a.getGraph() == um_b.getGraph())){

        const UnitigColors* color_set_a = um_a.getData()->getUnitigColors(um_a);
        const UnitigColors* color_set_b = um_b.getData()->getUnitigColors(um_b);

        if ((color_set_a != nullptr) || (color_set_b != nullptr)){

            size_t prev_color_id = 0xffffffffffffffff;
            size_t prev_km_dist = 0xffffffffffffffff;

            const size_t k = um_a.getGraph()->getK();
            const size_t um_a_km_sz = um_a.size - k + 1;
            const size_t um_b_km_sz = um_b.size - k + 1;

            if (color_set_a != nullptr){

                UnitigColors csd_rev;

                UnitigColorMap<U> new_um_a(0, 0, um_a.size + um_b.size - k + 1, um_a.strand);

                if (!um_a.strand){

                    csd_rev = color_set_a->reverse(um_a);
                    color_set_a = &csd_rev;
                }

                UnitigColors::const_iterator it(color_set_a->begin(0, um_a_km_sz, um_a_km_sz));
                const UnitigColors::const_iterator it_end(color_set_a->end());

                if (it != it_end){

                    prev_km_dist = it.getKmerPosition();
                    prev_color_id = it.getColorID();

                    new_um_a.dist = prev_km_dist;
                    new_um_a.len = 1;

                    ++it;
                }

                // Insert colors layer by layer
                for (; it != it_end; ++it){

                    const size_t km_dist = it.getKmerPosition();
                    const size_t color_id = it.getColorID();

                    if ((color_id != prev_color_id) || (km_dist != prev_km_dist + 1)){

                        new_cs.add(new_um_a, prev_color_id);

                        new_um_a.dist = km_dist;
                        new_um_a.len = 1;
                    }
                    else ++(new_um_a.len);

                    prev_color_id = color_id;
                    prev_km_dist = km_dist;
                }

                if (new_um_a.dist + new_um_a.len != 0) new_cs.add(new_um_a, prev_color_id);
            }

            if (color_set_b != nullptr){

                UnitigColors css_rev;

                UnitigColorMap<U> new_um_b(0, 0, um_a.size + um_b.size - k + 1, um_b.strand);

                if (!um_b.strand){

                    css_rev = color_set_b->reverse(um_b);
                    color_set_b = &css_rev;
                }

                UnitigColors::const_iterator it(color_set_b->begin(0, um_b_km_sz, um_b_km_sz));
                const UnitigColors::const_iterator it_end(color_set_b->end());

                if (it != it_end){

                    prev_km_dist = it.getKmerPosition();
                    prev_color_id = it.getColorID();

                    new_um_b.dist = prev_km_dist + um_a.size - k + 1;
                    new_um_b.len = 1;

                    ++it;
                }

                // Insert colors layer by layer
                for (; it != it_end; ++it){

                    const size_t km_dist = it.getKmerPosition();
                    const size_t color_id = it.getColorID();

                    if ((color_id != prev_color_id) || (km_dist != prev_km_dist + 1)){

                        new_cs.add(new_um_b, prev_color_id);

                        new_um_b.dist = km_dist + um_a.size - k + 1;
                        new_um_b.len = 1;
                    }
                    else ++(new_um_b.len);

                    prev_color_id = color_id;
                    prev_km_dist = km_dist;
                }

                if (new_um_b.dist + new_um_b.len != 0) new_cs.add(new_um_b, prev_color_id);
            }
        }
    }

    return new_cs;
}

template<typename U>
bool DataAccessor<U>::mergeUnitigColors(const UnitigColorMap<U>& um_a, const const_UnitigColorMap<U>& um_b) {

    if (!um_a.isEmpty && !um_b.isEmpty && (um_b.len == um_a.len) && (um_a.getGraph() == um_b.getGraph())){

        const size_t k = um_a.getGraph()->getK();

        UnitigColors& cs_a = um_a.getData()->getUnitigColors(um_a);

        const UnitigColors& cs_b = um_b.getData()->getUnitigColors(um_b);

        UnitigColorMap<U> um_tmp(0, 1, um_a.size, um_a.strand);

        if ((um_b.len * nb_colors_b * 16) < cs_b.size(um_b)){

            const size_t pos_end_b = um_b.dist + um_b.len;
            const size_t sz_km_b = um_b.size - k + 1;

            for (size_t colorID = 0; colorID < nb_colors_b; ++colorID){

                for (size_t km_dist = um_b.dist; km_dist < pos_end_b; ++km_dist){

                    if (cs_b->contains(colorID * sz_km_b + km_dist)){

                        if (um_a.strand != um_b.strand) um_tmp.dist = (um_b.len - 1 - (km_dist - um_b.dist)) + um_a.dist;
                        else um_tmp.dist = (km_dist - um_b.dist) + um_a.dist;

                        cs_a.add(um_tmp, colorID);
                    }
                }
            }
        }
        else {

            UnitigColors::const_iterator it(cs_b.begin(um_b));
            const UnitigColors::const_iterator it_end(cs_b.end());

            while (it != it_end){

                const size_t km_dist_b = it.getKmerPosition();

                if (um_a.strand != um_b.strand) um_tmp.dist = (um_b.len - 1 - (km_dist_b - um_b.dist)) + um_a.dist;
                else um_tmp.dist = (km_dist_b - um_b.dist) + um_a.dist;

                cs_a.add(um_tmp, it.getColorID());

                ++it;
            }
        }

        return true;
    }

    return false;
}

#endif
