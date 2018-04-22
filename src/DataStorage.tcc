#ifndef BFG_DATASTORAGE_TCC
#define BFG_DATASTORAGE_TCC

template<typename U>
DataStorage<U>::DataStorage() : color_sets(nullptr), unitig_cs_link(nullptr), data(nullptr), nb_seeds(0), nb_color_sets(0), nb_elem(0), nb_free_elem(0) {

    std::random_device rd; //Seed
    std::default_random_engine generator(rd()); //Random number generator
    std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF); //Distribution on which to apply the generator

    for (size_t i = 0; i != 256; ++i) seeds[i] = distribution(generator); //Initialize the hash function seeds
}

template<typename U>
DataStorage<U>::DataStorage(const size_t nb_seeds_, const size_t nb_elem_, const vector<string>& color_names_) :
                            color_sets(nullptr), unitig_cs_link(nullptr), data(nullptr), nb_seeds(nb_seeds_), nb_elem(nb_elem_),
                            nb_free_elem(0), nb_color_sets(nb_elem_), color_names(color_names_) {

    std::random_device rd; //Seed
    std::default_random_engine generator(rd()); //Random number generator
    std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF); //Distribution on which to apply the generator

    const size_t sz_unitig_cs_link = (nb_elem >> 6) + ((nb_elem & 0x3F) != 0);

    for (size_t i = 0; i != 256; ++i) seeds[i] = distribution(generator); //Initialize the hash function seeds

    color_sets = new UnitigColors[nb_elem];
    unitig_cs_link = new atomic<uint64_t>[sz_unitig_cs_link];
    data = new U[nb_elem];

    for (size_t i = 0; i != sz_unitig_cs_link; ++i) unitig_cs_link[i] = 0;
}

template<>
inline DataStorage<void>::DataStorage(const size_t nb_seeds_, const size_t nb_elem_, const vector<string>& color_names_) :
                                color_sets(nullptr), unitig_cs_link(nullptr), data(nullptr), nb_seeds(nb_seeds_), nb_elem(nb_elem_),
                                nb_free_elem(0), nb_color_sets(nb_elem_), color_names(color_names_) {

    std::random_device rd; //Seed
    std::default_random_engine generator(rd()); //Random number generator
    std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF); //Distribution on which to apply the generator

    const size_t sz_unitig_cs_link = (nb_elem >> 6) + ((nb_elem & 0x3F) != 0);

    for (int i = 0; i < 256; ++i) seeds[i] = distribution(generator); //Initialize the hash function seeds

    color_sets = new UnitigColors[nb_elem];
    unitig_cs_link = new atomic<uint64_t>[sz_unitig_cs_link];

    for (size_t i = 0; i != sz_unitig_cs_link; ++i) unitig_cs_link[i] = 0;
}

template<typename U>
DataStorage<U>::DataStorage(const DataStorage& o) : nb_seeds(o.nb_seeds), nb_color_sets(o.nb_color_sets), nb_elem(o.nb_elem),
                                                    nb_free_elem(o.nb_free_elem), color_names(o.color_names), data(nullptr),
                                                    color_sets(nullptr), unitig_cs_link(nullptr), overflow(o.overflow) {

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    if ((o.color_sets != nullptr) && (o.nb_elem != 0)){

        color_sets = new UnitigColors[nb_elem];
        copy(o.color_sets, o.color_sets + nb_elem, color_sets);

        const size_t sz_link = (nb_elem >> 6) + ((nb_elem & 0x3F) != 0);

        unitig_cs_link = new atomic<uint64_t>[sz_link];

        for (size_t i = 0; i != sz_link; ++i) unitig_cs_link[i] = o.sz_link[i].load();
    }

    if ((o.data != nullptr) && (o.nb_elem != 0)){

        data = new U[nb_elem];

        copy(o.data, o.data + nb_elem, data);
    }
}

template<>
inline DataStorage<void>::DataStorage(const DataStorage& o) :  nb_seeds(o.nb_seeds), nb_color_sets(o.nb_color_sets), nb_elem(o.nb_elem),
                                                        nb_free_elem(o.nb_free_elem), color_names(o.color_names), data(nullptr),
                                                        color_sets(nullptr), unitig_cs_link(nullptr), overflow(o.overflow) {

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    if ((o.color_sets != nullptr) && (o.nb_elem != 0)){

        color_sets = new UnitigColors[nb_elem];
        copy(o.color_sets, o.color_sets + nb_elem, color_sets);

        const size_t sz_link = (nb_elem >> 6) + ((nb_elem & 0x3F) != 0);

        unitig_cs_link = new atomic<uint64_t>[sz_link];

        for (size_t i = 0; i != sz_link; ++i) unitig_cs_link[i] = o.unitig_cs_link[i].load();
    }
}

template<typename U>
DataStorage<U>::DataStorage(DataStorage&& o) :  nb_seeds(o.nb_seeds), nb_color_sets(o.nb_color_sets), nb_elem(o.nb_elem),
                                                nb_free_elem(o.nb_free_elem), color_sets(o.color_sets), unitig_cs_link(o.unitig_cs_link),
                                                data(o.data), color_names(move(o.color_names)), overflow(move(o.overflow)) {

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    o.color_sets = nullptr;
    o.unitig_cs_link = nullptr;
    o.data = nullptr;

    o.clear();
}

template<typename U>
DataStorage<U>::~DataStorage() {

    empty();
}

template<typename U>
void DataStorage<U>::clear() {

    nb_seeds = 0;
    nb_color_sets = 0;
    nb_elem = 0;
    nb_free_elem = 0;

    empty();
}

template<typename U>
void DataStorage<U>::empty() {

    if (color_sets != nullptr){

        delete[] color_sets;
        color_sets = nullptr;
    }

    if (unitig_cs_link != nullptr){

        delete[] unitig_cs_link;
        unitig_cs_link = nullptr;
    }

    if (data != nullptr){

        delete[] data;
        data = nullptr;
    }

    color_names.clear();
    overflow.clear();
}

template<>
inline void DataStorage<void>::empty() {

    if (color_sets != nullptr){

        delete[] color_sets;
        color_sets = nullptr;
    }

    if (unitig_cs_link != nullptr){

        delete[] unitig_cs_link;
        unitig_cs_link = nullptr;
    }

    data = nullptr;

    color_names.clear();
    overflow.clear();
}

template<typename U>
DataStorage<U>& DataStorage<U>::operator=(const DataStorage& o) {

    empty();

    nb_seeds = o.nb_seeds;
    nb_color_sets = o.nb_color_sets;
    nb_elem = o.nb_elem;
    nb_free_elem = o.nb_free_elem;

    color_names = o.color_names;

    overflow = o.overflow;

    color_sets = nullptr;
    unitig_cs_link = nullptr;
    data = nullptr;

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    if ((o.color_sets != nullptr) && (o.nb_elem != 0)){

        color_sets = new UnitigColors[nb_elem];
        copy(o.color_sets, o.color_sets + nb_elem, color_sets);

        const size_t sz_link = (nb_elem >> 6) + ((nb_elem & 0x3F) != 0);

        unitig_cs_link = new atomic<uint64_t>[sz_link];
        for (size_t i = 0; i != sz_link; ++i) unitig_cs_link[i] = o.unitig_cs_link[i].load();
    }

    if ((o.data != nullptr) && (o.nb_elem != 0)){

        data = new U[nb_elem];

        copy(o.data, o.data + nb_elem, data);
    }

    return *this;
}

template<>
inline DataStorage<void>& DataStorage<void>::operator=(const DataStorage& o) {

    empty();

    nb_seeds = o.nb_seeds;
    nb_color_sets = o.nb_color_sets;
    nb_elem = o.nb_elem;
    nb_free_elem = o.nb_free_elem;

    color_names = o.color_names;

    overflow = o.overflow;

    color_sets = nullptr;
    data = nullptr;

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    if ((o.color_sets != nullptr) && (o.nb_elem != 0)){

        color_sets = new UnitigColors[nb_elem];
        copy(o.color_sets, o.color_sets + nb_elem, color_sets);

        const size_t sz_link = (nb_elem >> 6) + ((nb_elem & 0x3F) != 0);

        unitig_cs_link = new atomic<uint64_t>[sz_link];

        for (size_t i = 0; i != sz_link; ++i) unitig_cs_link[i] = o.unitig_cs_link[i].load();
    }

    return *this;
}

template<typename U>
DataStorage<U>& DataStorage<U>::operator=(DataStorage&& o) {

    if (this != &o) {

        empty();

        nb_seeds = o.nb_seeds;
        nb_color_sets = o.nb_color_sets;
        nb_elem = o.nb_elem;
        nb_free_elem = o.nb_free_elem;

        color_names = move(o.color_names);

        overflow = move(o.overflow);

        color_sets = o.color_sets;
        unitig_cs_link = o.unitig_cs_link;
        data = o.data;

        memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

        o.color_sets = nullptr;
        o.unitig_cs_link = nullptr;
        o.data = nullptr;

        o.clear();
    }

    return *this;
}

template<>
inline DataStorage<void>& DataStorage<void>::operator=(DataStorage&& o) {

    if (this != &o) {

        empty();

        nb_seeds = o.nb_seeds;
        nb_color_sets = o.nb_color_sets;
        nb_elem = o.nb_elem;
        nb_free_elem = o.nb_free_elem;

        color_names = move(o.color_names);

        overflow = move(o.overflow);

        color_sets = o.color_sets;
        unitig_cs_link = o.unitig_cs_link;

        memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

        o.color_sets = nullptr;
        o.unitig_cs_link = nullptr;

        o.clear();
    }

    return *this;
}

template<typename U>
const UnitigColors* DataStorage<U>::getUnitigColors(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (color_sets != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            typename KmerHashTable<size_t>::const_iterator it = overflow.find(head);

            if (it != overflow.end()) return &color_sets[*it];
        }
        else return &(color_sets[head.hash(seeds[da_id - 1]) % nb_color_sets]);
    }

    return nullptr;
}

template<typename U>
UnitigColors* DataStorage<U>::getUnitigColors(const UnitigColorMap<U>& um) {

    if (!um.isEmpty && (color_sets != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            typename KmerHashTable<size_t>::iterator it = overflow.find(head);

            if (it != overflow.end()) return &color_sets[*it];
        }
        else return &(color_sets[head.hash(seeds[da_id - 1]) % nb_color_sets]);
    }

    return nullptr;
}

template<typename U>
const U* DataStorage<U>::getData(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (data != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            typename KmerHashTable<size_t>::const_iterator it = overflow.find(head);

            if (it != overflow.end()) return &data[*it];
        }
        else return &(data[head.hash(seeds[da_id - 1]) % nb_color_sets]);
    }

    return nullptr;
}

template<> inline const void* DataStorage<void>::getData(const const_UnitigColorMap<void>& um) const {

    return nullptr;
}

template<typename U>
U* DataStorage<U>::getData(const UnitigColorMap<U>& um) {

    if (!um.isEmpty && (data != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            typename KmerHashTable<size_t>::iterator it = overflow.find(head);

            if (it != overflow.end()) return &data[*it];
        }
        else return &(data[head.hash(seeds[da_id - 1]) % nb_color_sets]);
    }

    return nullptr;
}

template<> inline void* DataStorage<void>::getData(const UnitigColorMap<void>& um) {

    return nullptr;
}

template<typename U>
bool DataStorage<U>::joinUnitigColors(const UnitigColorMap<U>& um_dest, const UnitigColorMap<U>& um_src) {

    if (!um_dest.isEmpty && !um_src.isEmpty && (color_sets != nullptr)){

        UnitigColors* color_set_dest = getUnitigColors(um_dest);
        UnitigColors* color_set_src = getUnitigColors(um_src);

        if ((color_set_dest != nullptr) && (color_set_src != nullptr)){

            UnitigColors new_cs, csd_rev, css_rev;

            size_t prev_color_id = 0xffffffffffffffff;
            size_t prev_km_dist = 0xffffffffffffffff;

            const size_t um_dest_km_sz = um_dest.size - um_dest.getCompactedDBG()->getK() + 1;
            const size_t um_src_km_sz = um_src.size - um_src.getCompactedDBG()->getK() + 1;

            UnitigColorMap<U> new_um_dest(0, 0, um_dest.size + um_src.size, um_dest.strand);

            if (!um_dest.strand){

                csd_rev = color_set_dest->reverse(um_dest);
                color_set_dest = &csd_rev;
            }

            UnitigColors::const_iterator it = color_set_dest->begin(), it_end = color_set_dest->end();

            if (it != it_end){

                prev_color_id = it->getColorID(um_dest_km_sz);
                prev_km_dist = it->getKmerPosition(um_dest_km_sz);

                new_um_dest.dist = prev_km_dist;
                new_um_dest.len = 1;

                ++it;
            }

            // Insert colors layer by layer
            for (; it != it_end; ++it){

                const size_t color_id = it->getColorID(um_dest_km_sz);
                const size_t km_dist = it->getKmerPosition(um_dest_km_sz);

                if ((color_id != prev_color_id) || (km_dist != prev_km_dist + 1)){

                    new_cs.add(new_um_dest, color_id);

                    new_um_dest.dist = km_dist;
                    new_um_dest.len = 1;
                }
                else ++new_um_dest.len;

                prev_color_id = color_id;
                prev_km_dist = km_dist;
            }

            if (new_um_dest.dist + new_um_dest.len != 0) new_cs.add(new_um_dest, prev_color_id);

            UnitigColorMap<U> new_um_src(0, 0, um_dest.size + um_src.size, um_src.strand);

            if (!um_src.strand){

                css_rev = color_set_src->reverse(um_src);
                color_set_src = &css_rev;
            }

            it = color_set_src->begin(), it_end = color_set_src->end();

            if (it != it_end){

                prev_color_id = it->getColorID(um_src_km_sz);
                prev_km_dist = it->getKmerPosition(um_src_km_sz);

                new_um_src.dist = prev_km_dist + um_dest.size;
                new_um_src.len = 1;

                ++it;
            }

            // Insert colors layer by layer
            for (; it != it_end; ++it){

                const size_t color_id = it->getColorID(um_src_km_sz);
                const size_t km_dist = it->getKmerPosition(um_src_km_sz);

                if ((color_id != prev_color_id) || (km_dist != prev_km_dist + 1)){

                    new_cs.add(new_um_src, color_id);

                    new_um_src.dist = km_dist + um_dest.size;
                    new_um_src.len = 1;
                }
                else ++new_um_src.len;

                prev_color_id = color_id;
                prev_km_dist = km_dist;
            }

            if (new_um_src.dist + new_um_src.len != 0) new_cs.add(new_um_src, prev_color_id);

            *color_set_dest = new_cs;

            return true;
        }
    }

    return false;
}

template<typename U>
UnitigColors DataStorage<U>::getSubUnitigColors(const UnitigColorMap<U>& um) const {

    UnitigColors new_cs;

    if (!um.isEmpty && (color_sets != nullptr)){

        const UnitigColors* cs = getUnitigColors(um);

        if (cs != nullptr){

            const size_t end = um.dist + um.len;
            const size_t um_km_sz = um.size - um.getCompactedDBG()->getK() + 1;

            UnitigColorMap<U> um_tmp(0, 1, um.len, um.strand);

            for (UnitigColors::const_iterator it = cs->begin(), it_end = cs->end(); it != it_end; ++it){

                const size_t color_id = it->getColorID(um_km_sz);
                const size_t km_dist = it->getKmerPosition(um_km_sz);

                if ((km_dist >= um.dist) && (km_dist < end)){

                    um_tmp.dist = km_dist - um.dist;

                    new_cs.add(um_tmp, color_id);
                }
            }
        }
    }

    return new_cs;
}

template<typename U>
vector<string> DataStorage<U>::getSubUnitigColorNames(const UnitigColorMap<U>& um) const {

    vector<string> v_out;

    if (!um.isEmpty && (color_sets != nullptr)){

        const UnitigColors* cs = getUnitigColors(um);

        if (cs != nullptr){

            size_t prev_color_id = 0xffffffffffffffff;

            const size_t end = um.dist + um.len;
            const size_t um_km_sz = um.size - um.getCompactedDBG()->getK() + 1;

            for (UnitigColors::const_iterator it = cs->begin(), it_end = cs->end(); it != it_end; ++it){

                const size_t color_id = it->getColorID(um_km_sz);
                const size_t km_dist = it->getKmerPosition(um_km_sz);

                if ((km_dist >= um.dist) && (km_dist < end) && (color_id != prev_color_id)){

                    v_out.push_back(color_names[color_id]);
                }

                prev_color_id = color_id;
            }
        }
    }

    return v_out;
}

template<typename U>
UnitigColors* DataStorage<U>::insert() {

    if (nb_free_elem == 0){

        UnitigColors* old_color_sets = color_sets;
        atomic<uint64_t>* old_unitig_cs_link = unitig_cs_link;
        U* old_data = data;

        const size_t old_nb_elem = nb_elem;
        const size_t old_sz_link = (old_nb_elem >> 6) + ((old_nb_elem & 0x3F) != 0);

        nb_free_elem = nb_elem * 0.1;
        nb_elem += nb_free_elem;

        const size_t sz_link = (nb_elem >> 6) + ((nb_elem & 0x3F) != 0);

        color_sets = new UnitigColors[nb_elem];

        move(color_sets, color_sets + old_nb_elem, old_color_sets);
        delete[] old_color_sets;

        unitig_cs_link = new atomic<uint64_t>[sz_link];

        for (size_t i = 0; i != old_sz_link; ++i) unitig_cs_link[i] = old_unitig_cs_link[i].load();
        for (size_t i = old_sz_link; i != sz_link; ++i) unitig_cs_link[i] = 0;

        delete[] old_unitig_cs_link;

        data = new U[nb_elem];

        move(data, data + old_nb_elem, old_data);
        delete[] old_data;
    }

    return &color_sets[nb_elem - nb_free_elem--];
}

template<>
inline UnitigColors* DataStorage<void>::insert() {

    if (nb_free_elem == 0){

        UnitigColors* old_color_sets = color_sets;
        atomic<uint64_t>* old_unitig_cs_link = unitig_cs_link;

        const size_t old_nb_elem = nb_elem;
        const size_t old_sz_link = (old_nb_elem >> 6) + ((old_nb_elem & 0x3F) != 0);

        nb_free_elem = nb_elem * 0.1;
        nb_elem += nb_free_elem;

        const size_t sz_link = (nb_elem >> 6) + ((nb_elem & 0x3F) != 0);

        color_sets = new UnitigColors[nb_elem];

        move(color_sets, color_sets + old_nb_elem, old_color_sets);
        delete[] old_color_sets;

        unitig_cs_link = new atomic<uint64_t>[sz_link];

        for (size_t i = 0; i != old_sz_link; ++i) unitig_cs_link[i] = old_unitig_cs_link[i].load();
        for (size_t i = old_sz_link; i != sz_link; ++i) unitig_cs_link[i] = 0;

        delete[] old_unitig_cs_link;
    }

    return &color_sets[nb_elem - nb_free_elem--];
}

template<typename U>
bool DataStorage<U>::write(const string prefix_output_filename, const size_t nb_threads, const bool verbose) const {

    if (verbose) cout << endl << "DataStorage::write(): Writing colors to disk" << endl;

    const string out = prefix_output_filename + ".bfg_colors";

    FILE* fp = fopen(out.c_str(), "wb");

    if (fp == NULL) {

        cerr << "DataStorage::write(): Could not open file " << out << " for writing color sets" << endl;
        return false;
    }
    else {

        fclose(fp);

        if (std::remove(out.c_str()) != 0) cerr << "DataStorage::write(): Could not remove temporary file " << out << endl;
    }

    ofstream colorsfile_out;
    ostream colors_out(nullptr);

    colorsfile_out.open(out.c_str(), ios_base::out | ios_base::binary);
    colors_out.rdbuf(colorsfile_out.rdbuf());

    const size_t format_version = BFG_COLOREDCDBG_FORMAT_VERSION;
    const size_t overflow_sz = overflow.size();
    const size_t nb_colors = color_names.size();

    //Write the file format version number
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&format_version), sizeof(size_t));
    //Write number of different seeds for hash function
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_seeds), sizeof(size_t));
   //Write number of colors in the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_colors), sizeof(size_t));
    //Write number of color sets in the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_color_sets), sizeof(size_t));
    //Write number of elements allocated
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_elem), sizeof(size_t));
    //Write number of free elements allocated
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_free_elem), sizeof(size_t));
    //Write number of (kmer, color set) overflowing
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&overflow_sz), sizeof(size_t));

    for (size_t i = 0; (i < nb_seeds) && colors_out.good(); ++i){
        //Write the hash function seeds of the graph
        colors_out.write(reinterpret_cast<const char*>(&seeds[i]), sizeof(uint64_t));
    }

    for (size_t i = 0; (i < nb_colors) && colors_out.good(); ++i){
        //Write the color names of the graph
        colors_out.write(color_names[i].c_str(), color_names[i].length() + 1);
    }

    for (uint64_t i = 0, j = ((nb_elem >> 6) + ((nb_elem & 0x3F) != 0)), e; (i != j) && colors_out.good(); ++i){

        e = unitig_cs_link[i].load();
        colors_out.write(reinterpret_cast<const char*>(&e), sizeof(uint64_t));
    }

    for (size_t i = 0; (i < nb_elem) && colors_out.good(); ++i) color_sets[i].write(colors_out); //Write the color sets

    typename KmerHashTable<size_t>::const_iterator it = overflow.begin();
    typename KmerHashTable<size_t>::const_iterator it_end = overflow.end();

    for (; (it != it_end) && colors_out.good(); ++it){

        it.getKey().write(colors_out);
        colors_out.write(reinterpret_cast<const char*>(&(*it)), sizeof(size_t));
    }

    const bool ret = colors_out.good();

    colorsfile_out.close();

    return ret;
}

template<typename U>
bool DataStorage<U>::read(const string& filename_colors, bool verbose) {

    if (verbose) cout << endl << "DataStorage::read(): Reading color sets from disk" << endl;

    FILE* fp = fopen(filename_colors.c_str(), "rb");

    if (fp == NULL) {

        cerr << "DataStorage::read(): Could not open file " << filename_colors << " for reading color sets" << endl;
        return false;
    }
    else fclose(fp);

    Kmer km;

    size_t format_version, overflow_sz, nb_colors;

    ifstream colorsfile_in;
    istream colors_in(nullptr);

    char* buffer = nullptr;

    clear();

    colorsfile_in.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
    colors_in.rdbuf(colorsfile_in.rdbuf());

    //Read the file format version number
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&format_version), sizeof(size_t));
    //Read number of different seeds for hash function
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_seeds), sizeof(size_t));
    //Read number of colors in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_colors), sizeof(size_t));
    //Read number of color sets in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_color_sets), sizeof(size_t));
    //Read number of elements allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_elem), sizeof(size_t));
    //Read number of free elements allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_free_elem), sizeof(size_t));
    //Read number of (kmer, color set) overflowing
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&overflow_sz), sizeof(size_t));

    if (nb_seeds >= 256){

        cerr << "DataStorage::read(): Does not support more than 255 hash seeds" << endl;
        return false;
    }

    const size_t sz_unitig_cs_link = (nb_elem >> 6) + ((nb_elem & 0x3F) != 0);

    overflow = KmerHashTable<size_t>(overflow_sz);

    color_sets = new UnitigColors[nb_elem];
    unitig_cs_link = new atomic<uint64_t>[sz_unitig_cs_link];
    data = new U[nb_elem];

    for (size_t i = 0; (i < nb_seeds) && colors_in.good(); ++i){
        //Read the hash function seeds of the graph
        colors_in.read(reinterpret_cast<char*>(&seeds[i]), sizeof(uint64_t));
    }

    buffer = new char[1000];

    for (size_t i = 0; (i < nb_colors) && colors_in.good(); ++i){
        //Read the hash function seeds of the graph
        colors_in.getline(buffer, 1000);
        color_names.push_back(string(buffer));
    }

    delete[] buffer;

    for (uint64_t i = 0, e; (i != sz_unitig_cs_link) && colors_in.good(); ++i){

        colors_in.read(reinterpret_cast<char*>(&e), sizeof(uint64_t));
        unitig_cs_link[i] = e;
    }

    for (size_t i = 0; (i < nb_elem) && colors_in.good(); ++i) color_sets[i].read(colors_in);

    for (size_t i = 0, pos; (i < overflow_sz) && colors_in.good(); ++i){

        km.read(colors_in);
        colors_in.read(reinterpret_cast<char*>(&pos), sizeof(size_t));

        std::pair<typename KmerHashTable<size_t>::iterator, bool> p = overflow.insert(km, pos);
    }

    const bool ret = colors_in.good();

    colorsfile_in.close();

    return true;
}

template<>
inline bool DataStorage<void>::read(const string& filename_colors, bool verbose) {

    if (verbose) cout << endl << "DataStorage::read(): Reading color sets from disk" << endl;

    FILE* fp = fopen(filename_colors.c_str(), "rb");

    if (fp == NULL) {

        cerr << "DataStorage::read(): Could not open file " << filename_colors << " for reading color sets" << endl;
        return false;
    }
    else fclose(fp);

    Kmer km;

    size_t format_version, overflow_sz, nb_colors;

    ifstream colorsfile_in;
    istream colors_in(nullptr);

    char* buffer = nullptr;

    clear();

    colorsfile_in.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
    colors_in.rdbuf(colorsfile_in.rdbuf());

    //Read the file format version number
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&format_version), sizeof(size_t));
    //Read number of different seeds for hash function
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_seeds), sizeof(size_t));
    //Read number of colors in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_colors), sizeof(size_t));
    //Read number of color sets in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_color_sets), sizeof(size_t));
    //Read number of elements allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_elem), sizeof(size_t));
    //Read number of free elements allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_free_elem), sizeof(size_t));
    //Read number of (kmer, color set) overflowing
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&overflow_sz), sizeof(size_t));

    if (nb_seeds >= 256){

        cerr << "DataStorage::read(): Does not support more than 255 hash seeds" << endl;
        return false;
    }

    const size_t sz_unitig_cs_link = (nb_elem >> 6) + ((nb_elem & 0x3F) != 0);

    overflow = KmerHashTable<size_t>(overflow_sz);

    color_sets = new UnitigColors[nb_elem];
    unitig_cs_link = new atomic<uint64_t>[sz_unitig_cs_link];

    for (size_t i = 0; (i < nb_seeds) && colors_in.good(); ++i){
        //Read the hash function seeds of the graph
        colors_in.read(reinterpret_cast<char*>(&seeds[i]), sizeof(uint64_t));
    }

    buffer = new char[1000];

    for (size_t i = 0; (i < nb_colors) && colors_in.good(); ++i){
        //Read the hash function seeds of the graph
        colors_in.getline(buffer, 1000);
        color_names.push_back(string(buffer));
    }

    delete[] buffer;

    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(unitig_cs_link), sz_unitig_cs_link * sizeof(atomic<uint64_t>));

    for (size_t i = 0; (i < nb_elem) && colors_in.good(); ++i) color_sets[i].read(colors_in);

    for (size_t i = 0, pos; (i < overflow_sz) && colors_in.good(); ++i){

        km.read(colors_in);
        colors_in.read(reinterpret_cast<char*>(&pos), sizeof(size_t));

        std::pair<typename KmerHashTable<size_t>::iterator, bool> p = overflow.insert(km, pos);
    }

    const bool ret = colors_in.good();

    colorsfile_in.close();

    return true;
}

template<typename U>
uint64_t DataStorage<U>::getHash(const UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (color_sets != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            typename KmerHashTable<size_t>::const_iterator it = overflow.find(head);

            if (it != overflow.end()) return *it;
        }
        else return head.hash(seeds[da_id - 1]) % nb_color_sets;
    }

    return 0;
}

#endif
