#ifndef BIFROST_DATA_STORAGE_TCC
#define BIFROST_DATA_STORAGE_TCC

DataStorage::DataStorage(const vector<string>& color_names_) : color_names(color_names_) {}
DataStorage::DataStorage(const DataStorage& o) : color_names(o.color_names) {}
DataStorage::DataStorage(DataStorage&& o) : color_names(move(o.color_names)) {}

void DataStorage::clear() {

    color_names.clear();
}

DataStorage& DataStorage::operator=(const DataStorage& o) {

    color_names = o.color_names;
}

DataStorage& DataStorage::operator=(DataStorage&& o) {

    if (this != &o) color_names = move(o.color_names);

    return *this;
}


bool DataStorage::write(const string& prefix_output_filename, const bool verbose) const {

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
    //colors_out.sync_with_stdio(false);

    const size_t format_version = BFG_COLOREDCDBG_FORMAT_VERSION;
    const size_t overflow_sz = overflow.size();
    const size_t nb_colors = color_names.size();

    const size_t block_sz = 1024;

    const char nl = '\n';

    streampos pos_f_cs;

    vector<streampos> v_pos_f_cs;

    //Write the file format version number
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&format_version), sizeof(size_t));
    //Write number of different seeds for hash function
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_seeds), sizeof(size_t));
   //Write number of colors in the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_colors), sizeof(size_t));
    //Write number of color sets in the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_cs), sizeof(size_t));
    //Write number of elements allocated
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&sz_cs), sizeof(size_t));
    //Write number of SharedUnitigColors allocated
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&sz_shared_cs), sizeof(size_t));
    //Write number of (kmer, color set) overflowing
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&overflow_sz), sizeof(size_t));
    //Write the hash function seeds of the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(seeds), nb_seeds * sizeof(uint64_t));
    // Write the block size of the color sets
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&block_sz), sizeof(size_t));

    pos_f_cs = colors_out.tellp();

    const size_t nb_pos_shared_cs = (sz_shared_cs / block_sz) + static_cast<size_t>((sz_shared_cs % block_sz) != 0);
    const size_t nb_pos_cs = (sz_cs / block_sz) + static_cast<size_t>((sz_cs % block_sz) != 0);

    for (size_t i = 0; (i < nb_pos_shared_cs) && colors_out.good(); ++i){
        // Reserve space to write positions in file of shared colorsets blocks
        colors_out.write(reinterpret_cast<const char*>(&pos_f_cs), sizeof(streampos));
    }

    for (size_t i = 0; (i < nb_pos_cs) && colors_out.good(); ++i){
        // Reserve space to write positions in file of non-shared colorsets blocks
        colors_out.write(reinterpret_cast<const char*>(&pos_f_cs), sizeof(streampos));
    }

    for (size_t i = 0; (i < nb_colors) && colors_out.good(); ++i){
        //Write the color names of the graph
        colors_out.write(color_names[i].c_str(), color_names[i].size() * sizeof(char));
        colors_out.write(&nl, sizeof(char));
    }

    for (uint64_t i = 0, j = ((sz_cs >> 6) + ((sz_cs & 0x3F) != 0)), e; (i != j) && colors_out.good(); ++i){

        e = unitig_cs_link[i].load();
        colors_out.write(reinterpret_cast<const char*>(&e), sizeof(uint64_t));
    }

    for (size_t i = 0; (i < sz_shared_cs) && colors_out.good(); ++i){

        if (i % block_sz == 0) v_pos_f_cs.push_back(colors_out.tellp());

        if (shared_color_sets[i].first.write(colors_out)){

            colors_out.write(reinterpret_cast<const char*>(&(shared_color_sets[i].second)), sizeof(size_t));
        }
    }

    for (size_t i = 0; (i < sz_cs) && colors_out.good(); ++i){

        if (i % block_sz == 0) v_pos_f_cs.push_back(colors_out.tellp());

        color_sets[i].write(colors_out, false); //Write the color sets
    }

    unordered_map<pair<Kmer, size_t>, size_t>::const_iterator it(overflow.begin());
    const unordered_map<pair<Kmer, size_t>, size_t>::const_iterator it_end(overflow.end());

    for (; (it != it_end) && colors_out.good(); ++it){

        it->first.first.write(colors_out); // Write the k-mer

        colors_out.write(reinterpret_cast<const char*>(&(it->first.second)), sizeof(size_t)); // Write the unitig length
        colors_out.write(reinterpret_cast<const char*>(&(it->second)), sizeof(size_t)); // Write the position
    }

    if (colors_out.good()){

        colors_out.seekp(pos_f_cs); // Re-position cursor to array of position of shared color sets at the beginning

        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&v_pos_f_cs[0]), v_pos_f_cs.size() * sizeof(streampos));
    }

    const bool ret = colors_out.good();

    colorsfile_out.close();

    return ret;
}

template<>
inline bool DataStorage<void>::write(const string& prefix_output_filename, const bool verbose) const {

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
    //colors_out.sync_with_stdio(false);

    const size_t format_version = BFG_COLOREDCDBG_FORMAT_VERSION;
    const size_t overflow_sz = overflow.size();
    const size_t nb_colors = color_names.size();

    const size_t block_sz = 1024;

    const char nl = '\n';

    streampos pos_f_cs;

    vector<streampos> v_pos_f_cs;

    //Write the file format version number
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&format_version), sizeof(size_t));
    //Write number of different seeds for hash function
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_seeds), sizeof(size_t));
   //Write number of colors in the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_colors), sizeof(size_t));
    //Write number of color sets in the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_cs), sizeof(size_t));
    //Write number of elements allocated
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&sz_cs), sizeof(size_t));
    //Write number of SharedUnitigColors allocated
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&sz_shared_cs), sizeof(size_t));
    //Write number of (kmer, color set) overflowing
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&overflow_sz), sizeof(size_t));
    //Write the hash function seeds of the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(seeds), nb_seeds * sizeof(uint64_t));
    // Write the block size of the color sets
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&block_sz), sizeof(size_t));

    pos_f_cs = colors_out.tellp();

    const size_t nb_pos_shared_cs = (sz_shared_cs / block_sz) + static_cast<size_t>((sz_shared_cs % block_sz) != 0);
    const size_t nb_pos_cs = (sz_cs / block_sz) + static_cast<size_t>((sz_cs % block_sz) != 0);

    for (size_t i = 0; (i < nb_pos_shared_cs) && colors_out.good(); ++i){
        // Reserve space to write positions in file of shared colorsets blocks
        colors_out.write(reinterpret_cast<const char*>(&pos_f_cs), sizeof(streampos));
    }

    for (size_t i = 0; (i < nb_pos_cs) && colors_out.good(); ++i){
        // Reserve space to write positions in file of non-shared colorsets blocks
        colors_out.write(reinterpret_cast<const char*>(&pos_f_cs), sizeof(streampos));
    }

    for (size_t i = 0; (i < nb_colors) && colors_out.good(); ++i){
        //Write the color names of the graph
        colors_out.write(color_names[i].c_str(), color_names[i].size() * sizeof(char));
        colors_out.write(&nl, sizeof(char));
    }

    for (uint64_t i = 0, j = ((sz_cs >> 6) + ((sz_cs & 0x3F) != 0)), e; (i != j) && colors_out.good(); ++i){

        e = unitig_cs_link[i].load();
        colors_out.write(reinterpret_cast<const char*>(&e), sizeof(uint64_t));
    }

    for (size_t i = 0; (i < sz_shared_cs) && colors_out.good(); ++i){

        if (i % block_sz == 0) v_pos_f_cs.push_back(colors_out.tellp());

        if (shared_color_sets[i].first.write(colors_out)){

            colors_out.write(reinterpret_cast<const char*>(&(shared_color_sets[i].second)), sizeof(size_t));
        }
    }

    for (size_t i = 0; (i < sz_cs) && colors_out.good(); ++i){

        if (i % block_sz == 0) v_pos_f_cs.push_back(colors_out.tellp());

        color_sets[i].write(colors_out, false); //Write the color sets
    }

    unordered_map<pair<Kmer, size_t>, size_t>::const_iterator it(overflow.begin());
    const unordered_map<pair<Kmer, size_t>, size_t>::const_iterator it_end(overflow.end());

    for (; (it != it_end) && colors_out.good(); ++it){

        it->first.first.write(colors_out); // Write the k-mer

        colors_out.write(reinterpret_cast<const char*>(&(it->first.second)), sizeof(size_t)); // Write the unitig length
        colors_out.write(reinterpret_cast<const char*>(&(it->second)), sizeof(size_t)); // Write the position
    }

    if (colors_out.good()){

        colors_out.seekp(pos_f_cs); // Re-position cursor to array of position of shared color sets at the beginning

        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&v_pos_f_cs[0]), v_pos_f_cs.size() * sizeof(streampos));
    }

    const bool ret = colors_out.good();

    colorsfile_out.close();

    return ret;
}


bool DataStorage::read(const string& filename_colors, const size_t nb_threads, const bool verbose) {

    if (verbose) cout << endl << "DataStorage::read(): Reading color sets from disk" << endl;

    FILE* fp = fopen(filename_colors.c_str(), "rb");

    if (fp == NULL) {

        cerr << "DataStorage::read(): Could not open file " << filename_colors << " for reading color sets" << endl;
        return false;
    }
    else fclose(fp);

    auto readSharedColorSets = [](UnitigColors::SharedUnitigColors* shared_color_sets, istream& colors_in, const size_t sz){

        for (size_t i = 0; (i < sz) && colors_in.good(); ++i){

            if (shared_color_sets[i].first.read(colors_in)){

                colors_in.read(reinterpret_cast<char*>(&(shared_color_sets[i].second)), sizeof(size_t));
            }
        }
    };

    auto readColorSets = [](UnitigColors* color_sets, istream& colors_in, const size_t sz){

        for (size_t i = 0; (i < sz) && colors_in.good(); ++i) color_sets[i].read(colors_in);
    };

    Kmer km;

    size_t format_version, overflow_sz, nb_colors;

    ifstream colorsfile_in;
    istream colors_in(nullptr);

    colorsfile_in.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
    colors_in.rdbuf(colorsfile_in.rdbuf());

    clear();

    //Read the file format version number
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&format_version), sizeof(size_t));
    //Read number of different seeds for hash function
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_seeds), sizeof(size_t));
    //Read number of colors in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_colors), sizeof(size_t));
    //Read number of color sets in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_cs), sizeof(size_t));
    //Read number of elements allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&sz_cs), sizeof(size_t));
    //Write number of SharedUnitigColors allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&sz_shared_cs), sizeof(size_t));
    //Read number of (kmer, color set) overflowing
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&overflow_sz), sizeof(size_t));

    if (nb_seeds >= 256){

        cerr << "DataStorage::read(): Does not support more than 255 hash seeds" << endl;
        return false;
    }

    const size_t sz_unitig_cs_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

    overflow = unordered_map<pair<Kmer, size_t>, size_t>(overflow_sz);

    color_sets = new UnitigColors[sz_cs];
    shared_color_sets = new UnitigColors::SharedUnitigColors[sz_shared_cs];
    unitig_cs_link = new atomic<uint64_t>[sz_unitig_cs_link];
    data = new U[sz_cs];

    //Read the hash function seeds of the graph
    colors_in.read(reinterpret_cast<char*>(seeds), nb_seeds * sizeof(uint64_t));

    if (format_version == 1){

        for (size_t i = 0; (i < nb_colors) && colors_in.good(); ++i){
            //Read the hash function seeds of the graph
            color_names.push_back(string());
            getline(colors_in, color_names[i]);
        }

        for (uint64_t i = 0, e; (i != sz_unitig_cs_link) && colors_in.good(); ++i){

            colors_in.read(reinterpret_cast<char*>(&e), sizeof(uint64_t));
            unitig_cs_link[i] = e;
        }

        readSharedColorSets(shared_color_sets, colors_in, sz_shared_cs);
        readColorSets(color_sets, colors_in, sz_cs);
    }
    else {

        size_t block_sz = 0;

        streampos* pos_f_cs = nullptr;

        //Read the hash function seeds of the graph
        if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&block_sz), sizeof(size_t));

        const size_t nb_pos_shared_cs = (sz_shared_cs / block_sz) + static_cast<size_t>((sz_shared_cs % block_sz) != 0);
        const size_t nb_pos_cs = (sz_cs / block_sz) + static_cast<size_t>((sz_cs % block_sz) != 0);
        const size_t pos_f_cs_sz = nb_pos_shared_cs + nb_pos_cs;

        if (pos_f_cs_sz != 0){

            pos_f_cs = new streampos[pos_f_cs_sz];

            if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(pos_f_cs), pos_f_cs_sz * sizeof(streampos));
        }

        for (size_t i = 0; (i < nb_colors) && colors_in.good(); ++i){
            //Read the hash function seeds of the graph
            color_names.push_back(string());
            getline(colors_in, color_names[i]);
        }

        for (uint64_t i = 0, e; (i != sz_unitig_cs_link) && colors_in.good(); ++i){

            colors_in.read(reinterpret_cast<char*>(&e), sizeof(uint64_t));
            unitig_cs_link[i] = e;
        }

        if ((nb_threads == 1) || (pos_f_cs_sz == 0)) {

            readSharedColorSets(shared_color_sets, colors_in, sz_shared_cs);
            readColorSets(color_sets, colors_in, sz_cs);
        }
        else {

            streampos colors_in_pos = colors_in.tellg();

            colorsfile_in.close();

            mutex m_colors_in_pos;

            vector<thread> workers; // need to keep track of threads so we can join them

            std::atomic<size_t> i;

            i = 0;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        ifstream colorsfile_in_t;
                        istream colors_in_t(nullptr);

                        colorsfile_in_t.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
                        colors_in_t.rdbuf(colorsfile_in_t.rdbuf());

                        while (true) {

                            const size_t l_i = i++;

                            if (l_i >= nb_pos_shared_cs){

                                const streampos colors_in_t_pos = colors_in_t.tellg();

                                {
                                    unique_lock<mutex> lock(m_colors_in_pos);

                                    colors_in_pos = max(colors_in_pos, colors_in_t_pos);
                                }

                                colorsfile_in_t.close();

                                break;
                            }

                            colors_in_t.seekg(pos_f_cs[l_i]);
                            readSharedColorSets(shared_color_sets + (l_i * block_sz), colors_in_t, min(block_sz, sz_shared_cs - (l_i * block_sz)));
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();

            workers.clear();

            i = nb_pos_shared_cs;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        ifstream colorsfile_in_t;
                        istream colors_in_t(nullptr);

                        colorsfile_in_t.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
                        colors_in_t.rdbuf(colorsfile_in_t.rdbuf());

                        while (true) {

                            size_t l_i = i++;

                            if (l_i >= pos_f_cs_sz){

                                const streampos colors_in_t_pos = colors_in_t.tellg();

                                {
                                    unique_lock<mutex> lock(m_colors_in_pos);

                                    colors_in_pos = max(colors_in_pos, colors_in_t_pos);
                                }

                                colorsfile_in_t.close();

                                break;
                            }

                            colors_in_t.seekg(pos_f_cs[l_i]);

                            l_i -= nb_pos_shared_cs;

                            readColorSets(color_sets + (l_i * block_sz), colors_in_t, min(block_sz, sz_cs - (l_i * block_sz)));
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();

            colorsfile_in.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
            colors_in.rdbuf(colorsfile_in.rdbuf());
            colors_in.seekg(colors_in_pos);
        }

        if (pos_f_cs != nullptr) delete[] pos_f_cs;
    }

    for (size_t i = 0, sz, pos; (i < overflow_sz) && colors_in.good(); ++i){

        km.read(colors_in);

        colors_in.read(reinterpret_cast<char*>(&sz), sizeof(size_t));
        colors_in.read(reinterpret_cast<char*>(&pos), sizeof(size_t));

        overflow.insert({{km, sz}, pos});
    }

    const bool ret = colors_in.good();

    colorsfile_in.close();

    return ret;
}

template<>
inline bool DataStorage<void>::read(const string& filename_colors, const size_t nb_threads, const bool verbose) {

    if (verbose) cout << endl << "DataStorage::read(): Reading color sets from disk" << endl;

    FILE* fp = fopen(filename_colors.c_str(), "rb");

    if (fp == NULL) {

        cerr << "DataStorage::read(): Could not open file " << filename_colors << " for reading color sets" << endl;
        return false;
    }
    else fclose(fp);

    auto readSharedColorSets = [](UnitigColors::SharedUnitigColors* shared_color_sets, istream& colors_in, const size_t sz){

        for (size_t i = 0; (i < sz) && colors_in.good(); ++i){

            if (shared_color_sets[i].first.read(colors_in)){

                colors_in.read(reinterpret_cast<char*>(&(shared_color_sets[i].second)), sizeof(size_t));
            }
        }
    };

    auto readColorSets = [](UnitigColors* color_sets, istream& colors_in, const size_t sz){

        for (size_t i = 0; (i < sz) && colors_in.good(); ++i) color_sets[i].read(colors_in);
    };

    Kmer km;

    size_t format_version, overflow_sz, nb_colors;

    ifstream colorsfile_in;
    istream colors_in(nullptr);

    colorsfile_in.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
    colors_in.rdbuf(colorsfile_in.rdbuf());

    clear();

    //Read the file format version number
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&format_version), sizeof(size_t));
    //Read number of different seeds for hash function
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_seeds), sizeof(size_t));
    //Read number of colors in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_colors), sizeof(size_t));
    //Read number of color sets in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_cs), sizeof(size_t));
    //Read number of elements allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&sz_cs), sizeof(size_t));
    //Write number of SharedUnitigColors allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&sz_shared_cs), sizeof(size_t));
    //Read number of (kmer, color set) overflowing
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&overflow_sz), sizeof(size_t));

    if (nb_seeds >= 256){

        cerr << "DataStorage::read(): Does not support more than 255 hash seeds" << endl;
        return false;
    }

    const size_t sz_unitig_cs_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

    overflow = unordered_map<pair<Kmer, size_t>, size_t>(overflow_sz);

    color_sets = new UnitigColors[sz_cs];
    shared_color_sets = new UnitigColors::SharedUnitigColors[sz_shared_cs];
    unitig_cs_link = new atomic<uint64_t>[sz_unitig_cs_link];

    //Read the hash function seeds of the graph
    colors_in.read(reinterpret_cast<char*>(seeds), nb_seeds * sizeof(uint64_t));

    if (format_version == 1){

        for (size_t i = 0; (i < nb_colors) && colors_in.good(); ++i){
            //Read the hash function seeds of the graph
            color_names.push_back(string());
            getline(colors_in, color_names[i]);
        }

        for (uint64_t i = 0, e; (i != sz_unitig_cs_link) && colors_in.good(); ++i){

            colors_in.read(reinterpret_cast<char*>(&e), sizeof(uint64_t));
            unitig_cs_link[i] = e;
        }

        readSharedColorSets(shared_color_sets, colors_in, sz_shared_cs);
        readColorSets(color_sets, colors_in, sz_cs);
    }
    else {

        size_t block_sz = 0;

        streampos* pos_f_cs = nullptr;

        //Read the hash function seeds of the graph
        if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&block_sz), sizeof(size_t));

        const size_t nb_pos_shared_cs = (sz_shared_cs / block_sz) + static_cast<size_t>((sz_shared_cs % block_sz) != 0);
        const size_t nb_pos_cs = (sz_cs / block_sz) + static_cast<size_t>((sz_cs % block_sz) != 0);
        const size_t pos_f_cs_sz = nb_pos_shared_cs + nb_pos_cs;

        if (pos_f_cs_sz != 0){

            pos_f_cs = new streampos[pos_f_cs_sz];

            if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(pos_f_cs), pos_f_cs_sz * sizeof(streampos));
        }

        for (size_t i = 0; (i < nb_colors) && colors_in.good(); ++i){
            //Read the hash function seeds of the graph
            color_names.push_back(string());
            getline(colors_in, color_names[i]);
        }

        for (uint64_t i = 0, e; (i != sz_unitig_cs_link) && colors_in.good(); ++i){

            colors_in.read(reinterpret_cast<char*>(&e), sizeof(uint64_t));
            unitig_cs_link[i] = e;
        }

        if ((nb_threads == 1) || (pos_f_cs_sz == 0)) {

            readSharedColorSets(shared_color_sets, colors_in, sz_shared_cs);
            readColorSets(color_sets, colors_in, sz_cs);
        }
        else {

            streampos colors_in_pos = colors_in.tellg();

            colorsfile_in.close();

            mutex m_colors_in_pos;

            vector<thread> workers; // need to keep track of threads so we can join them

            std::atomic<size_t> i;

            i = 0;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        ifstream colorsfile_in_t;
                        istream colors_in_t(nullptr);

                        colorsfile_in_t.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
                        colors_in_t.rdbuf(colorsfile_in_t.rdbuf());

                        while (true) {

                            const size_t l_i = i++;

                            if (l_i >= nb_pos_shared_cs){

                                const streampos colors_in_t_pos = colors_in_t.tellg();

                                {
                                    unique_lock<mutex> lock(m_colors_in_pos);

                                    colors_in_pos = max(colors_in_pos, colors_in_t_pos);
                                }

                                colorsfile_in_t.close();

                                break;
                            }

                            colors_in_t.seekg(pos_f_cs[l_i]);
                            readSharedColorSets(shared_color_sets + (l_i * block_sz), colors_in_t, min(block_sz, sz_shared_cs - (l_i * block_sz)));
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();

            workers.clear();

            i = nb_pos_shared_cs;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        ifstream colorsfile_in_t;
                        istream colors_in_t(nullptr);

                        colorsfile_in_t.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
                        colors_in_t.rdbuf(colorsfile_in_t.rdbuf());

                        while (true) {

                            size_t l_i = i++;

                            if (l_i >= pos_f_cs_sz){

                                const streampos colors_in_t_pos = colors_in_t.tellg();

                                {
                                    unique_lock<mutex> lock(m_colors_in_pos);

                                    colors_in_pos = max(colors_in_pos, colors_in_t_pos);
                                }

                                colorsfile_in_t.close();

                                break;
                            }

                            colors_in_t.seekg(pos_f_cs[l_i]);

                            l_i -= nb_pos_shared_cs;

                            readColorSets(color_sets + (l_i * block_sz), colors_in_t, min(block_sz, sz_cs - (l_i * block_sz)));
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();

            colorsfile_in.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
            colors_in.rdbuf(colorsfile_in.rdbuf());
            colors_in.seekg(colors_in_pos);
        }

        if (pos_f_cs != nullptr) delete[] pos_f_cs;
    }

    for (size_t i = 0, sz, pos; (i < overflow_sz) && colors_in.good(); ++i){

        km.read(colors_in);

        colors_in.read(reinterpret_cast<char*>(&sz), sizeof(size_t));
        colors_in.read(reinterpret_cast<char*>(&pos), sizeof(size_t));

        overflow.insert({{km, sz}, pos});
    }

    const bool ret = colors_in.good();

    colorsfile_in.close();

    return ret;
}

#endif
