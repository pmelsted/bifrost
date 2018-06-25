#ifndef BFG_COLOREDCDBG_TCC
#define BFG_COLOREDCDBG_TCC

template<typename U>
ColoredCDBG<U>::ColoredCDBG(int kmer_length, int minimizer_length) : CompactedDBG<DataAccessor<U>, DataStorage<U>>(kmer_length, minimizer_length){

    invalid = CompactedDBG<DataAccessor<U>, DataStorage<U>>::isInvalid();
}

template<typename U>
ColoredCDBG<U>::ColoredCDBG(const ColoredCDBG& o) : CompactedDBG<DataAccessor<U>, DataStorage<U>>(o), invalid(o.invalid) {}

template<typename U>
ColoredCDBG<U>::ColoredCDBG(ColoredCDBG&& o) :  CompactedDBG<DataAccessor<U>, DataStorage<U>>(o), invalid(o.invalid) {}

template<typename U>
void ColoredCDBG<U>::clear(){

    invalid = true;

    CompactedDBG<DataAccessor<U>, DataStorage<U>>::getData()->clear();
    CompactedDBG<DataAccessor<U>, DataStorage<U>>::clear();
}

template<typename U>
void ColoredCDBG<U>::empty(){

    invalid = true;

    CompactedDBG<DataAccessor<U>, DataStorage<U>>::getData()->empty();
    CompactedDBG<DataAccessor<U>, DataStorage<U>>::empty();
}

template<typename U>
ColoredCDBG<U>& ColoredCDBG<U>::operator=(const ColoredCDBG& o) {

    CompactedDBG<DataAccessor<U>, DataStorage<U>>::operator=(o);

    invalid = o.invalid;

    return *this;
}

template<typename U>
ColoredCDBG<U>& ColoredCDBG<U>::operator=(ColoredCDBG&& o) {

    if (this != &o) {

        CompactedDBG<DataAccessor<U>, DataStorage<U>>::operator=(o);

        invalid = o.invalid;

        o.clear();
    }

    return *this;
}

template<typename U>
bool ColoredCDBG<U>::buildGraph(const CCDBG_Build_opt& opt){

    if (!invalid){

        CDBG_Build_opt opt_ = opt.getCDBG_Build_opt();

        invalid = !CompactedDBG<DataAccessor<U>, DataStorage<U>>::build(opt_);
    }
    else cerr << "ColoredCDBG::buildGraph(): Graph is invalid and cannot be built." << endl;

    return !invalid;
}

template<typename U>
bool ColoredCDBG<U>::buildColors(const CCDBG_Build_opt& opt){

    if (!invalid){

        if (opt.filename_colors_in.size() == 0){

            initUnitigColors(opt);
            buildUnitigColors(opt.nb_threads);
        }
        else invalid = !readColorSets(opt);
    }
    else cerr << "ColoredCDBG::buildColors(): Graph is invalid (maybe not built yet?) and colors cannot be mapped." << endl;

    return !invalid;
}

template<typename U>
bool ColoredCDBG<U>::write(const string& prefix_output_filename, const size_t nb_threads, const bool verbose) const {

    if (!CompactedDBG<DataAccessor<U>, DataStorage<U>>::write(prefix_output_filename, nb_threads, true, verbose)) return false; // Write graph

    return CompactedDBG<DataAccessor<U>, DataStorage<U>>::getData()->write(prefix_output_filename, nb_threads, verbose); // Write colors
}

template<typename U>
bool ColoredCDBG<U>::read(const string& prefix_input_filename, const size_t nb_threads, const bool verbose) {

    if (verbose) cout << "ColoredCDBG::read(): Reading graph." << endl;

    if (!CompactedDBG<DataAccessor<U>, DataStorage<U>>::read(prefix_input_filename + ".gfa", verbose)) return false; // Read graph

    if (verbose) cout << "ColoredCDBG::read(): Reading colors." << endl;

    if (!CompactedDBG<DataAccessor<U>, DataStorage<U>>::getData()->read(prefix_input_filename + ".bfg_colors", verbose)) return false; // Read colors

    if (verbose) cout << "ColoredCDBG::read(): Joining unitigs to their color sets." << endl;

    GFA_Parser graph(prefix_input_filename + ".gfa");

    graph.open_read();

    auto reading_function = [&graph](vector<pair<Kmer, uint8_t>>& unitig_tags, const size_t chunk_size) {

        size_t i = 0;
        size_t graph_file_id = 0;

        bool new_file_opened = false;

        GFA_Parser::GFA_line r = graph.read(graph_file_id, new_file_opened, true);

        while ((i < chunk_size) && ((r.first != nullptr) || (r.second != nullptr))){

            if (r.first != nullptr){ // It is a sequence

                if (r.first->tags.empty()){

                    cerr << "ColoredCDBG::read(): One sequence line in GFA file has no DataAccessor tag. Operation aborted." << endl;
                    return false;
                }

                size_t i = 0;

                for (; i < r.first->tags.size(); ++i){

                    if (r.first->tags[i].substr(0, 5) == "DA:Z:") break;
                }

                if (i == r.first->tags.size()){

                    cerr << "ColoredCDBG::read(): One sequence line in GFA file has no DataAccessor tag. Operation aborted." << endl;
                    return false;
                }

                unitig_tags.push_back({Kmer(r.first->seq.c_str()), atoi(r.first->tags[i].c_str() + 5)});

                ++i;
            }

            r = graph.read(graph_file_id, new_file_opened, true);
        }

        return ((r.first != nullptr) || (r.second != nullptr));
    };

    auto join_function = [this](const vector<pair<Kmer, uint8_t>>& unitig_tags) {

        for (const auto& p : unitig_tags){

            UnitigColorMap<U> ucm(this->find(p.first, true));

            if (ucm.isEmpty){

                cerr << "ColoredCDBG::read(): Internal error, operation aborted." << endl;
                cerr << "ColoredCDBG::read(): A unitig from GFA file is not found in the in-memory graph." << endl;
                cerr << "ColoredCDBG::read(): Graph from GFA file possibly incorrectly compacted." << endl;

                return false;
            }

            DataAccessor<U>* da = ucm.getData();

            *da = DataAccessor<U>(p.second);

            if (!ucm.strand){ // Unitig has been inserted in reverse-complement, need to reverse order of color sets

                UnitigColors* uc = da->getUnitigColors(ucm);

                UnitigColors r_uc = uc->reverse(ucm);

                *uc = move(r_uc);
            }
        }

        return true;
    };

    {
        const size_t chunk = 10000;

        vector<thread> workers; // need to keep track of threads so we can join them
        vector<vector<pair<Kmer, uint8_t>>> v(nb_threads);

        mutex mutex_file;

        bool file_valid_for_read = true;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_file);

                            if (!file_valid_for_read) return;

                            file_valid_for_read = reading_function(v[t], chunk);

                        }

                        join_function(v[t]);
                        v[t].clear();
                    }
                }
            );
        }

        for (auto& t : workers) t.join();
    }

    return true;
}

template<typename U>
void ColoredCDBG<U>::initUnitigColors(const CCDBG_Build_opt& opt, const size_t max_nb_hash){

    size_t last_empty_pos = 0;

    mutex mutex_cs_overflow;

    DataStorage<U>* ds = this->getData();
    DataStorage<U> new_ds = DataStorage<U>(max_nb_hash, this->size(), opt.filename_seq_in);

    *ds = move(new_ds);

    auto worker_function = [&](typename ColoredCDBG<U>::iterator a, typename ColoredCDBG<U>::iterator b){

        int i;

        uint64_t h_v, id_link_mod;

        for (auto& unitig = a; unitig != b; ++unitig){

            const Kmer head = unitig->getUnitigHead();

            for (i = 0; i < max_nb_hash; ++i){

                h_v = head.hash(ds->seeds[i]) % ds->nb_cs; // Hash to which we can possibly put our colorset for current kmer
                id_link_mod = 1ULL << (h_v & 0x3F);

                if ((ds->unitig_cs_link[h_v >> 6].fetch_or(id_link_mod) & id_link_mod) == 0) break;
            }

            if (i == max_nb_hash){ // IF we couldn't find a hash matching an unoccupied color set for current k-mer

                unique_lock<mutex> lock(mutex_cs_overflow);

                while (true){

                    id_link_mod = 1ULL << (last_empty_pos & 0x3F);

                    if ((ds->unitig_cs_link[last_empty_pos >> 6].fetch_or(id_link_mod) & id_link_mod) == 0) break;

                    last_empty_pos = ((last_empty_pos + 1) == ds->nb_cs ? 0 : last_empty_pos + 1);
                }

                ds->overflow.insert(head, last_empty_pos); // Insertion
            }

            *(unitig->getData()) = DataAccessor<U>(static_cast<uint8_t>(i == max_nb_hash ? 0 : i + 1));
        }
    };

    {
        const size_t chunk = 1000;

        vector<thread> workers; // need to keep track of threads so we can join them

        typename ColoredCDBG<U>::iterator g_a = this->begin();
        typename ColoredCDBG<U>::iterator g_b = this->end();

        mutex mutex_it;

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    typename ColoredCDBG<U>::iterator l_a, l_b;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_it);

                            if (g_a == g_b) return;

                            l_a = g_a;
                            l_b = g_a;

                            for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                            g_a = l_b;
                        }

                        worker_function(l_a, l_b);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();
    }

    //cout << "Number of unitigs not hashed is " << ds->overflow.size() << " on " << ds->nb_cs << " unitigs." << endl;
}

template<typename U>
void ColoredCDBG<U>::buildUnitigColors(const size_t nb_threads){

    DataStorage<U>* ds = this->getData();

    const int k_ = this->getK();

    const size_t nb_locks = nb_threads * 1024;
    const size_t chunk_size = 64;

    size_t prev_file_id = 0;

    bool next_file = true;

    FileParser fp(ds->color_names);

    std::atomic_flag* cs_locks = new std::atomic_flag[nb_locks];

    for (size_t i = 0; i < nb_locks; ++i) cs_locks[i].clear();

    // Main worker thread
    auto worker_function = [&](const vector<pair<string, size_t>>& v_read_color) {

        // for each input
        for (const auto& read_color : v_read_color) {

            for (KmerIterator it_km(read_color.first.c_str()), it_km_end; it_km != it_km_end; ++it_km) {

                UnitigColorMap<U> um = this->find(it_km->first);

                if (!um.isEmpty) {

                    if (um.strand || (um.dist != 0)){

                        um.len = 1 + um.lcp(read_color.first.c_str(), it_km->second + k_, um.strand ? um.dist + k_ : um.dist - 1, !um.strand);

                        if ((um.size != k_) && !um.strand) um.dist -= um.len - 1;

                        it_km += um.len - 1;
                    }

                    const uint64_t id_lock = ds->getHash(um) % nb_locks;
                    UnitigColors* uc = ds->getUnitigColors(um);

                    while (cs_locks[id_lock].test_and_set(std::memory_order_acquire)); // Set the corresponding lock

                    uc->add(um, read_color.second);

                    cs_locks[id_lock].clear(std::memory_order_release);
                }
            }
        }
    };

    size_t pos_read = k_ - 1;
    size_t len_read = 0;

    string s;

    auto reading_function = [&](vector<pair<string, size_t>>& v_read_color) {

        size_t reads_now = 0;
        size_t file_id = prev_file_id;

        while ((pos_read < len_read) && (reads_now < chunk_size)){

            pos_read -= k_ - 1;

            v_read_color.emplace_back(make_pair(s.substr(pos_read, 1000), file_id));

            pos_read += 1000;

            ++reads_now;
        }

        while (reads_now < chunk_size) {

            if (fp.read(s, file_id)) {

                len_read = s.length();
                pos_read = len_read;

                if (len_read > 1000){

                    pos_read = k_ - 1;

                    while ((pos_read < len_read) && (reads_now < chunk_size)){

                        pos_read -= k_ - 1;

                        v_read_color.emplace_back(make_pair(s.substr(pos_read, 1000), file_id));

                        pos_read += 1000;

                        ++reads_now;
                    }
                }
                else {

                    v_read_color.emplace_back(make_pair(s, file_id));
                    ++reads_now;
                }
            }
            else {

                next_file = false;

                for (auto& p : v_read_color) std::transform(p.first.begin(), p.first.end(), p.first.begin(), ::toupper);

                return true;
            }
        }

        for (auto& p : v_read_color) std::transform(p.first.begin(), p.first.end(), p.first.begin(), ::toupper);

        const bool ret = (file_id != prev_file_id);

        next_file = true;
        prev_file_id = file_id;

        return ret;
    };

    {
        bool stop = false;

        vector<thread> workers; // need to keep track of threads so we can join them
        vector<vector<pair<string, size_t>>> reads_colors(nb_threads);

        mutex mutex_file;

        size_t prev_uc_sz = getCurrentRSS();

        while (next_file){

            stop = false;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        while (true) {

                            {
                                unique_lock<mutex> lock(mutex_file);

                                if (stop) return;

                                stop = reading_function(reads_colors[t]);
                            }

                            worker_function(reads_colors[t]);
                            reads_colors[t].clear();
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();

            workers.clear();

            for (size_t t = 0; t < nb_threads; ++t) reads_colors[t].clear();

            const size_t curr_uc_sz = getCurrentRSS();

            if ((curr_uc_sz - prev_uc_sz) >= 1073741824ULL){

                const size_t chunk = 1000;

                typename ColoredCDBG<U>::iterator g_a = this->begin();
                typename ColoredCDBG<U>::iterator g_b = this->end();

                mutex mutex_it;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&, t]{

                            typename ColoredCDBG<U>::iterator l_a, l_b;

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_it);

                                    if (g_a == g_b) return;

                                    l_a = g_a;
                                    l_b = g_a;

                                    for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                                    g_a = l_b;
                                }

                                while (l_a != l_b){

                                    l_a->getData()->getUnitigColors(*l_a)->optimizeFullColors(*l_a);
                                    ++l_a;
                                }
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();

                workers.clear();

                prev_uc_sz = getCurrentRSS();
            }
        }
    }

    fp.close();

    //checkColors(ds->color_names);

    /*typedef std::unordered_map<uint64_t, pair<int64_t, size_t>> uc_unordered_map;

    uc_unordered_map u_map;

    mutex mutex_u_map;

    vector<thread> workers;

    auto add_hash_function = [&](typename ColoredCDBG<U>::iterator it_a, typename ColoredCDBG<U>::iterator it_b) {

        while (it_a != it_b) {

            const const_UnitigColorMap<U> unitig(*it_a);
            const UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
            const UnitigColors uc_full = uc->makeFullColors(unitig);
            const UnitigColors* uc_full_array = uc_full.getFullColorsPtr();

            if (uc_full_array[0].size() != 0){

                const pair<int64_t, size_t> pv(0 - static_cast<int64_t>(uc_full_array[0].getSizeInBytes()) - static_cast<int64_t>(sizeof(size_t)), 0);

                const int64_t to_add = static_cast<int64_t>(uc->getSizeInBytes());
                const int64_t to_rm = (static_cast<int64_t>(uc_full_array[1].getSizeInBytes() + 2 * sizeof(UnitigColors)));

                {
                    unique_lock<mutex> lock(mutex_u_map);

                    pair<uc_unordered_map::iterator, bool> p = u_map.insert(make_pair(uc_full_array[0].hash(), pv));

                    p.first->second.first += to_add - to_rm;
                }
            }

            ++it_a;
        }
    };

    auto add_shared_function = [&](typename ColoredCDBG<U>::iterator it_a, typename ColoredCDBG<U>::iterator it_b) {

        while (it_a != it_b) {

            const UnitigColorMap<U> unitig(*it_a);

            UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
            UnitigColors uc_full = uc->makeFullColors(unitig);
            UnitigColors* uc_full_array = uc_full.getFullColorsPtr();

            if (uc_full_array[0].size() != 0){

                uc_unordered_map::const_iterator it = u_map.find(uc_full_array[0].hash());

                if (it->second.first > 0){ // If there is some sharing

                    const size_t id_shared = it->second.second;
                    const uint64_t id_lock = id_shared % nb_locks;

                    bool move_full = false;

                    while (cs_locks[id_lock].test_and_set(std::memory_order_acquire)); // Set the corresponding lock

                    if (ds->shared_color_sets[id_shared].second == 0){

                        ds->shared_color_sets[id_shared].first = move(uc_full_array[0]);
                        ds->shared_color_sets[id_shared].second = 0;

                        uc_full_array[0] = ds->shared_color_sets[id_shared];
                        move_full = true;
                    }
                    else if (uc_full_array[0] == ds->shared_color_sets[id_shared]) {

                        uc_full_array[0] = ds->shared_color_sets[id_shared];
                        move_full = true;
                    }

                    cs_locks[id_lock].clear(std::memory_order_release);

                    if (move_full) *uc = move(uc_full);
                }
            }

            ++it_a;
        }
    };

    {
        const size_t chunk = 1000;

        typename ColoredCDBG<U>::iterator g_a = this->begin();
        typename ColoredCDBG<U>::iterator g_b = this->end();

        mutex mutex_it;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    typename ColoredCDBG<U>::iterator l_a, l_b;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_it);

                            if (g_a == g_b) return;

                            l_a = g_a;
                            l_b = g_a;

                            for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                            g_a = l_b;
                        }

                        add_hash_function(l_a, l_b);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        workers.clear();
    }

    for (auto& p : u_map){

        if (p.second.first > 0){

            p.second.second = ds->sz_shared_cs;
            ++(ds->sz_shared_cs);
        }
    }

    ds->shared_color_sets = new UnitigColors::SharedUnitigColors[ds->sz_shared_cs];

    {
        const size_t chunk = 1000;

        typename ColoredCDBG<U>::iterator g_a = this->begin();
        typename ColoredCDBG<U>::iterator g_b = this->end();

        mutex mutex_it;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    typename ColoredCDBG<U>::iterator l_a, l_b;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_it);

                            if (g_a == g_b) return;

                            l_a = g_a;
                            l_b = g_a;

                            for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                            g_a = l_b;
                        }

                        add_shared_function(l_a, l_b);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        workers.clear();
    }*/

    delete[] cs_locks;
}

template<typename U>
string ColoredCDBG<U>::getColorName(const size_t color_id) const {

    if (invalid){

        cerr << "ColoredCDBG::getColorName(): Graph is invalid or colors are not yet mapped to unitigs." << endl;
        return string("ColoredCDBG::getColorName(): Graph is invalid or colors are not yet mapped to unitigs.");
    }

    const DataStorage<U>* ds = this->getData();

    if (color_id >= ds->color_names.size()){

        cerr << "ColoredCDBG::getColorName(): Color ID " << color_id << " is invalid, graph only has " <<
        ds->color_names.size() << " colors." << endl;

        return string("ColoredCDBG::getColorName(): Given color ID is invalid.");
    }

    return ds->color_names[color_id];
}

template<typename U>
void ColoredCDBG<U>::checkColors(const vector<string>& filename_seq_in) const {

    cout << "ColoredCDBG::checkColors(): Start" << endl;

    size_t file_id = 0;

    string s;

    KmerHashTable<tiny_vector<size_t, 1>> km_h;

    FastqFile FQ(filename_seq_in);

    while (FQ.read_next(s, file_id) >= 0){

        for (KmerIterator it_km(s.c_str()), it_km_end; it_km != it_km_end; ++it_km) {

            pair<KmerHashTable<tiny_vector<size_t, 1>>::iterator, bool> it = km_h.insert(it_km->first.rep(), tiny_vector<size_t, 1>());

            tiny_vector<size_t, 1>& tv = *(it.first);

            const size_t id = file_id / 64;

            while (tv.size() < (id + 1)) tv.push_back(0);

            tv[id] |= (1ULL << (file_id % 64));
        }
    }

    FQ.close();

    cout << "ColoredCDBG::checkColors(): All k-mers in the hash table with their colors" << endl;

    for (typename KmerHashTable<tiny_vector<size_t, 1>>::const_iterator it_km = km_h.begin(), it_km_end = km_h.end(); it_km != it_km_end; ++it_km){

        const Kmer km = it_km.getKey();
        const const_UnitigColorMap<U> ucm = this->find(km);

        if (ucm.isEmpty){

            cerr << "ColoredCDBG::checkColors(): K-mer " << km.toString() << " is not found in the graph" << endl;
            exit(1);
        }

        const UnitigColors* cs = ucm.getData()->getUnitigColors(ucm);

        if (cs == nullptr){

            cerr << "ColoredCDBG::checkColors(): K-mer " << km.toString() << " has no color set associated" << endl;
            exit(1);
        }

        const tiny_vector<size_t, 1>& tv = *it_km;
        const size_t tv_nb_max_elem = tv.size() * 64;

        for (size_t i = 0; i < std::min(filename_seq_in.size(), tv_nb_max_elem); ++i){

            const bool color_pres_graph = cs->contains(ucm, i);
            const bool color_pres_hasht = ((tv[i/64] >> (i%64)) & 0x1) == 0x1;

            if (color_pres_graph != color_pres_hasht){

                cerr << "ColoredCDBG::checkColors(): Current color is " << i << ": " << filename_seq_in[i] << endl;
                cerr << "ColoredCDBG::checkColors(): K-mer " << km.toString() << " for color " << i << ": " << filename_seq_in[i] << endl;
                cerr << "ColoredCDBG::checkColors(): Full unitig: " << ucm.toString() << endl;
                cerr << "ColoredCDBG::checkColors(): Present in graph: " << color_pres_graph << endl;
                cerr << "ColoredCDBG::checkColors(): Present in hash table: " << color_pres_hasht << endl;

                exit(1);
            }
        }
    }

    cout << "ColoredCDBG::checkColors(): Checked all colors of all k-mers: everything is fine" << endl;
    cout << "ColoredCDBG::checkColors(): Number of k-mers in the graph: " << km_h.size() << endl;
}

#endif
