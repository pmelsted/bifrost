#ifndef BFG_COLOREDCDBG_TCC
#define BFG_COLOREDCDBG_TCC

template<typename U>
ColoredCDBG<U>::ColoredCDBG(int kmer_length, int minimizer_length) : CompactedDBG<DataAccessor<U>, DataStorage<U>>(kmer_length, minimizer_length){

    invalid = this->isInvalid();
}

template<typename U>
ColoredCDBG<U>::ColoredCDBG(const ColoredCDBG& o) : CompactedDBG<DataAccessor<U>, DataStorage<U>>(o), invalid(o.invalid) {}

template<typename U>
ColoredCDBG<U>::ColoredCDBG(ColoredCDBG&& o) :  CompactedDBG<DataAccessor<U>, DataStorage<U>>(o), invalid(o.invalid) {}

template<typename U>
void ColoredCDBG<U>::clear(){

    invalid = true;

    CompactedDBG<DataAccessor<U>, DataStorage<U>>::clear();
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
bool ColoredCDBG<U>::build(const CCDBG_Build_opt& opt){

    if (!invalid){

        CDBG_Build_opt opt_ = opt.getCDBG_Build_opt();

        invalid = !CompactedDBG<DataAccessor<U>, DataStorage<U>>::build(opt_);
    }
    else cerr << "ColoredCDBG::build(): Graph is invalid and cannot be built." << endl;

    return !invalid;
}

template<typename U>
bool ColoredCDBG<U>::mapColors(const CCDBG_Build_opt& opt){

    if (!invalid){

        if (opt.filename_colors_in.size() == 0){

            initColorSets(opt);
            buildColorSets(opt.nb_threads);
        }
        else invalid = !readColorSets(opt);
    }
    else cerr << "ColoredCDBG::mapColors(): Graph is invalid (maybe not built yet?) and colors cannot be mapped." << endl;

    return !invalid;
}

template<typename U>
bool ColoredCDBG<U>::write(const string prefix_output_filename, const size_t nb_threads, const bool verbose){

    if (!CompactedDBG<DataAccessor<U>, DataStorage<U>>::write(prefix_output_filename, nb_threads, true, verbose)) return false;

    return this->getData()->write(prefix_output_filename, nb_threads, verbose);
}

template<typename U>
void ColoredCDBG<U>::initColorSets(const CCDBG_Build_opt& opt, const size_t max_nb_hash){

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

                h_v = head.hash(ds->seeds[i]) % ds->nb_color_sets; // Hash to which we can possibly put our colorset for current kmer
                id_link_mod = 1ULL << (h_v & 0x3F);

                if ((ds->unitig_cs_link[h_v >> 6].fetch_or(id_link_mod) & id_link_mod) == 0) break;
            }

            if (i == max_nb_hash){ // IF we couldn't find a hash matching an unoccupied color set for current k-mer

                unique_lock<mutex> lock(mutex_cs_overflow);

                while (true){

                    id_link_mod = 1ULL << (last_empty_pos & 0x3F);

                    if ((ds->unitig_cs_link[last_empty_pos >> 6].fetch_or(id_link_mod) & id_link_mod) == 0) break;

                    last_empty_pos = ((last_empty_pos + 1) == ds->nb_color_sets ? 0 : last_empty_pos + 1);
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

    cout << "Number of unitigs not hashed is " << ds->overflow.size() << " on " << ds->nb_color_sets << " unitigs." << endl;
}

template<typename U>
void ColoredCDBG<U>::buildColorSets(const size_t nb_threads){

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
        }
    }

    fp.close();

    delete[] cs_locks;

    //checkColors(ds->color_names);

    /*size_t nb_unitigs = 0;
    size_t nb_neighbors_modified = 0;

    size_t nb_km_colors = 0;
    size_t nb_km_full_colors = 0;

    size_t cs_sz_before = 0;
    size_t cs_sz_after = 0;

    size_t cs_card_before = 0;
    size_t cs_card_after = 0;

    vector<pair<UnitigColorMap<U>, pair<size_t, size_t>>> v_unitigs;

    vector<TinyBitmap> v_new_uc;

    unordered_set<Kmer, KmerHash> unitigs_modified;

    auto sort_unitigs = [](const pair<UnitigColorMap<U>, pair<size_t, size_t>>& p1,
                           const pair<UnitigColorMap<U>, pair<size_t, size_t>>& p2) {

        if (p1.second.first == p2.second.first) return (p1.second.second < p2.second.second);
        return (p1.second.first < p2.second.first);

        //return ((p1.second.first * p1.second.second) < (p2.second.first * p2.second.second));
    };

    for (auto& unitig : *this){

        const DataAccessor<U>* da = unitig.getData();
        const UnitigColors* uc = da->getUnitigColors(unitig);

        pair<size_t, size_t> degree(0, 0); // (edge degree, full color degree)

        degree.first = unitig.getPredecessors().cardinality() + unitig.getSuccessors().cardinality();

        for (UnitigColors::const_iterator it = uc->begin(); it != uc->end(); it.nextColor(unitig.len)){

            if (uc->size(unitig, it->getColorID(unitig.len)) == unitig.len) ++(degree.second);
        }

        v_unitigs.push_back(make_pair(unitig, degree));
    }

    sort(v_unitigs.begin(), v_unitigs.end(), sort_unitigs);

    for (auto& p : v_unitigs){

        UnitigColorMap<U>& unitig = p.first;

        UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);

        // If unitig has not been used already and if it has at least one successor or predecessor
        if ((p.second.first != 0) && (unitigs_modified.find(unitig.getUnitigHead()) == unitigs_modified.end())){

            const Kmer unitig_head = unitig.getUnitigHead();

            size_t nb_neighbor_inter = 0;

            vector<size_t> id_unitig, id_unitig_inter;

            for (UnitigColors::const_iterator it = uc->begin(); it != uc->end(); it.nextColor(unitig.len)){

                const size_t color_id = it->getColorID(unitig.len);

                if (uc->size(unitig, color_id) == unitig.len) id_unitig.push_back(color_id);
            }

            id_unitig_inter = id_unitig;

            for (const auto& pred : unitig.getPredecessors()){

                if ((pred.getUnitigHead() != unitig_head) && (unitigs_modified.find(pred.getUnitigHead()) == unitigs_modified.end())){

                    const UnitigColors* uc_p = pred.getData()->getUnitigColors(pred);

                    vector<size_t> id_unitig_p, intersection;

                    for (UnitigColors::const_iterator it_p = uc_p->begin(); it_p != uc_p->end(); it_p.nextColor(pred.len)){

                        const size_t color_id = it_p->getColorID(pred.len);

                        if (uc_p->size(pred, color_id) == pred.len) id_unitig_p.push_back(color_id);
                    }

                    set_intersection(id_unitig_inter.begin(), id_unitig_inter.end(), id_unitig_p.begin(), id_unitig_p.end(), back_inserter(intersection));

                    id_unitig_inter = intersection;

                    ++nb_neighbor_inter;
                }
            }

            for (const auto& succ : unitig.getSuccessors()){

                if ((succ.getUnitigHead() != unitig_head) && (unitigs_modified.find(succ.getUnitigHead()) == unitigs_modified.end())){

                    const UnitigColors* uc_p = succ.getData()->getUnitigColors(succ);

                    vector<size_t> id_unitig_p, intersection;

                    for (UnitigColors::const_iterator it_p = uc_p->begin(); it_p != uc_p->end(); it_p.nextColor(succ.len)){

                        const size_t color_id = it_p->getColorID(succ.len);

                        if (uc_p->size(succ, color_id) == succ.len) id_unitig_p.push_back(color_id);
                    }

                    set_intersection(id_unitig_inter.begin(), id_unitig_inter.end(), id_unitig_p.begin(), id_unitig_p.end(), back_inserter(intersection));

                    id_unitig_inter = intersection;

                    ++nb_neighbor_inter;
                }
            }

            if ((nb_neighbor_inter != 0) && (id_unitig_inter.size() >= 0.5 * id_unitig.size())){

                unitigs_modified.insert(unitig.getUnitigHead());

                const UnitigMapBase umb(0, 1, Kmer::k, true);

                UnitigColors new_uc;

                for (const auto color_id : id_unitig_inter){

                    uc->remove(unitig, color_id);
                    new_uc.add(umb, color_id);
                }

                v_new_uc.push_back(move(new_uc));

                for (auto& pred : unitig.getPredecessors()){

                    if (unitigs_modified.find(pred.getUnitigHead()) == unitigs_modified.end()){

                        UnitigColors* uc_p = pred.getData()->getUnitigColors(pred);

                        for (const auto color_id : id_unitig_inter) uc_p->remove(pred, color_id);

                        unitigs_modified.insert(pred.getUnitigHead());

                        ++nb_neighbors_modified;
                    }
                }

                for (auto& succ : unitig.getSuccessors()){

                    if (unitigs_modified.find(succ.getUnitigHead()) == unitigs_modified.end()){

                        UnitigColors* uc_p = succ.getData()->getUnitigColors(succ);

                        for (const auto color_id : id_unitig_inter) uc_p->remove(succ, color_id);

                        unitigs_modified.insert(succ.getUnitigHead());

                        ++nb_neighbors_modified;
                    }
                }

                ++nb_unitigs;
            }
        }
    }*/
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
