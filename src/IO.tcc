#ifndef BIFROST_IO_CDBG_TCC
#define BIFROST_IO_CDBG_TCC

template<typename U, typename G>
bool CompactedDBG<U, G>::write(const string& output_filename, const size_t nb_threads, const bool GFA_output, const bool write_meta_file, const bool verbose) const {

    if (invalid){

        cerr << "CompactedDBG::write(): Graph is invalid and cannot be written to disk" << endl;
        return false;
    }

    if (nb_threads <= 0){

        cerr << "CompactedDBG::write(): Number of threads cannot be less than 0" << endl;
        return false;
    }

    if (nb_threads > std::thread::hardware_concurrency()){

        cerr << "CompactedDBG::write(): Number of threads cannot exceed " << std::thread::hardware_concurrency() << "threads" << endl;
        return false;
    }

    bool write_success = true;

    {
        if (verbose) cout << endl << "CompactedDBG::write(): Writing graph to disk" << endl;

        const string out = output_filename + (GFA_output ? ".gfa" : ".fasta");

        FILE* fp = fopen(out.c_str(), "w");

        if (fp == NULL) {

            cerr << "CompactedDBG::write(): Could not open file " << out << " for writing graph" << endl;

            return false;
        }
        else {

            fclose(fp);

            if (std::remove(out.c_str()) != 0) cerr << "CompactedDBG::write(): Could not remove temporary file " << out << endl;
        }

        write_success = (GFA_output ? writeGFA(out, nb_threads) : writeFASTA(out));
    }

    if (write_success && write_meta_file) {

        if (verbose) cout << endl << "CompactedDBG::write(): Writing meta file to disk" << endl;

        const string out = output_filename + ".meta.bfg";

        FILE* fp = fopen(out.c_str(), "w");

        if (fp == NULL) {

            cerr << "CompactedDBG::write(): Could not open file " << out << " for writing meta file" << endl;
            
            return false;
        }
        else {

            fclose(fp);

            if (std::remove(out.c_str()) != 0) cerr << "CompactedDBG::write(): Could not remove temporary file " << out << endl;
        }

        write_success = writeBinaryMeta(out, checksum(), nb_threads);
    }

    return write_success;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::read(const string& input_graph_filename, const size_t nb_threads, const bool verbose){

    if (verbose) cout << endl << "CompactedDBG::read(): Reading graph from disk" << endl;

    const int format = FileParser::getFileFormat(input_graph_filename.c_str());

    if (format == -1){

        cerr << "CompactedDBG::read(): Input graph file " << input_graph_filename << " does not exist, is ill-formed or is not a valid graph file format." << endl;

        return false;
    }
    else if ((format != 0) && (format != 2) && (format != 3)){

        cerr << "CompactedDBG::read(): Input graph file must be in FASTA, GFA or GRAPH.BFG format." << endl;

        return false;
    }

    if (format == 0) { // FASTA input

        const int k = k_;
        const int g = g_;

        clear();

        {
            KmerStream_Build_opt kms_opt;

            kms_opt.threads = nb_threads;
            kms_opt.verbose = verbose;
            kms_opt.k = k;
            kms_opt.g = g;
            kms_opt.q = 0;

            kms_opt.files.push_back(input_graph_filename);

            KmerStream kms(kms_opt);

            MinimizerIndex hmap_min_unitigs_tmp(max(1UL, kms.MinimizerF0()) * 1.05);

            hmap_min_unitigs = std::move(hmap_min_unitigs_tmp);
        }

        setKmerGmerLength(k, g);

        makeGraphFromFASTA(input_graph_filename, nb_threads);
    }
    else if (format == 2){ // GFA format

        FILE* fp = fopen(input_graph_filename.c_str(), "r");

        if (fp == NULL) {

            cerr << "CompactedDBG::read(): Could not open file " << input_graph_filename << " for reading graph" << endl;
            return false;
        }

        fclose(fp);

        char buffer[4096];

        ifstream graphfile_in(input_graph_filename);
        istream graph_in(graphfile_in.rdbuf());

        graph_in.getline(buffer, 4096); // Read and discard header
        graphfile_in.close();

        const string header(buffer);

        if (header[0] != 'H'){

            cerr << "CompactedDBG::read(): An error occurred while reading input GFA file." << endl;
            return false;
        }

        stringstream hs(&buffer[2]); // Skip the first 2 char. of the line "H\t"
        string sub;

        int k = k_;
        int g = g_;

        while (hs.good()){ // Split line based on tabulation

            getline(hs, sub, '\t');

            const string tag = sub.substr(0, 5);

            if (tag == "KL:Z:") k = atoi(sub.c_str() + 5);
            else if (tag == "ML:Z:") g = atoi(sub.c_str() + 5);
        }

        clear();

        {
            KmerStream_Build_opt kms_opt;

            kms_opt.threads = nb_threads;
            kms_opt.verbose = verbose;
            kms_opt.k = k;
            kms_opt.g = g;
            kms_opt.q = 0;

            kms_opt.files.push_back(input_graph_filename);

            KmerStream kms(kms_opt);

            MinimizerIndex hmap_min_unitigs_tmp(max(1UL, kms.MinimizerF0()) * 1.05);
            hmap_min_unitigs = std::move(hmap_min_unitigs_tmp);
        }

        setKmerGmerLength(k, g);

        if (!invalid) makeGraphFromGFA(input_graph_filename, nb_threads);
        if (verbose) cout << endl << "CompactedDBG::read(): Finished reading graph from disk" << endl;

        return !invalid;
    }

    // Set coverages
    {
        setFullCoverage(1);

        for (auto& unitig : *this) unitig.setFullCoverage();
    }

    invalid = false;

    if (verbose) cout << endl << "CompactedDBG::read(): Finished reading graph from disk" << endl;

    return true;
}

/*template<typename U, typename G>
bool CompactedDBG<U, G>::read(const string& input_graph_filename, const string& input_meta_filename, const size_t nb_threads, const bool verbose){

    if (verbose) cout << endl << "CompactedDBG::read(): Reading graph from disk" << endl;

    const int format_graph = FileParser::getFileFormat(input_graph_filename.c_str());
    const int format_meta = FileParser::getFileFormat(input_meta_filename.c_str());

    if (format_graph == -1){

        cerr << "CompactedDBG::read(): Input graph file " << input_graph_filename << " does not exist, is ill-formed or is not a valid graph file format." << endl;

        return false;
    }
    else if ((format_graph != 0) && (format_graph != 2) && (format_graph != 3)){

        cerr << "CompactedDBG::read(): Input graph file must be in FASTA, GFA or GRAPH.BFG format." << endl;

        return false;
    }

    if (format_meta != 4) {

         cerr << "CompactedDBG::read(): Input meta file " << input_meta_filename << " does not exist, is ill-formed or is not a valid meta file format." << endl;

        return false;
    }

    if (format == 0) { // FASTA input

        const int k = k_;
        const int g = g_;

        clear();

        {
            KmerStream_Build_opt kms_opt;

            kms_opt.threads = nb_threads;
            kms_opt.verbose = verbose;
            kms_opt.k = k;
            kms_opt.g = g;
            kms_opt.q = 0;

            kms_opt.files.push_back(input_graph_filename);

            KmerStream kms(kms_opt);

            MinimizerIndex hmap_min_unitigs_tmp(max(1UL, kms.MinimizerF0()) * 1.05);

            hmap_min_unitigs = std::move(hmap_min_unitigs_tmp);
        }

        setKmerGmerLength(k, g);

        makeGraphFromFASTA(input_graph_filename, nb_threads);
    }
    else if (format == 2){ // GFA format

        FILE* fp = fopen(input_graph_filename.c_str(), "r");

        if (fp == NULL) {

            cerr << "CompactedDBG::read(): Could not open file " << input_graph_filename << " for reading graph" << endl;
            return false;
        }

        fclose(fp);

        char buffer[4096];

        ifstream graphfile_in(input_graph_filename);
        istream graph_in(graphfile_in.rdbuf());

        graph_in.getline(buffer, 4096); // Read and discard header
        graphfile_in.close();

        const string header(buffer);

        if (header[0] != 'H'){

            cerr << "CompactedDBG::read(): An error occurred while reading input GFA file." << endl;
            return false;
        }

        stringstream hs(&buffer[2]); // Skip the first 2 char. of the line "H\t"
        string sub;

        int k = k_;
        int g = g_;

        while (hs.good()){ // Split line based on tabulation

            getline(hs, sub, '\t');

            const string tag = sub.substr(0, 5);

            if (tag == "KL:Z:") k = atoi(sub.c_str() + 5);
            else if (tag == "ML:Z:") g = atoi(sub.c_str() + 5);
        }

        clear();

        {
            KmerStream_Build_opt kms_opt;

            kms_opt.threads = nb_threads;
            kms_opt.verbose = verbose;
            kms_opt.k = k;
            kms_opt.g = g;
            kms_opt.q = 0;

            kms_opt.files.push_back(input_graph_filename);

            KmerStream kms(kms_opt);

            MinimizerIndex hmap_min_unitigs_tmp(max(1UL, kms.MinimizerF0()) * 1.05);
            hmap_min_unitigs = std::move(hmap_min_unitigs_tmp);
        }

        setKmerGmerLength(k, g);

        if (!invalid) makeGraphFromGFA(input_graph_filename, nb_threads);
        if (verbose) cout << endl << "CompactedDBG::read(): Finished reading graph from disk" << endl;

        return !invalid;
    }

    // Set coverages
    {
        setFullCoverage(1);

        for (auto& unitig : *this) unitig.setFullCoverage();
    }

    invalid = false;

    if (verbose) cout << endl << "CompactedDBG::read(): Finished reading graph from disk" << endl;

    return true;
}*/

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<!is_void, void>::type CompactedDBG<U, G>::writeGFA_sequence_(GFA_Parser& graph, KmerHashTable<size_t>& idmap) const {

    size_t labelA = 1;

    for (const auto& unitig : *this){

        const string seq(unitig.referenceUnitigToString());

        graph.write_sequence(std::to_string(labelA), seq.size(), seq, unitig.getData()->serialize(unitig));

        if (unitig.isAbundant) idmap.insert(Kmer(unitig.referenceUnitigToString().c_str()), labelA);

        ++labelA;
    }
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<is_void, void>::type CompactedDBG<U, G>::writeGFA_sequence_(GFA_Parser& graph, KmerHashTable<size_t>& idmap) const {

    size_t labelA = 1;

    for (const auto& unitig : *this){

        const string seq(unitig.referenceUnitigToString());

        graph.write_sequence(std::to_string(labelA), seq.size(), seq, "");

        if (unitig.isAbundant) idmap.insert(Kmer(unitig.referenceUnitigToString().c_str()), labelA);

        ++labelA;
    }
}

// It is very important to write unitigs to disk in the same following order:
// 1 - All unitigs with length > k
// 2 - All unitigs with length == k which do not have abundant minimizers
// 3 - All unitigs with length == k which have abundant minimizers
// The binary graph file is written in that order
// and the checksum is stored in the meta file is computed for that order
template<typename U, typename G>
bool CompactedDBG<U, G>::writeGFA(const string& fn, const size_t nb_threads) const {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = km_unitigs.size();

    size_t labelA, labelB, id = v_unitigs_sz + v_kmers_sz + 1;

    const string header_tag("BV:Z:" + string(BFG_VERSION) + "\t" + "KL:Z:" + to_string(k_) + "\t" + "ML:Z:" + to_string(g_));

    KmerHashTable<size_t> idmap(h_kmers_ccov.size());

    GFA_Parser graph(fn);

    graph.open_write(1, header_tag);

    writeGFA_sequence_<is_void<U>::value>(graph, idmap);

    if (nb_threads == 1){

        for (labelA = 1; labelA <= v_unitigs_sz; labelA++) {

            const Unitig<U>* unitig = v_unitigs[labelA - 1];
            const Kmer head = unitig->getSeq().getKmer(0);
            const Kmer tail = unitig->getSeq().getKmer(unitig->length() - k_);

            const vector<const_UnitigMap<U, G>> pred = findPredecessors(head, true);
            const vector<const_UnitigMap<U, G>> succ = findSuccessors(tail, 4, true);

            for (const auto& unitig : pred){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, 0, k_-1, false, slabelB, pos, pos + k_ - 1, !unitig.strand);
                }
            }

            for (const auto& unitig : succ){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, unitig.size - k_ + 1, unitig.size, true, slabelB, pos, pos + k_ - 1, unitig.strand);
                }
            }
        }

        for (labelA = v_unitigs_sz + 1; labelA <= v_kmers_sz + v_unitigs_sz; labelA++) {

            const Kmer km_unitig = km_unitigs.getKmer(labelA - v_unitigs_sz - 1);

            const vector<const_UnitigMap<U, G>> pred = findPredecessors(km_unitig, true);
            const vector<const_UnitigMap<U, G>> succ = findSuccessors(km_unitig, 4, true);

            for (const auto& unitig : pred){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, 0, k_-1, false, slabelB, pos, pos + k_ - 1, !unitig.strand);
                }
            }

            for (const auto& unitig : succ){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, unitig.size - k_ + 1, unitig.size, true, slabelB, pos, pos + k_ - 1, unitig.strand);
                }
            }
        }

        for (KmerHashTable<size_t>::iterator it = idmap.begin(); it != idmap.end(); it++) {

            labelA = *it;

            const vector<const_UnitigMap<U, G>> pred = findPredecessors(it.getKey(), true);
            const vector<const_UnitigMap<U, G>> succ = findSuccessors(it.getKey(), 4, true);

            for (const auto& unitig : pred){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, 0, k_-1, false, slabelB, pos, pos + k_ - 1, !unitig.strand);
                }
            }

            for (const auto& unitig : succ){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, unitig.size - k_ + 1, unitig.size, true, slabelB, pos, pos + k_ - 1, unitig.strand);
                }
            }
        }
    }
    else {

        const size_t chunk_size = 1024;

        auto worker_v_unitigs = [v_unitigs_sz, &idmap, this](const size_t labelA_start, const size_t labelA_end,
                                                             vector<pair<pair<size_t, bool>, pair<size_t, bool>>>* v_out){

            // We need to deal with the tail of long unitigs
            for (size_t labelA = labelA_start; labelA < labelA_end; ++labelA) {

                const Unitig<U>* unitig = v_unitigs[labelA - 1];

                const Kmer head = unitig->getSeq().getKmer(0);
                const Kmer tail = unitig->getSeq().getKmer(unitig->length() - k_);

                const vector<const_UnitigMap<U, G>> pred = this->findPredecessors(head, true);
                const vector<const_UnitigMap<U, G>> succ = this->findSuccessors(tail, 4, true);

                for (const auto& um : pred) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(labelA, false), make_pair(labelB, !um.strand)));
                    }
                }

                for (const auto& um : succ) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(labelA, true), make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        auto worker_v_kmers = [v_unitigs_sz, &idmap, this](const size_t labelA_start, const size_t labelA_end,
                                                             vector<pair<pair<size_t, bool>, pair<size_t, bool>>>* v_out){

            // We need to deal with the tail of long unitigs
            for (size_t labelA = labelA_start; labelA < labelA_end; ++labelA) {

                const Kmer km_unitig = km_unitigs.getKmer(labelA - v_unitigs_sz - 1);

                const vector<const_UnitigMap<U, G>> pred = this->findPredecessors(km_unitig, true);
                const vector<const_UnitigMap<U, G>> succ = this->findSuccessors(km_unitig, 4, true);

                for (const auto& um : pred) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(labelA, false), make_pair(labelB, !um.strand)));
                    }
                }

                for (const auto& um : succ) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(labelA, true), make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        auto worker_v_abundant = [v_unitigs_sz, chunk_size, &idmap, this](  KmerHashTable<size_t>::iterator* l_it,
                                                                            vector<pair<pair<size_t, bool>, pair<size_t, bool>>>* v_out){

            KmerHashTable<size_t>::iterator& it = *l_it;

            // We need to deal with the tail of long unitigs
            for (size_t i = 0; (it != idmap.end()) && (i < chunk_size); ++i, ++it) {

                const vector<const_UnitigMap<U, G>> pred = this->findPredecessors(it.getKey(), true);
                const vector<const_UnitigMap<U, G>> succ = this->findSuccessors(it.getKey(), 4, true);

                for (const auto& um : pred) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(*it, false), make_pair(labelB, !um.strand)));
                    }
                }

                for (const auto& um : succ) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(*it, true), make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        {
            atomic<size_t> label(1);

            vector<vector<pair<pair<size_t, bool>, pair<size_t, bool>>>> v_out(nb_threads);
            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mutex_file;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        while (true) {

                            const size_t old_labelA = label.fetch_add(chunk_size);

                            if (old_labelA <= v_unitigs_sz){

                                if (old_labelA + chunk_size <= v_unitigs_sz) worker_v_unitigs(old_labelA, old_labelA + chunk_size, &v_out[t]);
                                else worker_v_unitigs(old_labelA, v_unitigs_sz + 1, &v_out[t]);

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    for (const auto& p : v_out[t]){

                                        const string slabelA = std::to_string(p.first.first);
                                        const string slabelB = std::to_string(p.second.first);

                                        graph.write_edge(slabelA, 0, k_-1, p.first.second, slabelB, 0, k_-1, p.second.second);
                                    }
                                }

                                v_out[t].clear();
                            }
                            else return;
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }

        {
            const size_t v_kmers_unitigs_sz = v_kmers_sz + v_unitigs_sz;

            atomic<size_t> label(v_unitigs_sz + 1);

            vector<vector<pair<pair<size_t, bool>, pair<size_t, bool>>>> v_out(nb_threads);
            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mutex_file;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        while (true) {

                            const size_t old_labelA = label.fetch_add(chunk_size);

                            if (old_labelA <= v_kmers_unitigs_sz){

                                if (old_labelA + chunk_size <= v_kmers_unitigs_sz) worker_v_kmers(old_labelA, old_labelA + chunk_size, &v_out[t]);
                                else worker_v_kmers(old_labelA, v_kmers_unitigs_sz + 1, &v_out[t]);

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    for (const auto& p : v_out[t]){

                                        const string slabelA = std::to_string(p.first.first);
                                        const string slabelB = std::to_string(p.second.first);

                                        graph.write_edge(slabelA, 0, k_-1, p.first.second, slabelB, 0, k_-1, p.second.second);
                                    }
                                }

                                v_out[t].clear();
                            }
                            else return;
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }

        {
            KmerHashTable<size_t>::iterator it = idmap.begin(), it_end = idmap.end();

            vector<vector<pair<pair<size_t, bool>, pair<size_t, bool>>>> v_out(nb_threads);
            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mutex_file, mutex_it;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        KmerHashTable<size_t>::iterator l_it;

                        bool stop;

                        while (true) {

                            {
                                unique_lock<mutex> lock(mutex_it);

                                l_it = it;

                                for (size_t i = 0; (it != it_end) && (i < chunk_size); ++i, ++it){}

                                stop = (l_it == it_end) && (it == it_end);
                            }

                            if (!stop){

                                worker_v_abundant(&l_it, &v_out[t]);

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    for (const auto& p : v_out[t]){

                                        const string slabelA = std::to_string(p.first.first);
                                        const string slabelB = std::to_string(p.second.first);

                                        graph.write_edge(slabelA, 0, k_-1, p.first.second, slabelB, 0, k_-1, p.second.second);
                                    }
                                }

                                v_out[t].clear();
                            }
                            else return;
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }
    }

    graph.close();

    return true;
}

// It is very important to write unitigs to disk in the same following order:
// 1 - All unitigs with length > k
// 2 - All unitigs with length == k which do not have abundant minimizers
// 3 - All unitigs with length == k which have abundant minimizers
// The binary graph file is written in that order
// and the checksum is stored in the meta file is computed for that order
template<typename U, typename G>
bool CompactedDBG<U, G>::writeFASTA(const string& fn) const {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = km_unitigs.size();

    size_t i = 0;

    ofstream graphfile;
    ostream graph(0);

    graphfile.open(fn.c_str());
    graph.rdbuf(graphfile.rdbuf());
    //graph.sync_with_stdio(false);

    for (size_t j = 0; !graph.fail() && (j < v_unitigs_sz); ++j, ++i) graph << ">" << i << "\n" << v_unitigs[j]->getSeq().toString() << "\n";
    for (size_t j = 0; !graph.fail() && (j < v_kmers_sz); ++j, ++i) graph << ">" << i << "\n" << km_unitigs.getKmer(j).toString() << "\n";

    for (typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(); !graph.fail() && (it != h_kmers_ccov.end()); ++it, ++i) {

        graph << ">" << i << "\n" << it.getKey().toString() << "\n";
    }

    const bool write_success = !graph.fail();

    graphfile.close();

    return write_success;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinary(const string& fn) {

    if ((fn.length() == 0) || !check_file_exists(fn)) return false;

    ifstream infile;
    istream in(0);

    infile.open(fn.c_str());
    in.rdbuf(infile.rdbuf());

    return readBinary(in);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinary(istream& in) {

    //return (!in.fail() && readBinaryGraph(in) && readBinaryMeta(in, checksum()));

    if (!in.fail()) {

    	const pair<uint64_t, bool> p_readSuccess_checksum = readBinaryGraph(in);

    	if (p_readSuccess_checksum.second) return readBinaryMeta(in, p_readSuccess_checksum.first);
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinaryMeta(const string& fn, const uint64_t checksum) {

    if ((fn.length() == 0) || !check_file_exists(fn)) return false;

    ifstream infile;
    istream in(0);

    infile.open(fn.c_str());
    in.rdbuf(infile.rdbuf());

    return readBinaryMeta(in, checksum);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinaryMetaHead(const string& fn, size_t& file_format_version, size_t& v_unitigs_sz, size_t& km_unitigs_sz,
											size_t& h_kmers_ccov_sz, size_t& hmap_min_unitigs_sz, uint64_t& read_checksum) const {

    if ((fn.length() == 0) || !check_file_exists(fn)) return false;

    ifstream infile;
    istream in(0);

    infile.open(fn.c_str());
    in.rdbuf(infile.rdbuf());

    return readBinaryMetaHead(in, file_format_version, v_unitigs_sz, km_unitigs_sz, h_kmers_ccov_sz, hmap_min_unitigs_sz, read_checksum);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinaryMetaHead(istream& in, size_t& file_format_version, size_t& v_unitigs_sz, size_t& km_unitigs_sz,
											size_t& h_kmers_ccov_sz, size_t& hmap_min_unitigs_sz, uint64_t& read_checksum) const {

	if (in.fail()) return false;

    in.read(reinterpret_cast<char*>(&file_format_version), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&read_checksum), sizeof(uint64_t));

    in.read(reinterpret_cast<char*>(&v_unitigs_sz), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&km_unitigs_sz), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&h_kmers_ccov_sz), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&hmap_min_unitigs_sz), sizeof(size_t));

    return !in.fail();
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinaryMeta(istream& in, const uint64_t checksum) {

    bool read_success = !in.fail();

    // 0 - Write file format version, checksum and number of minimizers
    if (read_success) {

        size_t file_format_version = 0, v_unitigs_sz = 0, km_unitigs_sz = 0, h_kmers_ccov_sz = 0, hmap_min_unitigs_sz = 0;
        uint64_t read_checksum = 0;

        read_success = readBinaryMetaHead(in, file_format_version, v_unitigs_sz, km_unitigs_sz, h_kmers_ccov_sz, hmap_min_unitigs_sz, read_checksum);

        if (!read_success || ((file_format_version >> 32) != BFG_METABIN_FORMAT_HEADER)) return false;
        if (!read_success || (read_checksum != checksum)) return false;

        hmap_min_unitigs = MinimizerIndex(hmap_min_unitigs_sz);
    }

    if (read_success) {

        size_t nb_bmp_unitigs = 0;

        in.read(reinterpret_cast<char*>(&nb_bmp_unitigs), sizeof(size_t));

        vector<BitContainer> v_bmp_unitigs(nb_bmp_unitigs);

        for (size_t i = 0; read_success && (i < nb_bmp_unitigs); ++i) read_success = v_bmp_unitigs[i].read(in);

        if (read_success && !v_bmp_unitigs.empty() && !v_bmp_unitigs.front().isEmpty()) {

            size_t unitig_id = 0, tot_unitig_len = 0, curr_unitig_len = v_unitigs[0]->getSeq().size() - g_ + 1;

            for (size_t i = 0; i < nb_bmp_unitigs; ++i) {

                const size_t id_bmp = i << 32;

                for (const auto pos_bmp : v_bmp_unitigs[i]) {

                    const size_t pos = id_bmp + pos_bmp;

                    while (pos >= (tot_unitig_len + curr_unitig_len)) {

                        ++unitig_id;
                        tot_unitig_len += curr_unitig_len;

                        curr_unitig_len = v_unitigs[unitig_id]->getSeq().size() - g_ + 1;
                    }

                    const size_t relative_pos = pos - tot_unitig_len;
                    const size_t pos_id_unitig = (unitig_id << 32) | relative_pos;

                    const Minimizer minz_rep = v_unitigs[unitig_id]->getSeq().getMinimizer(relative_pos).rep();

                    std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);

                    packed_tiny_vector& v = p.first.getVector();
                    uint8_t& flag_v = p.first.getVectorSize();

                    flag_v = v.push_back(pos_id_unitig, flag_v);
                }

                v_bmp_unitigs[i].clear();
            }
        }
    }

    if (read_success) {

        size_t nb_bmp_km = 0;

        in.read(reinterpret_cast<char*>(&nb_bmp_km), sizeof(size_t));

        vector<BitContainer> v_bmp_km(nb_bmp_km);

        for (size_t i = 0; (i < nb_bmp_km) && read_success; ++i) read_success = v_bmp_km[i].read(in);

        if (read_success && !v_bmp_km.empty() && !v_bmp_km.front().isEmpty()) {

            const size_t km_glen = k_ - g_ + 1;

            for (size_t i = 0; i < nb_bmp_km; ++i) {

                const size_t id_bmp = i << 32;

                for (const auto pos_bmp : v_bmp_km[i]) {

                    const size_t pos = id_bmp + pos_bmp;

                    const size_t km_id = pos / km_glen;
                    const size_t km_pos = pos % km_glen;

                    const size_t pos_id_unitig = (km_id << 32) | MASK_CONTIG_TYPE | km_pos;

                    const Minimizer minz_rep = km_unitigs.getMinimizer(km_id, km_pos).rep();

                    std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);

                    packed_tiny_vector& v = p.first.getVector();
                    uint8_t& flag_v = p.first.getVectorSize();

                    flag_v = v.push_back(pos_id_unitig, flag_v);
                }

                v_bmp_km[i].clear();
            }
        }
    }

    if (read_success) {

        size_t nb_minz_abundant_overcrowded = 0;

        in.read(reinterpret_cast<char*>(&nb_minz_abundant_overcrowded), sizeof(size_t));

        read_success = !in.fail();

        if (read_success) {

            vector<BitContainer> v_bmp_abundant((nb_minz_abundant_overcrowded >> 32) + 1);
            vector<BitContainer> v_bmp_overcrowded((nb_minz_abundant_overcrowded >> 32) + 1);

            for (size_t i = 0; (i < v_bmp_abundant.size()) && read_success; ++i) read_success = v_bmp_abundant[i].read(in);
            for (size_t i = 0; (i < v_bmp_overcrowded.size()) && read_success; ++i) read_success = v_bmp_overcrowded[i].read(in);

            for (size_t i = 0; (i < nb_minz_abundant_overcrowded) && read_success; ++i) {

                Minimizer minz_rep;

                const bool isAbundant = v_bmp_abundant[i >> 32].contains(i & 0x00000000ffffffffULL);
                const bool isOvercrowded = v_bmp_overcrowded[i >> 32].contains(i & 0x00000000ffffffffULL);

                size_t pos_id_unitig = MASK_CONTIG_ID;

                pos_id_unitig |= (static_cast<size_t>(!isOvercrowded) - 1) & MASK_CONTIG_TYPE;

                read_success = minz_rep.read(in);

                if (isAbundant && read_success) {

                    uint32_t count_abundant = 0;

                    in.read(reinterpret_cast<char*>(&count_abundant), sizeof(uint32_t));

                    pos_id_unitig |= (static_cast<size_t>(count_abundant));
                    read_success = !in.fail();
                }

                if (read_success) {

                    std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);

                    packed_tiny_vector& v = p.first.getVector();
                    uint8_t& flag_v = p.first.getVectorSize();

                    flag_v = v.push_back(pos_id_unitig, flag_v);
                }
            }
        }
    }

    if (read_success) {

        size_t nb_mismatch_minz = 0;

        in.read(reinterpret_cast<char*>(&nb_mismatch_minz), sizeof(size_t));

        read_success = !in.fail();

		for (size_t i = 0; (i < nb_mismatch_minz) && read_success; ++i) {

            Minimizer minz_rep;

            size_t pos_id_unitig = 0;

            minz_rep.read(in);

            in.read(reinterpret_cast<char*>(&pos_id_unitig), sizeof(size_t));

            read_success = !in.fail();

            if (read_success) {

                std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);

                packed_tiny_vector& v = p.first.getVector();
                uint8_t& flag = p.first.getVectorSize();
                size_t v_sz = v.size(flag);

                if (v_sz == 0 || ((v(v_sz-1, flag) & MASK_CONTIG_ID) != MASK_CONTIG_ID)) flag = v.push_back(pos_id_unitig, flag);
                else flag = v.insert(pos_id_unitig, v_sz-1, flag);
            }
        }
    }

    return read_success;
}

template<typename U, typename G>
pair<uint64_t, bool> CompactedDBG<U, G>::readBinaryGraph(const string& fn) {

    if ((fn.length() == 0) || !check_file_exists(fn)) return false;

    ifstream infile;
    istream in(0);

    infile.open(fn.c_str());
    in.rdbuf(infile.rdbuf());

    return readBinaryGraph(in);
}

/*template<typename U, typename G>
bool CompactedDBG<U, G>::readGraphFromMetaFASTA(const string& graph_fn, const string& meta_fn) {

    size_t graph_file_id = 0;

    FastqFile ff(vector<string>(1, fn));

    string seq;

	size_t file_format_version = 0, v_unitigs_sz = 0, km_unitigs_sz = 0, h_kmers_ccov_sz = 0, hmap_min_unitigs_sz = 0;
	uint64_t read_checksum = 0;

    if (!readBinaryMetaHead(meta_fn, file_format_version, v_unitigs_sz, km_unitigs_sz, h_kmers_ccov_sz, hmap_min_unitigs_sz, read_checksum)) return false;

    while (ff.read_next(seq, graph_file_id) != -1) {

    	if (seq.length())
    }
}*/

template<typename U, typename G>
pair<uint64_t, bool> CompactedDBG<U, G>::readBinaryGraph(istream& in) {

    bool read_success = !in.fail();

    uint64_t graph_checksum = 0;

    clear();

    if (read_success) {

        size_t file_format_version = 0;
        int rk = 0, rg = 0;

        in.read(reinterpret_cast<char*>(&file_format_version), sizeof(size_t));
        in.read(reinterpret_cast<char*>(&rk), sizeof(int));
        in.read(reinterpret_cast<char*>(&rg), sizeof(int));

        read_success = (read_success && !in.fail());

        if (read_success) {

	        const size_t k = static_cast<size_t>(rk);
	        const size_t g = static_cast<size_t>(rg);

	        graph_checksum = wyhash(&k, sizeof(size_t), 0, _wyp);
	        graph_checksum = wyhash(&g, sizeof(size_t), graph_checksum, _wyp);
        }

        if ((file_format_version >> 32) != BFG_GRAPHBIN_FORMAT_HEADER) return {graph_checksum, false};
        if (read_success) *this = CompactedDBG<U, G>(rk, rg);
        if (invalid) return {graph_checksum, false};
    }

    if (read_success) {

        size_t v_unitigs_sz = 0;

        in.read(reinterpret_cast<char*>(&v_unitigs_sz), sizeof(size_t));

        read_success = !in.fail();

        if (read_success) {

            v_unitigs.reserve(v_unitigs_sz);

            for (size_t i = 0; (i < v_unitigs_sz) && read_success; ++i) {

                CompressedSequence cs;
                CompressedCoverage cc;

                Unitig<U>* unitig;

                read_success = cs.read(in);
                graph_checksum = cs.hash(graph_checksum);
                cc = CompressedCoverage(cs.size() - k_ + 1, false);
                unitig = new Unitig<U>(move(cs), move(cc));

                v_unitigs.push_back(unitig);
            }
        }
    }

    if (read_success) {

    	read_success = km_unitigs.read(in);

    	for (size_t i = 0; i < km_unitigs.size(); ++i) graph_checksum = km_unitigs.getKmer(i).hash(graph_checksum);
    }

    if (read_success) {

        const CompressedCoverage cc(1, false);

        size_t h_kmers_ccov_sz = 0;

        in.read(reinterpret_cast<char*>(&h_kmers_ccov_sz), sizeof(size_t));

        read_success = !in.fail();

        if (read_success) {

            h_kmers_ccov.reserve(h_kmers_ccov_sz);

            for (size_t i = 0; read_success && (i < h_kmers_ccov_sz); ++i) {

                Kmer km;

                read_success = km.read(in);
                graph_checksum = km.hash(graph_checksum);

                h_kmers_ccov.insert(km, cc);
            }
        }
    }

    if (read_success) {

        setFullCoverage(1);

        for (auto& unitig : *this) unitig.setFullCoverage();
    }

    return {graph_checksum, read_success};
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinary(const string& fn, const size_t nb_threads) const {

    if (fn.length() == 0) return false;

    ofstream outfile;
    ostream out(0);

    outfile.open(fn.c_str());
    out.rdbuf(outfile.rdbuf());

    return writeBinary(out, nb_threads);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinary(ostream& out, const size_t nb_threads) const {

    return (!out.fail() && writeBinaryGraph(out, nb_threads) && writeBinaryMeta(out, checksum(), nb_threads));
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinaryGraph(const string& fn, const size_t nb_threads) const {

    if (fn.length() == 0) return false;

    ofstream outfile;
    ostream out(0);

    outfile.open(fn.c_str());
    out.rdbuf(outfile.rdbuf());

    return writeBinaryGraph(out, nb_threads);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinaryGraph(ostream& out, const size_t nb_threads) const {

    bool write_success = !out.fail();

    // 0- Write file format version, k-mer and g-mer lengths
    if (write_success) {

        const size_t fileformat_version = (static_cast<size_t>(BFG_GRAPHBIN_FORMAT_HEADER) << 32) | static_cast<size_t>(BFG_GRAPHBIN_FORMAT_VERSION);

        out.write(reinterpret_cast<const char*>(&fileformat_version), sizeof(size_t));
        out.write(reinterpret_cast<const char*>(&k_), sizeof(int));
        out.write(reinterpret_cast<const char*>(&g_), sizeof(int));

        write_success = (write_success && !out.fail());
    }

    // 1 - Write unitigs longer than k
    if (write_success) {

        const size_t v_unitigs_sz = v_unitigs.size();

        out.write(reinterpret_cast<const char*>(&v_unitigs_sz), sizeof(size_t));

        write_success = !out.fail();

        for (size_t i = 0; (i < v_unitigs_sz) && write_success; ++i) write_success = v_unitigs[i]->getSeq().write(out);
    }

    // 2 - Write unitigs with length k for which minimizers are not over abundant
    if (write_success) write_success = km_unitigs.write(out);

    // 3 - Write unitigs with length k for which minimizers are over abundant
    if (write_success) {

        const size_t h_kmers_ccov_sz = h_kmers_ccov.size();

        out.write(reinterpret_cast<const char*>(&h_kmers_ccov_sz), sizeof(size_t));

        write_success = !out.fail();

        for (typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(); (it != h_kmers_ccov.end()) && write_success; ++it) write_success = it.getKey().write(out);
    }

    return (write_success && !out.fail());
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinaryMeta(const string& fn, const uint64_t checksum, const size_t nb_threads) const {

    if (fn.length() == 0) return false;

    ofstream outfile;
    ostream out(0);

    outfile.open(fn.c_str());
    out.rdbuf(outfile.rdbuf());

    return writeBinaryMeta(out, checksum, nb_threads);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinaryMeta(ostream& out, const uint64_t checksum, const size_t nb_threads) const {

    bool write_success = !out.fail();

    // 0 - Write file format version, checksum and number of minimizers
    if (write_success) {

        const size_t fileformat_version = (static_cast<size_t>(BFG_METABIN_FORMAT_HEADER) << 32) | static_cast<size_t>(BFG_METABIN_FORMAT_VERSION);

        const size_t v_unitigs_sz = v_unitigs.size();
		const size_t km_unitigs_sz = km_unitigs.size();
		const size_t h_kmers_ccov_sz = h_kmers_ccov.size();
        const size_t hmap_min_unitigs_sz = hmap_min_unitigs.size();

        out.write(reinterpret_cast<const char*>(&fileformat_version), sizeof(size_t)); // Write header for binary meta file, including file format version
        out.write(reinterpret_cast<const char*>(&checksum), sizeof(uint64_t)); // Write graph checksum

        out.write(reinterpret_cast<const char*>(&v_unitigs_sz), sizeof(size_t)); // Write number of unitigs with length > k
        out.write(reinterpret_cast<const char*>(&km_unitigs_sz), sizeof(size_t)); // Write number of unitigs with length == k and non-abundant minimizer
        out.write(reinterpret_cast<const char*>(&h_kmers_ccov_sz), sizeof(size_t)); // Write number of unitigs with length == k and abundant minimizer

        out.write(reinterpret_cast<const char*>(&hmap_min_unitigs_sz), sizeof(size_t)); // Write number of minimizers

        write_success = (write_success && !out.fail());
    }

    if (write_success) {

        const size_t v_unitigs_sz = v_unitigs.size();
        const size_t nb_block_unitigs = max(static_cast<size_t>((v_unitigs_sz + 15) / 16) + 1, static_cast<size_t>(1));

        vector<size_t> v_block_len_unitigs(nb_block_unitigs, 0);

        for (size_t i = 0; i < v_unitigs_sz; ++i) v_block_len_unitigs[(i >> 4) + 1] += v_unitigs[i]->getSeq().size() - g_ + 1;
        for (size_t i = 1; i < nb_block_unitigs; ++i) v_block_len_unitigs[i] += v_block_len_unitigs[i-1];

        {
            const size_t nb_bmp_unitigs = (v_block_len_unitigs.back() >> 32) + 1;
            const size_t nb_bmp_km_short = ((km_unitigs.size() * (k_ - g_ + 1)) >> 32) + 1;

            vector<BitContainer> v_bmp_unitigs(nb_bmp_unitigs), v_bmp_km_short(nb_bmp_km_short), v_bmp_abundant(1), v_bmp_overcrowded(1);

            MinimizerIndex::const_iterator it = hmap_min_unitigs.begin();
            MinimizerIndex::const_iterator ite = hmap_min_unitigs.end();

            size_t nb_minz_abundant_overcrowded = 0;
            size_t nb_mismatch_minz = 0;

            {
                while (it != ite) { // Annotate in bitmap the position of every minimizer in the unitigs

                    const packed_tiny_vector& v = it.getVector();
                    const uint8_t flag_v = it.getVectorSize();
                    const int v_sz = v.size(flag_v);

                    const Minimizer& minz_key = it.getKey();

                    for (size_t i = 0; i < v_sz; ++i){

                        const size_t unitig_idx = v(i, flag_v);
                        const size_t unitig_id = unitig_idx >> 32;
                        const size_t unitig_pos = unitig_idx & MASK_CONTIG_POS;

                        const bool isShort = static_cast<bool>(unitig_idx & MASK_CONTIG_TYPE);

                        if (unitig_id == RESERVED_ID) {

                            const size_t id_bmp = nb_minz_abundant_overcrowded >> 32;
                            const size_t pos_bmp = nb_minz_abundant_overcrowded & 0x00000000ffffffffULL;

                            while (id_bmp > (v_bmp_abundant.size() - 1)) {

                                v_bmp_abundant.push_back(BitContainer());
                                v_bmp_overcrowded.push_back(BitContainer());
                            }

                            if (unitig_pos != 0) v_bmp_abundant[id_bmp].add(pos_bmp);
                            if (isShort) v_bmp_overcrowded[id_bmp].add(pos_bmp);

                            ++nb_minz_abundant_overcrowded;
                        }
                        else if (isShort) {

                        	const Minimizer minz = km_unitigs.getMinimizer(unitig_id, unitig_pos).rep();

                        	if (minz == minz_key) {

	                            const size_t pos = (unitig_id * (k_ - g_ + 1)) + unitig_pos;

	                            v_bmp_km_short[pos >> 32].add(pos & 0x00000000ffffffffULL);
                        	}
                        	else ++nb_mismatch_minz;
                        }
                        else {

                        	const Minimizer minz = v_unitigs[unitig_id]->getSeq().getMinimizer(unitig_pos).rep();

                        	if (minz == minz_key) {

	                            size_t pos = v_block_len_unitigs[unitig_id >> 4] + unitig_pos;

	                            for (size_t j = (unitig_id & 0xfffffffffffffff0ULL); j < unitig_id; ++j) pos += v_unitigs[j]->getSeq().size() - g_ + 1;

	                            v_bmp_unitigs[pos >> 32].add(pos & 0x00000000ffffffffULL);
                    		}
                    		else ++nb_mismatch_minz;
                        }
                    }

                    ++it;
                }

                {
                    for (auto& bmp : v_bmp_unitigs) bmp.runOptimize();

                    out.write(reinterpret_cast<const char*>(&nb_bmp_unitigs), sizeof(size_t));

                    write_success = !out.fail();

                    for (size_t i = 0; (i < nb_bmp_unitigs) && write_success; ++i) write_success = v_bmp_unitigs[i].write(out);

                    v_bmp_unitigs.clear();
                }

                {
                    for (auto& bmp : v_bmp_km_short) bmp.runOptimize();

                    out.write(reinterpret_cast<const char*>(&nb_bmp_km_short), sizeof(size_t));

                    write_success = !out.fail();

                    for (size_t i = 0; (i < nb_bmp_km_short) && write_success; ++i) write_success = v_bmp_km_short[i].write(out);

                    v_bmp_km_short.clear();
                }
            }

            if (write_success) {

                it = hmap_min_unitigs.begin();
                ite = hmap_min_unitigs.end();

                out.write(reinterpret_cast<const char*>(&nb_minz_abundant_overcrowded), sizeof(size_t)); // Pre-reserve space

                write_success = !out.fail();

                {
                    for (auto& bmp : v_bmp_abundant) bmp.runOptimize();
                    for (auto& bmp : v_bmp_overcrowded) bmp.runOptimize();

                    for (size_t i = 0; (i < v_bmp_abundant.size()) && write_success; ++i) write_success = v_bmp_abundant[i].write(out);
                    for (size_t i = 0; (i < v_bmp_overcrowded.size()) && write_success; ++i) write_success = v_bmp_overcrowded[i].write(out);

                    v_bmp_abundant.clear();
                    v_bmp_overcrowded.clear();
                }

                while ((it != ite) && write_success) { // Annotate in bitmap the position of every minimizer in the unitigs

                    const packed_tiny_vector& v = it.getVector();
                    const uint8_t flag_v = it.getVectorSize();
                    const int v_sz = v.size(flag_v);

                    const size_t unitig_idx = v(v_sz-1, flag_v);
                    const size_t unitig_id = unitig_idx >> 32;

                    const uint32_t unitig_pos = static_cast<uint32_t>(unitig_idx & MASK_CONTIG_POS);

                    if (unitig_id == RESERVED_ID) {

                        write_success = it.getKey().write(out);

                        if (write_success && (unitig_pos != 0)) {

                            out.write(reinterpret_cast<const char*>(&unitig_pos), sizeof(uint32_t)); // Pre-reserve space

                            write_success = !out.fail();
                        }
                    }

                    ++it;
                }
            }

            if (write_success) {

                it = hmap_min_unitigs.begin();
                ite = hmap_min_unitigs.end();

                out.write(reinterpret_cast<const char*>(&nb_mismatch_minz), sizeof(size_t)); // Pre-reserve space

                write_success = !out.fail();

                while ((it != ite) && write_success) { // Annotate in bitmap the position of every minimizer in the unitigs

                    const packed_tiny_vector& v = it.getVector();
                    const uint8_t flag_v = it.getVectorSize();
                    const int v_sz = v.size(flag_v);

                    const Minimizer& minz_key = it.getKey();

                    for (size_t i = 0; (i < v_sz) && write_success; ++i){

                        const size_t unitig_idx = v(i, flag_v);
                        const size_t unitig_id = unitig_idx >> 32;
                        const size_t unitig_pos = unitig_idx & MASK_CONTIG_POS;

                        const bool isShort = static_cast<bool>(unitig_idx & MASK_CONTIG_TYPE);

                        if (unitig_id != RESERVED_ID) {

                        	Minimizer minz;

                        	if (isShort) minz = km_unitigs.getMinimizer(unitig_id, unitig_pos).rep();
                        	else minz = v_unitigs[unitig_id]->getSeq().getMinimizer(unitig_pos).rep();

                        	if (minz != minz_key) {

                        		minz_key.write(out);

                        		out.write(reinterpret_cast<const char*>(&unitig_idx), sizeof(size_t)); // Pre-reserve space

                            	write_success = !out.fail();
                            }
                        }
                    }

                    ++it;
                }
            }
        }
    }

    return (write_success && !out.fail());
}

template<typename U, typename G>
void CompactedDBG<U, G>::makeGraphFromGFA(const string& fn, const size_t nb_threads) {

    size_t graph_file_id = 0;

    bool new_file_opened = false;

    GFA_Parser graph(fn);

    graph.open_read();

    GFA_Parser::GFA_line r = graph.read(graph_file_id, new_file_opened, true);

    if (nb_threads == 1){

        while ((r.first != nullptr) || (r.second != nullptr)){

            if (r.first != nullptr) addUnitig(r.first->seq, (r.first->seq.length() == k_) ? km_unitigs.size() : v_unitigs.size());

            r = graph.read(graph_file_id, new_file_opened, true);
        }
    }
    else {

        const size_t block_sz = 1024;

        std::atomic<size_t> v_kmers_sz;
        std::atomic<size_t> v_unitigs_sz;

        bool is_first = true;
        bool stop = false;

        SpinLock lck_unitig, lck_kmer;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_file;

        v_kmers_sz = 0;
        v_unitigs_sz = 0;

        hmap_min_unitigs.init_threads();

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&]{

                    vector<string> seq;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_file);

                            if (stop) return;

                            seq.clear();

                            for (size_t i = 0; (i < block_sz) && !stop; ++i){

                                if (!is_first) r = graph.read(graph_file_id, new_file_opened, true);
                                if (r.first != nullptr) seq.push_back(r.first->seq);

                                stop = ((r.first == nullptr) && (r.second == nullptr));
                                is_first = false;
                            }
                        }

                        for (const auto& s : seq) addUnitig(s, (s.length() == k_) ? v_kmers_sz++ : v_unitigs_sz++, lck_unitig, lck_kmer);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        hmap_min_unitigs.release_threads();

        moveToAbundant();
	}
}

template<typename U, typename G>
void CompactedDBG<U, G>::makeGraphFromFASTA(const string& fn, const size_t nb_threads) {

    size_t graph_file_id = 0;

    FastqFile ff(vector<string>(1, fn));

    string seq;

    if (nb_threads == 1){

        while (ff.read_next(seq, graph_file_id) != -1) addUnitig(seq, (seq.length() == k_) ? km_unitigs.size() : v_unitigs.size());
    }
    else {

        const size_t block_sz = 1024;

        bool stop = false;

        std::atomic<size_t> v_kmers_sz;
        std::atomic<size_t> v_unitigs_sz;

        SpinLock lck_unitig, lck_kmer;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_file;

        v_kmers_sz = 0;
        v_unitigs_sz = 0;

        hmap_min_unitigs.init_threads();

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&]{

                    vector<string> v_seq;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_file);

                            if (stop) return;

                            v_seq.clear();

                            for (size_t i = 0; (i < block_sz) && !stop; ++i){

                                stop = (ff.read_next(seq, graph_file_id) == -1);

                                if (!stop && !seq.empty()) v_seq.push_back(seq);
                            }
                        }

                        for (const auto& s : v_seq) addUnitig(s, (s.length() == k_) ? v_kmers_sz++ : v_unitigs_sz++, lck_unitig, lck_kmer);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        hmap_min_unitigs.release_threads();

        moveToAbundant();
    }
}

#endif