#include "CompactedDBG.hpp"

CompactedDBG::CompactedDBG(int kmer_length, int minimizer_length) : k_(kmer_length), g_(minimizer_length), invalid(false) {

    if (k_ == 0){

        cerr << "CompactedDBG::CompactedDBG(): Length k of k-mers cannot be less than or equal to 0" << endl;
        clear();

        return;
    }

    if (k_ >= MAX_KMER_SIZE){

        cerr << "CompactedDBG::CompactedDBG(): Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        clear();

        return;
    }

    if (g_ == 0){

        cerr << "CompactedDBG::CompactedDBG(): Length g of minimizers cannot be less than or equal to 0" << endl;
        clear();

        return;
    }

    if (g_ >= MAX_KMER_SIZE){

        cerr << "CompactedDBG::CompactedDBG(): Length g of minimizers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        clear();

        return;
    }

    Kmer::set_k(k_);
    Minimizer::set_g(g_);
}

CompactedDBG::~CompactedDBG() {

    empty();
}

void CompactedDBG::clear(){

    k_ = 0;
    g_ = 0;

    invalid = true;

    empty();
}

void CompactedDBG::empty(){

    for (auto unitig : v_unitigs) delete unitig;

    v_unitigs.clear();
    v_kmers.clear();

    hmap_min_unitigs.clear();
    h_kmers_ccov.clear();

    bf.clear();
}

bool CompactedDBG::build(const vector<string>& fastx_filename_in, const size_t nb_unique_kmers, const size_t nb_non_unique_kmers,
                        const bool reference_mode, const size_t nb_threads, const size_t nb_bits_unique_kmers_bf,
                        const size_t nb_bits_non_unique_kmers_bf, const string& inFilenameBBF, const string& outFilenameBBF,
                        const size_t read_chunksize, const size_t unitig_size, const bool verbose){

    size_t max_threads = std::thread::hardware_concurrency();

    if (invalid){

        cerr << "CompactedDBG::build(): Graph is invalid and cannot be built" << endl;
        return false;
    }

    if (nb_threads <= 0){

        cerr << "CompactedDBG::build(): Number of threads cannot be less than or equal to 0" << endl;
        return false;
    }

    if (nb_threads > max_threads){

        cerr << "CompactedDBG::build(): Number of threads cannot exceed " << max_threads << "threads" << endl;
        return false;
    }

    if (nb_bits_unique_kmers_bf <= 0){

        cerr << "CompactedDBG::build(): Number of Bloom filter bits per unique k-mer cannot be less than or equal to 0" << endl;
        return false;
    }

    if (!reference_mode && (nb_bits_non_unique_kmers_bf <= 0)){

        cerr << "CompactedDBG::build(): Number of Bloom filter bits per non unique k-mer cannot be less than or equal to 0" << endl;
        return false;
    }

    if (!reference_mode && (nb_non_unique_kmers > nb_unique_kmers)){

        cerr << "CompactedDBG::build(): The estimated number of non unique k-mers ";
        cerr << "cannot be greater than the estimated number of unique k-mers" << endl;
        return false;
    }

    if (read_chunksize <= 0){

        cerr << "CompactedDBG::build(): Chunk size of reads to share among threads cannot be less than or equal to 0" << endl;
        return false;
    }

    if (unitig_size <= 0){

        cerr << "CompactedDBG::build(): Maximum unitig size cannot be less than or equal to 0" << endl;
        return false;
    }

    if (fastx_filename_in.size() == 0){

        cerr << "CompactedDBG::build(): Number of FASTA/FASTQ files in input cannot be less than or equal to 0" << endl;
        return false;
    }

    if (inFilenameBBF.length() != 0){

        FILE* fp = fopen(inFilenameBBF.c_str(), "rb");

        if (fp == NULL) {

            cerr << "CompactedDBG::build(): Could not open input Blocked Bloom filter: " << inFilenameBBF << endl;
            return false;
        }
        else fclose(fp);
    }

    if (outFilenameBBF.length() != 0){

        FILE* fp = fopen(outFilenameBBF.c_str(), "wb");

        if (fp == NULL) {

            cerr << "CompactedDBG::build(): Could not open file for writing output Blocked Bloom filter: " << outFilenameBBF << endl;
            return false;
        }
        else fclose(fp);
    }

    for (vector<string>::const_iterator it = fastx_filename_in.begin(); it != fastx_filename_in.end(); it++){

        FILE* fp = fopen(it->c_str(), "r");

        if (fp == NULL) {

            cerr << "CompactedDBG::build(): Could not open input FASTA/FASTQ file " << *it << endl;
            return false;
        }
        else fclose(fp);
    }

    empty();

    bool construct_finished;

    if (inFilenameBBF.length() != 0){

        FILE* fp = fopen(inFilenameBBF.c_str(), "rb");

        if (fp == NULL) {

            cerr << "CompactedDBG::build(): Could not open input Blocked Bloom filter: " << inFilenameBBF << endl;
            return false;
        }

        construct_finished = bf.ReadBloomFilter(fp);

        fclose(fp);
    }
    else construct_finished = filter(fastx_filename_in, nb_unique_kmers, nb_non_unique_kmers, reference_mode,
                                        nb_threads, nb_bits_unique_kmers_bf, nb_bits_non_unique_kmers_bf,
                                        read_chunksize, unitig_size, verbose);

    if (construct_finished){

        if (outFilenameBBF.length() != 0){

            FILE* fp = fopen(outFilenameBBF.c_str(), "wb");

            if (fp == NULL) {

                cerr << "CompactedDBG::build(): Could not open file for writing output Blocked Bloom filter: " << outFilenameBBF << endl;
                return false;
            }

            bf.WriteBloomFilter(fp);

            fclose(fp);
        }

        // Construction step
        construct_finished = construct(fastx_filename_in, nb_threads, read_chunksize, verbose);

        bf.clear();
    }

    return construct_finished;
}

bool CompactedDBG::simplify(const bool delete_short_isolated_unitigs, const bool clip_short_tips, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::simplify(): Graph is invalid and cannot be simplified" << endl;
        return false;
    }

    if (delete_short_isolated_unitigs || clip_short_tips){

        if (verbose) cout << endl << "CompactedDBG::simplify(): Removing isolated unitigs and/or clipping tips" << endl;

        vector<Kmer> v_joins;
        size_t joined = 0;

        size_t removed = removeUnitigs(delete_short_isolated_unitigs, clip_short_tips, v_joins);

        if (clip_short_tips) joined = joinAllUnitigs(&v_joins);

        v_joins.clear();

        if (verbose){

            cout << "CompactedDBG::simplify(): After: " << size() << " unitigs" << endl;
            cout << "CompactedDBG::simplify(): Removed " << removed << " unitigs" << endl;
            cout << "CompactedDBG::simplify(): Joined " << joined << " unitigs" << endl;
        }

        return true;
    }

    return false;
}

bool CompactedDBG::join(const UnitigMap& um, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::join(): Graph is invalid and no unitigs can be joined in it" << endl;
        return false;
    }

    if (verbose) cout << "CompactedDBG::join(): Start joining unitigs" << endl;

    vector<Kmer> v_km;

    v_km.push_back(um.getHead());
    if (!um.isShort && !um.isAbundant) v_km.push_back(um.getTail());

    size_t joined = joinAllUnitigs(&v_km);

    if (verbose) cout << "CompactedDBG::join(): " << joined << " unitigs have been joined" << endl;

    return true;
}

bool CompactedDBG::join(const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::join(): Graph is invalid and no unitigs can be joined in it" << endl;
        return false;
    }

    if (verbose) cout << "CompactedDBG::join(): Start joining all unitigs of the graph" << endl;

    vector<Kmer> v_km;

    size_t joined = joinAllUnitigs(&v_km);

    if (verbose) cout << "CompactedDBG::join(): " << joined << " unitigs have been joined" << endl;

    return true;
}

bool CompactedDBG::write(const string output_filename, const bool verbose) const {

    if (invalid){

        cerr << "CompactedDBG::write(): Graph is invalid and cannot be written to disk" << endl;
        return false;
    }

    if (verbose) cout << endl << "CompactedDBG::write(): Writing graph to GFA file" << endl;

    string out = output_filename + ".gfa";

    FILE* fp = fopen(out.c_str(), "w");

    if (fp == NULL) {

        cerr << "CompactedDBG::write(): Could not open file " << out << " for writing graph" << endl;
        return false;
    }
    else fclose(fp);

    writeGFA(out);

    return true;
}

bool CompactedDBG::add(const string& seq, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::add(): Graph is invalid and no sequence can be added to it" << endl;
        return false;
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::add(): Input sequence length cannot be less than k = " << k_ << endl;
        return false;
    }

    size_t nxt_pos_insert_v_unitigs = v_unitigs.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t v_kmers_sz = v_kmers.size();
    size_t added = 0;

    vector<pair<int,int>> sp;
    vector<Kmer> v_joins;

    for (KmerIterator it_km(seq.c_str()), it_km_end; it_km != it_km_end; it_km++) { //non-ACGT char. are discarded

        const std::pair<Kmer, int>& p = *it_km;
        const UnitigMap cm = findUnitig(p.first, seq, p.second);

        if (cm.isEmpty){

            for (size_t i = 0; i != 4; i++){

                UnitigMap cm_bw = find(p.first.backwardBase(alpha[i]));

                if (!cm_bw.isEmpty && !cm_bw.isAbundant && !cm_bw.isShort){

                    if (cm_bw.strand) cm_bw.dist++;

                    if ((cm_bw.dist != 0) && (cm_bw.dist != cm_bw.size - k_ + 1)){

                        sp.push_back(make_pair(0, cm_bw.dist));
                        sp.push_back(make_pair(cm_bw.dist, cm_bw.size - k_ + 1));

                        splitUnitig(cm_bw.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

                        sp.clear();
                    }
                }
            }

            for (size_t i = 0; i != 4; i++){

                UnitigMap cm_fw = find(p.first.forwardBase(alpha[i]));

                if (!cm_fw.isEmpty && !cm_fw.isAbundant && !cm_fw.isShort){

                    if (!cm_fw.strand) cm_fw.dist++;

                    if ((cm_fw.dist != 0) && (cm_fw.dist != cm_fw.size - k_ + 1)){

                        sp.push_back(make_pair(0, cm_fw.dist));
                        sp.push_back(make_pair(cm_fw.dist, cm_fw.size - k_ + 1));

                        splitUnitig(cm_fw.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

                        sp.clear();
                    }
                }
            }

            if (!addUnitig(p.first.toString(), v_kmers_sz)){

                v_kmers[v_kmers_sz].second.setFull();

                v_kmers_sz++;
                added++;
            }
            else h_kmers_ccov.find(p.first.rep())->second.setFull();

            v_joins.push_back(p.first);
        }
        else {

            mapRead(cm);
            it_km += cm.len - 1;
        }
    }

    if (nxt_pos_insert_v_unitigs < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert_v_unitigs);
    if (v_kmers_sz < v_kmers.size()) v_kmers.resize(v_kmers_sz);

    size_t joined = joinAllUnitigs(&v_joins);

    if (verbose){

        cout << "CompactedDBG::add(): Added " << added << " k-mers to the graph" << endl;
        cout << "CompactedDBG::add(): Joined " << joined << " k-mers from the graph" << endl;
    }

    return true;
}

bool CompactedDBG::remove(const UnitigMap& um, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::remove(): Graph is invalid and no unitig can be removed from it" << endl;
        return false;
    }

    vector<Kmer> v_km;

    Kmer head, tail;

    if (um.isShort){

        head = v_kmers[um.pos_unitig].first;
        tail = head;
    }
    else if (um.isAbundant){

        head = h_kmers_ccov.find(um.pos_unitig)->first;
        tail = head;
    }
    else {

        const Unitig* unitig = v_unitigs[um.pos_unitig];

        head = unitig->seq.getKmer(0);
        tail = unitig->seq.getKmer(unitig->numKmers() - 1);
    }

    for (size_t i = 0; i != 4; i++){

        const Kmer bw = head.backwardBase(alpha[i]);

        if (!find(bw, true).isEmpty) v_km.push_back(bw);
    }

    for (size_t i = 0; i != 4; i++){

        const Kmer fw = tail.forwardBase(alpha[i]);

        if (!find(fw, true).isEmpty) v_km.push_back(fw);
    }

    const size_t swap_position = (um.isShort ? v_kmers.size() : v_unitigs.size()) - 1;

    if (!um.isAbundant && (um.pos_unitig != swap_position)) swapUnitigs(um.isShort, um.pos_unitig, swap_position);

    deleteUnitig(um.isShort, um.isAbundant, um.pos_unitig);

    if (um.isShort) v_kmers.resize(swap_position);
    else if (!um.isAbundant) v_unitigs.resize(swap_position);

    joinAllUnitigs(&v_km);

    return true;
}

bool CompactedDBG::filter(const vector<string>& fastx_filename_in, const size_t nb_unique_kmers, const size_t nb_non_unique_kmers,
                            const bool reference_mode, const size_t nb_threads, const size_t nb_bits_unique_kmers_bf,
                            const size_t nb_bits_non_unique_kmers_bf, const size_t read_chunksize, const size_t unitig_size,
                            const bool verbose) {

    if (invalid){

        cerr << "CompactedDBG::filter(): Graph is invalid and it cannot be built" << endl;
        return false;
    }

    BlockedBloomFilter bf_tmp;

    if (reference_mode){

        BlockedBloomFilter tmp(nb_unique_kmers, nb_bits_unique_kmers_bf);
        bf.get(tmp);
        tmp.clear();
    }
    else {

        BlockedBloomFilter tmp(nb_unique_kmers, nb_bits_unique_kmers_bf);
        bf_tmp.get(tmp);
        tmp.clear();

        BlockedBloomFilter tmp2(nb_non_unique_kmers, nb_bits_non_unique_kmers_bf);
        bf.get(tmp2);
        tmp2.clear();
    }

    bool done = false;

    const bool multi_threaded = (nb_threads != 1);

    char name[8192];

    string s;

    size_t name_len, len;

    uint64_t n_read = 0;

    atomic<uint64_t> num_kmers(0), num_ins(0);

    FastqFile FQ(fastx_filename_in);

    vector<string> readv;

    // Main worker thread
    auto worker_function = [&](vector<string>::iterator a, vector<string>::iterator b) {

        uint64_t l_num_kmers = 0, l_num_ins = 0;

        // for each input
        for (auto x = a; x != b; ++x) {

            const char* str = x->c_str();
            int len = x->length();

            KmerHashIterator<RepHash> it_kmer_h(str, len, k_), it_kmer_h_end;
            minHashIterator<RepHash> it_min(str, len, k_, g_, RepHash(), true);

            for (int last_pos = -1; it_kmer_h != it_kmer_h_end; ++it_kmer_h, ++it_min, ++l_num_kmers) {

                std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                // If one or more k-mer were jumped because contained non-ACGT char.
                if (p_.second != last_pos + 1)
                    it_min = minHashIterator<RepHash>(&str[p_.second], len - p_.second, k_, g_, RepHash(), true);

                last_pos = p_.second;
                uint64_t min_hr = it_min.getHash();

                if (!reference_mode) {

                    if (bf_tmp.search_and_insert(p_.first, min_hr, multi_threaded)) ++l_num_ins;
                    //else if (bf.search_and_insert(p_.first, min_hr, multi_threaded)) ++l_num_ins;
                    else bf.insert(p_.first, min_hr);
                }
                else {

                    bf.insert(p_.first, min_hr);
                    ++l_num_ins;
                }
            }
        }

        // atomic adds
        num_kmers += l_num_kmers;
        num_ins += l_num_ins;
    };

    while (!done) {

        readv.clear();
        size_t reads_now = 0;

        while (reads_now < read_chunksize) {

            if (FQ.read_next(name, &name_len, s, &len, NULL, NULL) >= 0) {

                readv.emplace_back(s);
                ++n_read;
                ++reads_now;
            }
            else {

                done = true;
                break;
            }
        }

        vector<thread> workers;

        auto rit = readv.begin();
        size_t batch_size = readv.size() / nb_threads;
        size_t leftover   = readv.size() % nb_threads;

        for (size_t i = 0; i < nb_threads; i++) {

            size_t jump = batch_size + ((i < leftover) ? 1 : 0);

            auto rit_end(rit);
            advance(rit_end, jump);
            workers.push_back(thread(worker_function, rit, rit_end));

            rit = rit_end;
        }

        assert(rit==readv.end());

        for (auto &t : workers) t.join();
    }

    FQ.close();

    if (verbose) {

        cout << "CompactedDBG::filter(): Closed all fasta/fastq files" << endl;
        cout << "CompactedDBG::filter(): Processed " << num_kmers << " k-mers in " << n_read  << " reads" << endl;
        cout << "CompactedDBG::filter(): Found " << num_ins << " unique k-mers" << endl;
        cout << "CompactedDBG::filter(): Number of blocks in Bloom filter is " << bf.getNbBlocks() << endl;
    }

    return true;
}

bool CompactedDBG::construct(const vector<string>& fastx_filename_in, const size_t nb_threads,
                             const size_t read_chunksize, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::construct(): Graph is invalid and cannot be built" << endl;
        return false;
    }

    KmerIterator iter, iterend;

    FastqFile FQ(fastx_filename_in);

    char name[8192];

    string s;

    size_t name_len, len;

    uint64_t n_read = 0;

    vector<string> readv;

    const bool multi_threaded = (nb_threads != 1);

    tiny_vector<Kmer, 3>* fp_candidate = new tiny_vector<Kmer, 3>[bf.getNbBlocks()];

    KmerHashTable<bool> ignored_km_tips;

    auto worker_function = [&](vector<string>::const_iterator a, vector<string>::const_iterator b,
                                vector<NewUnitig>* smallv,
                                vector<tuple<bool, uint64_t, Kmer>>* v_fp_cand,
                                vector<Kmer>* l_ignored_km_tip) {

        uint64_t it_min_h, last_it_min_h;
        std::pair<uint64_t*, uint64_t*> block_bf;

        Kmer km;
        RepHash rep;

        // for each input
        for (auto x = a; x != b; ++x) {

            const char* s_x = x->c_str();

            KmerHashIterator<RepHash> it_kmer_h(s_x, x->length(), k_), it_kmer_h_end;
            preAllocMinHashIterator<RepHash> it_min(s_x, x->length(), k_, g_, rep, true);

            for (int last_pos_km = -2; it_kmer_h != it_kmer_h_end; it_kmer_h++, it_min++) {

                std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                if (p_.second != last_pos_km + 1){ // If one or more k-mer were jumped because contained non-ACGT char.

                    km = Kmer(&s_x[p_.second]);

                    it_min += (last_pos_km == -2 ? p_.second : (p_.second - last_pos_km) - 1);
                    it_min_h = it_min.getHash();

                    block_bf = bf.getBlock(it_min_h);
                }
                else {

                    km.selfForwardBase(s_x[p_.second + k_ - 1]);

                    it_min_h = it_min.getHash();
                    if (it_min_h != last_it_min_h) block_bf = bf.getBlock(it_min_h);
                }

                last_pos_km = p_.second;
                last_it_min_h = it_min_h;

                size_t r = bf.contains_block(p_.first, block_bf);

                if (r != 0){

                    UnitigMap cm = findUnitig(km, *x, p_.second, it_min);

                    if (cm.isEmpty) { // kmer did not map, push into queue for next unitig generation round

                        bool selfLoop = false;
                        bool isIsolated = false;

                        string newseq;

                        size_t pos_match = findUnitigSequenceBBF(km, newseq, selfLoop, isIsolated, *l_ignored_km_tip); //Build unitig from Bloom filter

                        if (isIsolated){ // According to the BF, k-mer is isolated in the graph and is a potential false positive

                            const uint64_t block = ((r == 1 ? block_bf.first : block_bf.second) - bf.getTable_ptr()) / NB_ELEM_BLOCK;

                            Kmer km_rep = km.rep();
                            const tiny_vector<Kmer, 3>& v = fp_candidate[block];
                            tuple<bool, uint64_t, Kmer> t_fp_cand = make_tuple(true, block, km_rep);

                            for (auto km_tmp : v){ // Go over global list of existing FP candidates

                                if (km_tmp == km_rep){ // If already stored as a FP candidate, it is a TP

                                    std::get<0>(t_fp_cand) = false; // K-mer must be removed from list of FP candidate
                                    break;
                                }
                            }

                            if (!std::get<0>(t_fp_cand)){ // If k-mer was a FP candidate now turned into a TP, insert into data structure

                                smallv->emplace_back(km, *x, p_.second, selfLoop ? std::string() : newseq);
                                v_fp_cand->emplace_back(t_fp_cand); // TP will be removed from list of FP later
                            }
                            else { // Need to make sure the FP candidate was not already inserted locally as a FP candidate (else it is a TP)

                                bool found_fp = false;

                                for (auto& t_fp_cand_tmp : *v_fp_cand){ // Go over local (to the thread) list of existing FP candidates

                                    if (std::get<2>(t_fp_cand_tmp) == km_rep){ // If present as FP candidate

                                        std::get<0>(t_fp_cand_tmp) = false; // Indicate FP must be removed from list of false positive candidate
                                        found_fp = true;
                                        break;
                                    }
                                }

                                if (found_fp) smallv->emplace_back(km, *x, p_.second, selfLoop ? std::string() : newseq);
                                else v_fp_cand->emplace_back(t_fp_cand);
                            }
                        }
                        else {

                            smallv->emplace_back(km, *x, p_.second, selfLoop ? std::string() : newseq);
                            it_kmer_h += cstrMatch(&s_x[p_.second + k_], &newseq.c_str()[pos_match + k_]);
                        }
                    }
                    else {

                        mapRead(cm);
                        it_kmer_h += cm.len - 1;
                    }
                }
            }
        }
    };

    vector<vector<NewUnitig>> parray(nb_threads);
    vector<vector<tuple<bool, uint64_t, Kmer>>> v_fp_cand_threads(nb_threads);
    vector<vector<Kmer>> v_ignored_km_tip_thread(nb_threads);

    int round = 0;
    bool done = false;

    while (!done) {

        readv.clear();

        size_t reads_now = 0;

        while (reads_now < read_chunksize) {

            if (FQ.read_next(name, &name_len, s, &len, NULL, NULL) >= 0){

                readv.emplace_back(s);
                ++reads_now;
                ++n_read;
            }
            else {

                done = true;
                break;
            }
        }

        ++round;

        if (verbose) cout << "CompactedDBG::construct(): Starting round " << round << endl;

        // run parallel code
        vector<thread> workers;
        auto rit = readv.begin();
        size_t batch_size = readv.size() / nb_threads;
        size_t leftover   = readv.size() % nb_threads;

        for (size_t i = 0; i < nb_threads; i++) {

            size_t jump = batch_size + (i < leftover ? 1 : 0);
            auto rit_end(rit);

            advance(rit_end, jump);
            workers.push_back(thread(worker_function, rit, rit_end, &parray[i],
                                     &v_fp_cand_threads[i], &v_ignored_km_tip_thread[i]));
            rit = rit_end;
        }

        assert(rit == readv.end());

        for (auto& t : workers) t.join();

        for (auto &v : parray) { // for each thread

            for (auto &x : v) addUnitigSequenceBBF(x.km, x.read, x.pos, x.seq, v_ignored_km_tip_thread[0]); // add each unitig for this thread

            v.clear(); //clear the map
        }

        for (auto &v : v_fp_cand_threads) { // for each thread

            for (auto &x : v){ //for each FP cadndiate to add or delete

                Kmer km_fp = std::get<2>(x);

                tiny_vector<Kmer, 3>& v = fp_candidate[std::get<1>(x)];

                size_t i = 0;

                for (; i < v.size(); i++){ //Search vector for k-mer

                    if (v[i] == km_fp) break;
                }

                if (std::get<0>(x) == false){ // If FP is to delete

                    if (i < v.size()) v.remove(i); // If k-mer has not already been deleted

                    UnitigMap cm = find(km_fp); // Find FP to delete (was inserted into data structure as TP)

                    if (cm.isEmpty){

                        empty();
                        delete[] fp_candidate;

                        cerr << "CompactedDBG::construct(): Deleted false positive candidate was not inserted into data structure" << endl;
                        return false;
                    }

                    mapRead(cm); //Map the TP a second time;
                }
                else if (i >= v.size()){ //

                    if (!multi_threaded || find(km_fp).isEmpty) v.push_back(km_fp); //Add k-mer to vector if not already present
                }
                else { // K-mer already added as a FP candidate by other thread: it is now a TP

                    v.remove(i); //// K-mer deleted from list of FP candidates

                    const string str_km = km_fp.toString();

                    addUnitigSequenceBBF(km_fp, str_km, 0, str_km, v_ignored_km_tip_thread[0]);

                    UnitigMap cm = find(km_fp); // Find FP to delete (was inserted into data structure as TP)

                    if (cm.isEmpty){

                        empty();
                        delete[] fp_candidate;

                        cerr << "CompactedDBG::construct(): Deleted false positive candidate was not inserted into data structure" << endl;
                        return false;
                    }

                    mapRead(cm); //Map the TP a second time;
                }
            }

            v.clear(); //clear the map
        }

        for (auto &v : v_ignored_km_tip_thread) { // for each thread

            for (auto x : v) ignored_km_tips.insert(make_pair(x, false));

            v.clear();
        }

        if (read_chunksize > 1 && verbose) {

            cout << "CompactedDBG::construct(): End of round" << endl;
            cout << "CompactedDBG::construct(): Processed " << size() << " unitigs" << endl;
        }
    }

    FQ.close();

    parray.clear();
    v_fp_cand_threads.clear();
    v_ignored_km_tip_thread.clear();

    delete[] fp_candidate;

    if (verbose) cout << "CompactedDBG::construct(): Closed all fasta/fastq files" << endl;

    size_t unitigsBefore = size();

    if (verbose) cout << endl << "CompactedDBG::construct(): Splitting unitigs (1/2)" << endl;

    pair<size_t, size_t> unitigSplit = splitAllUnitigs();

    int unitigsAfter1 = size();

    if (verbose) cout << endl << "CompactedDBG::construct(): Splitting unitigs (2/2)" << endl;

    check_fp_tips(ignored_km_tips);
    ignored_km_tips.clear();

    int unitigsAfter2 = size();

    if (verbose) {

        cout << "CompactedDBG::construct(): Before split: " << unitigsBefore << " unitigs" << endl;
        cout << "CompactedDBG::construct(): After split (1/2): " << unitigsAfter1 << " unitigs" <<  endl;
        cout << "CompactedDBG::construct(): After split (2/2): " << unitigsAfter2 << " unitigs" <<  endl;
        cout << "CompactedDBG::construct(): Unitigs split: " << unitigSplit.first << endl;
        cout << "CompactedDBG::construct(): Unitigs deleted: " << unitigSplit.second << endl;
    }

    if (verbose) cout << endl << "CompactedDBG::construct(): Joining unitigs" << endl;

    size_t joined = joinAllUnitigs();

    int unitigsAfter3 = size();

    if (verbose) {

        cout << "CompactedDBG::construct(): After join: " << unitigsAfter3 << " unitigs" << endl;
        cout << "CompactedDBG::construct(): Joined " << joined << " unitigs" << endl;
    }

    return true;
}

// use:  cm.mapRead(km,pos,cc)
// pre:  cc is a reference to a current unitig in cm, km maps to cc
// post: the coverage information in cc has been updated
void CompactedDBG::mapRead(const UnitigMap& cc) {

    if (cc.isEmpty) return; // nothing maps, move on

    if (cc.isShort) v_kmers[cc.pos_unitig].second.cover(cc.dist, cc.dist + cc.len - 1);
    else if (cc.isAbundant) h_kmers_ccov.find(cc.pos_unitig)->second.cover(cc.dist, cc.dist + cc.len - 1);
    else v_unitigs[cc.pos_unitig]->ccov.cover(cc.dist, cc.dist + cc.len - 1);
}

// use: b = cm.addUnitig(km,read)
// pre:
// post: either unitig string containsin has been added and b == true
//       or it was present and the coverage information was updated, b == false
//       NOT Threadsafe!
//bool UnitigMapper::addUnitigSequenceBBF(Kmer km, const string& read, size_t pos, const string& seq) {
bool CompactedDBG::addUnitigSequenceBBF(Kmer km, const string& read, size_t pos, const string& seq, vector<Kmer>& l_ignored_km_tip) {

    string s;
    bool selfLoop = false;
    bool isIsolated = false;

    if (!seq.empty()) s = seq;
    else findUnitigSequenceBBF(km, s, selfLoop, isIsolated, l_ignored_km_tip);

    size_t k = Kmer::k;

    if (selfLoop) {

        bool foundAny = false;

        for (KmerIterator it(s.c_str()), it_end; it != it_end; ++it) {

            UnitigMap cm = find(it->first);

            if (!cm.isEmpty) {

                mapRead(cm);
                foundAny = true;
            }
        }

        if (!foundAny) {

            addUnitig(s, s.length() == k ? v_kmers.size() : v_unitigs.size());

            for (KmerIterator it(s.c_str()), it_end; it != it_end; ++it) mapRead(find(it->first));
        }

        return true;
    }

    UnitigMap cm = findUnitig(km, read, pos);

    if (cm.isEmpty){

        addUnitig(s, s.length() == k ? v_kmers.size() : v_unitigs.size());
        cm = findUnitig(km, read, pos);
    }

    mapRead(cm);

    return !cm.isEmpty;
}

// use:  cm.findUnitigSequenceBBF(km, s, selfLoop)
// pre:  km is in the bloom filter
// post: s is the unitig containing the kmer km
//       and the first k-mer in s is smaller (wrt. < operator)
//       than the last kmer
//       selfLoop is true of the unitig is a loop or hairpin
size_t CompactedDBG::findUnitigSequenceBBF(Kmer km, string& s, bool& selfLoop, bool& isIsolated, vector<Kmer>& l_ignored_km_tip) {

    string fw_s;

    Kmer end = km;
    Kmer last = end;
    Kmer twin = km.twin();

    char c;

    size_t j = 0;

    bool has_no_neighbor;

    selfLoop = false;
    isIsolated = false;

    while (fwStepBBF(end, end, c, has_no_neighbor, l_ignored_km_tip)) {

        j++;

        if (end == km) {
            selfLoop = true;
            break;
        }
        else if (end == twin) break;
        else if (end == last.twin()) break;

        fw_s.push_back(c);
        last = end;
    }

    string bw_s;
    Kmer front = km;
    Kmer first = front;

    if (!selfLoop) {

        isIsolated = (j == 0) && has_no_neighbor;
        j = 0;

        while (bwStepBBF(front, front, c, has_no_neighbor, l_ignored_km_tip)) {

            j++;

            if (front == km) {
                selfLoop = true;
                break;
            }
            else if (front == twin) break;
            else if (front == first.twin()) break;

            bw_s.push_back(c);
            first = front;
        }

        if (isIsolated) isIsolated = (j == 0) && has_no_neighbor;

        reverse(bw_s.begin(), bw_s.end());
    }

    s.reserve(Kmer::k + fw_s.size() + bw_s.size());
    s.append(bw_s);

    char tmp[Kmer::MAX_K];
    km.toString(tmp);
    s.append(tmp);

    s.append(fw_s);

    return bw_s.size();
}

// use:  cc = cm.findUnitig(km,s,pos)
// pre:  s[pos,pos+k-1] is the kmer km
// post: cc contains either the reference to the unitig position
//       or empty if none found
UnitigMap CompactedDBG::findUnitig(const Kmer& km, const string& s, size_t pos) const {

    // need to check if we find it right away, need to treat this common case
    UnitigMap cc = find(km);

    if (!cc.isEmpty && !cc.isShort && !cc.isAbundant){

        const CompressedSequence& seq = v_unitigs[cc.pos_unitig]->seq;
        size_t km_dist = cc.dist;
        size_t jlen = 0;
        size_t k = Kmer::k;

        if (cc.strand) jlen = seq.jump(s.c_str(), pos, cc.dist, false) - k + 1;
        else {

            jlen = seq.jump(s.c_str(), pos, cc.dist + k - 1, true) - k + 1; // match s_fw to comp(seq)_bw
            km_dist -= jlen - 1;
        }

        return UnitigMap(cc.pos_unitig, km_dist, jlen, cc.size, false, false, cc.strand, cc.dbg);
    }

    return cc;
}

UnitigMap CompactedDBG::findUnitig(const Kmer& km, const string& s, size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h) const {

    // need to check if we find it right away, need to treat this common case
    UnitigMap cc = find(km, it_min_h);

    if (!cc.isEmpty && !cc.isShort && !cc.isAbundant){

        const CompressedSequence& seq = v_unitigs[cc.pos_unitig]->seq;
        size_t km_dist = cc.dist;
        size_t jlen = 0;
        size_t k = Kmer::k;

        if (cc.strand) jlen = seq.jump(s.c_str(), pos, cc.dist, false) - k + 1;
        else {

            jlen = seq.jump(s.c_str(), pos, cc.dist + k - 1, true) - k + 1; // match s_fw to comp(seq)_bw
            km_dist -= jlen - 1;
        }

        return UnitigMap(cc.pos_unitig, km_dist, jlen, cc.size, false, false, cc.strand, cc.dbg);
    }

    return cc;
}

UnitigMap CompactedDBG::find(const Kmer& km, bool extremities_only) const {

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int k = Kmer::k;
    const int diff = k - Minimizer::g;

    char km_tmp[k + 1];
    km.toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k, k, Minimizer::g, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz = Minimizer(&km_tmp[min_h_res.pos]).rep();
        hmap_min_unitigs_t::const_iterator it = hmap_min_unitigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const tiny_vector<size_t,tiny_vector_sz>& v = it->second;
            const int v_sz = v.size();

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; i++){

                unitig_id_pos = v[i];
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()) return UnitigMap(it_km.getHash(), 0, 1, Kmer::k, false, true, km == km_rep, this);
                    }

                    if ((unitig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(&km_tmp[mhr.pos]).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (unitig_id_pos & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos &= MASK_CONTIG_POS;

                    if (isShort){

                        if (min_h_res.pos == unitig_id_pos){

                            if (v_kmers[unitig_id].first == km_rep) return UnitigMap(unitig_id, 0, 1, k, true, false, true, this);
                        }
                        else if ((min_h_res.pos == diff - unitig_id_pos) && (v_kmers[unitig_id].first == km_rep)){

                            return UnitigMap(unitig_id, 0, 1, k, true, false, false, this);
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k;
                        pos_match = unitig_id_pos - min_h_res.pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, km)){

                                return UnitigMap(unitig_id, pos_match, 1, len + k, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, km_twin)){

                                return UnitigMap(unitig_id, pos_match, 1, len + k, false, false, false, this);
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, km)){

                                return UnitigMap(unitig_id, pos_match, 1, len + k, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, km_twin)){

                                return UnitigMap(unitig_id, pos_match, 1, len + k, false, false, false, this);
                            }
                        }
                    }
                }
            }
        }

        it_it_min++;
    }

    return UnitigMap();
}

UnitigMap CompactedDBG::find(const Kmer& km, const preAllocMinHashIterator<RepHash>& it_min_h) const {

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort;

    size_t unitig_id, len;

    int64_t pos_match;

    const int k = Kmer::k;
    const int diff = k - Minimizer::g;

    preAllocMinHashIterator<RepHash> it_min = preAllocMinHashIterator<RepHash>(it_min_h, k);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz = Minimizer(&it_min.s[min_h_res.pos]).rep();
        hmap_min_unitigs_t::const_iterator it = hmap_min_unitigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const tiny_vector<size_t,tiny_vector_sz>& v = it->second;

            it = hmap_min_unitigs.end();

            for (auto unitig_id_pos : v){ // For each unitig associated with the minimizer

                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()){

                            return UnitigMap(it_km.getHash(), 0, 1, Kmer::k, false, true, km == km_rep, this);
                        }
                    }

                    if ((unitig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(&it_min.s[mhr.pos]).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (unitig_id_pos & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos &= MASK_CONTIG_POS;

                    if (isShort){

                        if (min_h_res.pos == unitig_id_pos){

                            if (v_kmers[unitig_id].first == km_rep){

                                return UnitigMap(unitig_id, 0, 1, k, true, false, true, this);
                            }
                        }
                        else if ((min_h_res.pos == diff - unitig_id_pos) && (v_kmers[unitig_id].first == km_rep)){

                            return UnitigMap(unitig_id, 0, 1, k, true, false, false, this);
                        }
                    }
                    else {

                        pos_match = unitig_id_pos - min_h_res.pos;
                        len = v_unitigs[unitig_id]->length() - k;

                        if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, km)){

                            return UnitigMap(unitig_id, pos_match, 1, len + k, false, false, true, this);
                        }

                        pos_match = unitig_id_pos - diff + min_h_res.pos;

                        if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, km_twin)){

                            return UnitigMap(unitig_id, pos_match, 1, len + k, false, false, false, this);
                        }
                    }
                }
            }
        }

        it_it_min++;
    }

    return UnitigMap();
}

bool CompactedDBG::bwStepBBF(Kmer km, Kmer& front, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand) const {

    size_t i, j = -1, j_tmp, k = Kmer::k, g = Minimizer::g;
    size_t nb_neigh = 0;

    char km_tmp[k + 1];

    front.backwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
    std::pair<uint64_t*, uint64_t*> block = bf.getBlock(it_min_h);

    RepHash rep_h(k - 1), rep_h_cpy;
    rep_h.init(km_tmp + 1);

    int found_fp_bw = 0;
    Kmer km_fp;

    bool pres_neigh_bw[4] = {false, false, false, false};

    for (i = 0; i < 4; ++i) {

        rep_h_cpy = rep_h;
        rep_h_cpy.extendBW(alpha[i]);

        if (bf.contains(rep_h_cpy.hash(), block)){

            j = i;
            pres_neigh_bw[i] = true;
            nb_neigh++;

            if (!check_fp_cand && (nb_neigh >= 2)) break;
        }
    }

    if (check_fp_cand && (nb_neigh >= 2)){

        for (i = 0; i < 4; ++i) {

            if (pres_neigh_bw[i]){

                char dummy;
                bool has_no_neighbor_tmp = false;

                km_tmp[0] = alpha[i];

                km_fp = Kmer(km_tmp);

                bwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                if (has_no_neighbor_tmp && fwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false)) found_fp_bw++;
                else {

                    j_tmp = i;
                    pres_neigh_bw[i] = false;
                }
            }
        }

        if (found_fp_bw != 0){

            if ((nb_neigh - found_fp_bw) != 0){

                j = j_tmp;
                nb_neigh -= found_fp_bw;
            }
            else found_fp_bw = 0;
        }
    }

    if (nb_neigh != 1) {

        has_no_neighbor = (nb_neigh == 0);
        return false;
    }
    else has_no_neighbor = false;

    if (check_fp_cand){

        nb_neigh = 0;

        int found_fp_fw = 0;
        bool pres_neigh_fw[4] = {false, false, false, false};

        Kmer bw = front.backwardBase(alpha[j]);

        bw.forwardBase('A').toString(km_tmp);

        it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
        block = bf.getBlock(it_min_h);

        rep_h.init(km_tmp);

        for (i = 0; i < 4; ++i) {

            rep_h_cpy = rep_h;
            rep_h_cpy.extendFW(alpha[i]);

            if (bf.contains(rep_h_cpy.hash(), block)){

                nb_neigh++;
                pres_neigh_fw[i] = true;
            }
        }

        if (nb_neigh >= 2){

            for (i = 0; i < 4; ++i) {

                if (pres_neigh_fw[i]){

                    char dummy;
                    bool add = true, has_no_neighbor_tmp = false;

                    km_tmp[k - 1] = alpha[i];

                    km_fp = Kmer(km_tmp);

                    fwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                    if (has_no_neighbor_tmp && bwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false)){

                        if (km_fp != km) found_fp_fw++;
                        else {

                            found_fp_fw = 0;
                            break;
                        }
                    }
                    else pres_neigh_fw[i] = false;
                }
            }

            if (found_fp_fw != 0){

                if ((nb_neigh - found_fp_fw) != 0) nb_neigh -= found_fp_fw;
                else found_fp_fw = 0;
            }
        }

        if (nb_neigh != 1) return false;

        if (bw != km) {
            // exactly one k-mer in bw, character used is c

            for (i = 0; (i < 4) && (found_fp_fw != 0); ++i) {

                if (pres_neigh_fw[i]){

                    km_tmp[k - 1] = alpha[i];
                    km_fp = Kmer(km_tmp).rep();

                    l_ignored_km_tip.push_back(km_fp);

                    found_fp_fw--;
                }
            }

            front.backwardBase('A').toString(km_tmp);

            for (i = 0; (i < 4) && (found_fp_bw != 0); ++i) {

                if (pres_neigh_bw[i]){

                    km_tmp[0] = alpha[i];
                    km_fp = Kmer(km_tmp).rep();

                    l_ignored_km_tip.push_back(km_fp);

                    found_fp_bw--;
                }
            }

            front = bw;
            c = alpha[j];

            return true;
        }

        return false;
    }

    return true;
}

bool CompactedDBG::fwStepBBF(Kmer km, Kmer& end, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand) const {

    size_t i, j = -1, j_tmp, k = Kmer::k, g = Minimizer::g;
    size_t nb_neigh = 0;

    char km_tmp[k + 1];

    end.forwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
    std::pair<uint64_t*, uint64_t*> block = bf.getBlock(it_min_h);

    RepHash rep_h(k - 1), rep_h_cpy;
    rep_h.init(km_tmp);

    int found_fp_fw = 0;
    Kmer km_fp;

    bool pres_neigh_fw[4] = {false, false, false, false};

    for (i = 0; i < 4; ++i) {

        //km_tmp[k - 1] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendFW(alpha[i]);

        if (bf.contains(rep_h_cpy.hash(), block)){

            j = i;
            pres_neigh_fw[i] = true;
            nb_neigh++;

            if (!check_fp_cand && (nb_neigh >= 2)) break;
        }
    }

    if (check_fp_cand && (nb_neigh >= 2)){

        for (i = 0; i < 4; ++i) {

            if (pres_neigh_fw[i]){

                km_tmp[k - 1] = alpha[i];

                char dummy;
                km_fp = Kmer(km_tmp);

                bool has_no_neighbor_tmp = false;

                fwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                if (has_no_neighbor_tmp && bwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false)) found_fp_fw++;
                else {

                    j_tmp = i;
                    pres_neigh_fw[i] = false;
                }
            }
        }

        if (found_fp_fw != 0){

            if ((nb_neigh - found_fp_fw) != 0){

                j = j_tmp;
                nb_neigh -= found_fp_fw;
            }
            else found_fp_fw = 0;
        }
    }

    if (nb_neigh != 1) {

        has_no_neighbor = (nb_neigh == 0);
        return false;
    }
    else has_no_neighbor = false;

    if (check_fp_cand){

        nb_neigh = 0;

        int found_fp_bw = 0;
        bool pres_neigh_bw[4] = {false, false, false, false};

        Kmer fw = end.forwardBase(alpha[j]);

        // check bw from fw link
        fw.backwardBase('A').toString(km_tmp);

        it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
        block = bf.getBlock(it_min_h);

        rep_h.init(km_tmp + 1);

        for (i = 0; i < 4; ++i) {

            //km_tmp[0] = alpha[i];

            rep_h_cpy = rep_h;
            rep_h_cpy.extendBW(alpha[i]);

            if (bf.contains(rep_h_cpy.hash(), block)){

                nb_neigh++;
                pres_neigh_bw[i] = true;
            }
        }

        if (nb_neigh >= 2){

            for (i = 0; i < 4; ++i) {

                if (pres_neigh_bw[i]){

                    char dummy;
                    bool has_no_neighbor_tmp = false;

                    km_tmp[0] = alpha[i];

                    km_fp = Kmer(km_tmp);

                    bwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                    if (has_no_neighbor_tmp && fwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false)){

                        if (km_fp != km) found_fp_bw++;
                        else {

                            found_fp_bw = 0;
                            break;
                        }
                    }
                    else pres_neigh_bw[i] = false;
                }
            }

            if (found_fp_bw != 0){

                if ((nb_neigh - found_fp_bw) != 0) nb_neigh -= found_fp_bw;
                else found_fp_bw = 0;
            }
        }

        if (nb_neigh != 1) return false;

        if (fw != km) {

            for (i = 0; (i < 4) && (found_fp_bw != 0); ++i) {

                if (pres_neigh_bw[i]){

                    km_tmp[0] = alpha[i];
                    km_fp = Kmer(km_tmp).rep();

                    l_ignored_km_tip.push_back(km_fp);

                    found_fp_bw--;
                }
            }

            end.forwardBase('A').toString(km_tmp);

            for (i = 0; (i < 4) && (found_fp_fw != 0); ++i) {

                if (pres_neigh_fw[i]){

                    km_tmp[k - 1] = alpha[i];
                    km_fp = Kmer(km_tmp).rep();

                    l_ignored_km_tip.push_back(km_fp);

                    found_fp_fw--;
                }
            }

            end = fw;
            c = alpha[j];
            return true;
        }

        return false;
    }

    return true;
}

bool CompactedDBG::addUnitig(const string& str_unitig, const size_t id_unitig){

    int k = Kmer::k, g = Minimizer::g, pos;

    size_t len = str_unitig.size();
    size_t pos_id_unitig = id_unitig << 32;

    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    bool isShort = false, isAbundant = false, isForbidden = false;

    char* c_str = const_cast<char*>(str_unitig.c_str());

    char km_tmp[k + 1];

    Kmer km_rep;

    if (len == k){ // Unitig to add is short, maybe abundant as well

        isShort = true;

        pos_id_unitig |= MASK_CONTIG_TYPE;

        km_rep = Kmer(c_str).rep();
        km_rep.toString(km_tmp);
        c_str = km_tmp;
    }

    minHashIterator<RepHash> it_min(c_str, len, k, Minimizer::g, RepHash(), true), it_min_end;

    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){

        //If current minimizer was not seen before
        if ((last_pos_min < it_min.getPosition()) || isForbidden){

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&c_str[min_h_res.pos]).rep(); //Get the minimizer to insert
                std::pair<hmap_min_unitigs_t::iterator, bool> p = hmap_min_unitigs.insert(make_pair(minz_rep, tiny_vector<size_t,tiny_vector_sz>()));
                size_t v_sz = p.first->second.size();

                pos = min_h_res.pos;
                pos_id_unitig = (pos_id_unitig & mask) | ((size_t) pos);

                if (!isShort){

                    mhr = min_h_res;

                    while ((v_sz >= max_abundance_lim) || ((v_sz > 0) && ((p.first->second[v_sz-1] & mask) == mask))){

                        mhr_tmp = it_min.getNewMin(mhr);
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            if ((p.first->second[v_sz-1] & mask) != mask){ // Minimizer was never signaled before as overcrowded

                                // If minimizer bin already contains abundant k-mer, just set flag for unitig overcrowding
                                if ((p.first->second[v_sz-1] & MASK_CONTIG_ID) == MASK_CONTIG_ID) p.first->second[v_sz-1] |= MASK_CONTIG_TYPE;
                                else p.first->second.push_back(mask);
                            }

                            mhr = mhr_tmp;
                            minz_rep = Minimizer(&c_str[mhr.pos]).rep();
                            p = hmap_min_unitigs.insert(make_pair(minz_rep, tiny_vector<size_t,tiny_vector_sz>()));
                            v_sz = p.first->second.size();
                        }
                        else break;
                    }
                }

                tiny_vector<size_t,tiny_vector_sz>& v = p.first->second;

                if (v_sz == 0) v.push_back(pos_id_unitig); //Newly created vector, just push unitig ID
                else if (isShort && (v_sz >= min_abundance_lim)){ //The minimizer (is/might be) too abundant

                    isShort = false;
                    isAbundant = true;

                    it_min = it_min_end;

                    break;
                }
                else if ((v[v_sz-1] & MASK_CONTIG_ID) == MASK_CONTIG_ID){ //If minimizer is abundant or crowded

                    if ((v_sz == 1) || (v[v_sz-2] != pos_id_unitig)) v.insert(pos_id_unitig, v_sz-1);
                }
                else if (v[v_sz-1] != pos_id_unitig) v.push_back(pos_id_unitig);

                last_pos_min = min_h_res.pos;
                it_it_min++;
            }
        }
    }

    if (isAbundant){

        if (id_unitig == v_kmers.size()) v_kmers.push_back(make_pair(km_rep, CompressedCoverage(1)));
        else v_kmers[id_unitig] = make_pair(km_rep, CompressedCoverage(1));

        deleteUnitig(true, false, id_unitig);
        if (id_unitig == v_kmers.size() - 1) v_kmers.resize(v_kmers.size() - 1);

        it_min = minHashIterator<RepHash>(c_str, len, k, Minimizer::g, RepHash(), true);

        for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){

            if (last_pos_min < it_min.getPosition()){ //If current minimizer was not seen before

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&c_str[min_h_res.pos]).rep(); //Get the minimizer to insert

                    std::pair<hmap_min_unitigs_t::iterator, bool> p = hmap_min_unitigs.insert(make_pair(minz_rep, tiny_vector<size_t,tiny_vector_sz>()));

                    tiny_vector<size_t,tiny_vector_sz>& v = p.first->second;
                    size_t v_sz = v.size();

                    if ((v_sz > 0) && ((v[v_sz-1] & MASK_CONTIG_ID) == MASK_CONTIG_ID)) v[v_sz-1]++;
                    else v.push_back(MASK_CONTIG_ID + 1);

                    last_pos_min = min_h_res.pos;
                    it_it_min++;
                }
            }
        }

        h_kmers_ccov.insert(make_pair(km_rep, CompressedCoverage(1)));
    }
    else if (isShort){

        if (id_unitig == v_kmers.size()) v_kmers.push_back(make_pair(km_rep, CompressedCoverage(1)));
        else v_kmers[id_unitig] = make_pair(km_rep, CompressedCoverage(1));
    }
    else if (id_unitig == v_unitigs.size()) v_unitigs.push_back(new Unitig(c_str)); //Push unitig to list of unitigs
    else v_unitigs[id_unitig] = new Unitig(c_str);

    return isAbundant;
}

void CompactedDBG::deleteUnitig(const bool isShort, const bool isAbundant, const size_t id_unitig){

    int k = Kmer::k;

    if (isAbundant){

        char km_str[k + 1];

        Kmer km = h_kmers_ccov.find(id_unitig)->first;

        km.toString(km_str);

        minHashIterator<RepHash> it_min(km_str, k, k, Minimizer::g, RepHash(), true), it_min_end;

        for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ // Iterate over minimizers of unitig to delete

            if (last_pos_min < it_min.getPosition()){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&km_str[min_h_res.pos]).rep(); // Get canonical minimizer

                    hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(minz_rep); // Look for the minimizer in the hash table

                    if (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        tiny_vector<size_t, tiny_vector_sz>& v = it_h->second;
                        size_t last_pos_v = v.size() - 1;

                        v[last_pos_v]--;

                        if (((v[last_pos_v] & RESERVED_ID) == 0) && ((v[last_pos_v] & MASK_CONTIG_TYPE) != MASK_CONTIG_TYPE)){

                            if (last_pos_v == 0) hmap_min_unitigs.erase(minz_rep);
                            else v.remove(v.size() - 1);
                        }
                    }

                    last_pos_min = min_h_res.pos;
                    it_it_min++;
                }
            }
        }

        h_kmers_ccov.erase(km);

        return;
    }

    bool isForbidden = false;

    size_t pos_id_unitig = id_unitig << 32;
    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    string str;

    if (isShort){

        str = v_kmers[id_unitig].first.toString();
        pos_id_unitig |= MASK_CONTIG_TYPE;
    }
    else str = v_unitigs[id_unitig]->seq.toString();

    const char* s = str.c_str();

    size_t len = str.size();

    minHashIterator<RepHash> it_min(s, len, k, Minimizer::g, RepHash(), true), it_min_end;

    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ // Iterate over minimizers of unitig to delete

        if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer hash is found in unitig to delete

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep(); // Get canonical minimizer
                hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(minz_rep); // Look for the minimizer in the hash table

                mhr = min_h_res;

                while (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                    tiny_vector<size_t, tiny_vector_sz>& v = it_h->second;
                    const int v_sz = v.size();

                    for (size_t i = 0; i < v_sz; i++){

                        if ((v[i] & mask) == pos_id_unitig){

                            v.remove(i);
                            break;
                        }
                    }

                    it_h = hmap_min_unitigs.end();

                    if (v.size() == 0) hmap_min_unitigs.erase(minz_rep);
                    else if (!isShort && ((v[v_sz-1] & mask) == mask)){ //Minimizer bin is overcrowded

                        mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz_rep = Minimizer(&s[mhr.pos]).rep();
                            it_h = hmap_min_unitigs.find(minz_rep);
                        }
                        else break;
                    }
                }

                last_pos_min = min_h_res.pos;
                it_it_min++;
            }
        }
    }

    // The unitig is deleted but its space in the unitig vector is not because:
    // 1 - It would change indices in the minimizer hash table
    if (isShort) v_kmers[id_unitig].first.set_deleted();
    else {

        delete v_unitigs[id_unitig];
        v_unitigs[id_unitig] = NULL;
    }
}

void CompactedDBG::swapUnitigs(const bool isShort, const size_t id_a, const size_t id_b){

    size_t shift_id_unitig_a = id_a << 32;
    size_t shift_id_unitig_b = id_b << 32;

    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    string str;

    // Swap the unitig pointers in v_unitigs
    if (isShort){

        std::swap(v_kmers[id_a], v_kmers[id_b]);

        shift_id_unitig_a |= MASK_CONTIG_TYPE;
        shift_id_unitig_b |= MASK_CONTIG_TYPE;

        str = v_kmers[id_a].first.toString();
    }
    else {

        std::swap(v_unitigs[id_a], v_unitigs[id_b]);

        str = v_unitigs[id_a]->seq.toString();
    }

    bool isForbidden = false;

    // Swap the unitig IDs in the minimizer hash table
    const char* s = str.c_str();

    size_t len = str.size();

    vector<Minimizer> v_min_a;

    minHashIterator<RepHash> it_min(s, len, Kmer::k, Minimizer::g, RepHash(), true), it_min_end;

    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ // Iterate over minimizers of unitig

        if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer is found in unitig

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep();

                if (!isShort){

                    hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(minz_rep);

                    if (it_h != hmap_min_unitigs.end()){

                        v_min_a.push_back(minz_rep); //Add minimizer to list of minimizers

                        size_t v_sz = it_h->second.size();

                        mhr = min_h_res;

                        while ((it_h->second[v_sz-1] & mask) == mask){

                            mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                            isForbidden = true;

                            if (mhr_tmp.hash != mhr.hash){

                                minz_rep = Minimizer(&s[mhr_tmp.pos]).rep();
                                it_h = hmap_min_unitigs.find(minz_rep);

                                if (it_h == hmap_min_unitigs.end()) break;

                                mhr = mhr_tmp;
                                v_sz = it_h->second.size();

                                v_min_a.push_back(minz_rep);
                            }
                            else break;
                        }
                    }
                }
                else v_min_a.push_back(minz_rep); //Add minimizer to list of minimizers

                last_pos_min = min_h_res.pos;
                it_it_min++;
            }
        }
    }

    sort(v_min_a.begin(), v_min_a.end()); //Sort the list of minimizers

    for (vector<Minimizer>::iterator it_v_min = v_min_a.begin(); it_v_min != v_min_a.end(); it_v_min++){ // Iterate over minimizers

        if ((it_v_min == v_min_a.begin()) || (*it_v_min != *(it_v_min-1))){ //If the minimizer is diff. from the previous one

            hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(*it_v_min); // Look for the minimizer in the hash table

            if (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                tiny_vector<size_t,tiny_vector_sz>& v_id_unitigs = it_h->second;
                tiny_vector<size_t,tiny_vector_sz>::iterator it_v_c = v_id_unitigs.begin();

                // Iterate over unitig IDs associated with this minimizer
                for (; it_v_c != v_id_unitigs.end(); it_v_c++){

                     // Swap the unitig ids but do not change positions;
                    if ((*it_v_c & mask) == shift_id_unitig_b) *it_v_c = shift_id_unitig_a | (*it_v_c & MASK_CONTIG_POS);
                    else if ((*it_v_c & mask) == shift_id_unitig_a) *it_v_c = shift_id_unitig_b | (*it_v_c & MASK_CONTIG_POS);
                }
            }
        }
    }

    vector<Minimizer> v_min_b;

    str = isShort ? v_kmers[id_b].first.toString() : v_unitigs[id_b]->seq.toString();
    s = str.c_str();
    len = str.size();

    isForbidden = false;

    it_min = minHashIterator<RepHash>(s, len, Kmer::k, Minimizer::g, RepHash(), true);

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ // Iterate over minimizers of unitig

        if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer is found in unitig

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep();

                if (!isShort){

                    hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(minz_rep);

                    if (it_h != hmap_min_unitigs.end()){

                        v_min_b.push_back(minz_rep); //Add minimizer to list of minimizers

                        size_t v_sz = it_h->second.size();

                        mhr = min_h_res;

                        while ((it_h->second[v_sz-1] & mask) == mask){

                            mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                            isForbidden = true;

                            if (mhr_tmp.hash != mhr.hash){

                                minz_rep = Minimizer(&s[mhr_tmp.pos]).rep();
                                it_h = hmap_min_unitigs.find(minz_rep);

                                if (it_h == hmap_min_unitigs.end()) break;

                                mhr = mhr_tmp;
                                v_sz = it_h->second.size();

                                v_min_b.push_back(minz_rep);
                            }
                            else break;
                        }
                    }
                }
                else v_min_b.push_back(minz_rep); //Add minimizer to list of minimizers

                last_pos_min = min_h_res.pos;
                it_it_min++;
            }
        }
    }

    sort(v_min_b.begin(), v_min_b.end()); //Sort the list of minimizers

    auto it_a = std::begin(v_min_a);
    auto it_a_end = std::end(v_min_a);

    //Remove Minimizers that we have already changed
    auto iter = std::remove_if (std::begin(v_min_b), std::end(v_min_b),
                                [&it_a, &it_a_end](Minimizer& minz) -> bool {
                                    while  (it_a != it_a_end && *it_a < minz) ++it_a;
                                    return (it_a != it_a_end && *it_a == minz);
                                });

    for (vector<Minimizer>::iterator it_v_min = v_min_b.begin(); it_v_min != iter; it_v_min++){ // Iterate over minimizers

        if ((it_v_min == v_min_b.begin()) || (*it_v_min != *(it_v_min-1))){ //If the minimizer is diff. from the previous one

            hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(*it_v_min); // Look for the minimizer in the hash table

            if (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                tiny_vector<size_t,tiny_vector_sz>& v_id_unitigs = it_h->second;
                tiny_vector<size_t,tiny_vector_sz>::iterator it_v_c = v_id_unitigs.begin();

                // Iterate over unitig IDs associated with this minimizer
                for (; it_v_c != v_id_unitigs.end(); it_v_c++){

                     // Swap the unitig ids;
                    if ((*it_v_c & mask) == shift_id_unitig_a) *it_v_c = shift_id_unitig_b | (*it_v_c & MASK_CONTIG_POS);
                }
            }
        }
    }
}


// use:  split, deleted = mapper.splitAllUnitigs()
// post: All unitigs with 1 coverage somewhere have been split where the coverage is 1
//       split is the number of unitigs splitted
//       deleted is the number of unitigs deleted
//       Now every unitig in mapper has coverage >= 2 everywhere
pair<size_t, size_t> CompactedDBG::splitAllUnitigs() {

    size_t i;
    size_t split = 0, deleted = 0;
    size_t v_kmers_sz = v_kmers.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t nxt_pos_insert = v_unitigs.size();

    for (h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); it++) {

        if (!it->second.isFull()){

            deleteUnitig(false, true, it.getHash());
            deleted++;
        }
    }

    for (i = 0; i < v_kmers_sz;) {

        if (!v_kmers[i].second.isFull()) {

            v_kmers_sz--;

            if (i != v_kmers_sz) swapUnitigs(true, i, v_kmers_sz);

            deleteUnitig(true, false, v_kmers_sz);

            deleted++;
        }
        else i++;
    }

    for (i = 0; i < v_unitigs_sz;) { // Iterate over unitigs created so far

        if (!v_unitigs[i]->ccov.isFull()) { //Coverage not full, unitig must be splitted

            vector<pair<int,int>> sp = v_unitigs[i]->ccov.splittingVector();

            if (splitUnitig(i, nxt_pos_insert, v_unitigs_sz, v_kmers_sz, sp)) deleted++;
            else {

                split++;
                sp.clear();
            }
        }
        else i++;
    }

    if (nxt_pos_insert < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert);
    if (v_kmers_sz < v_kmers.size()) v_kmers.resize(v_kmers_sz);

    return make_pair(split, deleted);
}

bool CompactedDBG::splitUnitig(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs, size_t& v_unitigs_sz, size_t& v_kmers_sz,
                               const vector<pair<int,int>>& sp){

    Unitig* unitig = v_unitigs[pos_v_unitigs];

    bool first_long_unitig = true;
    bool deleted = true;

    int k = Kmer::k;

    if (!sp.empty()){

        pair<size_t, size_t> lowpair = unitig->ccov.lowCoverageInfo();

        size_t totalcoverage = unitig->coveragesum - lowpair.second;
        size_t ccov_size = unitig->ccov.size();

        const string str = unitig->seq.toString();

        int k = Kmer::k;

        for (vector<pair<int,int>>::const_iterator sit = sp.begin(); sit != sp.end(); ++sit) { //Iterate over created split unitigs

            size_t pos = sit->first;
            size_t len = sit->second - pos;

            string split_str = str.substr(pos, len + k - 1); // Split unitig sequence
            uint64_t cov_tmp = (totalcoverage * len) / (ccov_size - lowpair.first); // Split unitig coverage

            if (split_str.length() == k){

                if (addUnitig(split_str, v_kmers_sz)) h_kmers_ccov.find(Kmer(split_str.c_str()).rep())->second.setFull();
                else {

                    v_kmers[v_kmers_sz].second.setFull(); // We don't care about the coverage per k-mer anymore
                    v_kmers_sz++;
                }
            }
            else if (first_long_unitig){

                // The unitig is deleted but its space in the unitig vector is not because:
                // 1 - It would change indices in the minimizer hash table
                // 2 - It is going to be reused for one split unitig (more efficient than deleting)
                deleteUnitig(false, false, pos_v_unitigs);

                addUnitig(split_str, pos_v_unitigs);

                v_unitigs[pos_v_unitigs]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                v_unitigs[pos_v_unitigs]->coveragesum = cov_tmp;

                first_long_unitig = false;
            }
            else {

                addUnitig(split_str, nxt_pos_insert_v_unitigs);

                v_unitigs[nxt_pos_insert_v_unitigs]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                v_unitigs[nxt_pos_insert_v_unitigs]->coveragesum = cov_tmp;

                nxt_pos_insert_v_unitigs++;
            }
        }

        deleted = false;
    }

    if (first_long_unitig){

        nxt_pos_insert_v_unitigs--; //Position of the last unitig in the vector which is not NULL

        if (pos_v_unitigs != nxt_pos_insert_v_unitigs){ // Do not proceed to swap if swap positions are the same

            swapUnitigs(false, pos_v_unitigs, nxt_pos_insert_v_unitigs); // Swap unitigs

            // If the swapped unitig, previously in position nxt_pos_insert, was a split unitig
            // created in this method, do not try to split it again
            if (nxt_pos_insert_v_unitigs >= v_unitigs_sz) pos_v_unitigs++;
            else v_unitigs_sz--;
        }
        else v_unitigs_sz--;

        deleteUnitig(false, false, nxt_pos_insert_v_unitigs);
    }
    else pos_v_unitigs++;

    return deleted;
}

// use:  joined = mapper.joinAllUnitigs()
// pre:  no short unitigs exist in sUnitigs.
// post: all unitigs that could be connected have been connected
//       joined is the number of joined unitigs
size_t CompactedDBG::joinAllUnitigs(vector<Kmer>* v_joins) {

    size_t i, k = Kmer::k;
    size_t joined = 0;
    size_t v_unitigs_size = v_unitigs.size();
    size_t v_kmers_size = v_kmers.size();

    // a and b are candidates for joining
    KmerHashTable<Kmer> joins;

    if (v_joins == NULL){

        for (h_kmers_ccov_t::iterator it_ccov = h_kmers_ccov.begin(); it_ccov != h_kmers_ccov.end(); it_ccov++) {

            Kmer& tail = it_ccov->first;
            Kmer head_twin = tail.twin();
            Kmer fw, bw;

            const UnitigMap cm(it_ccov.getHash(), 0, 1, k, false, true, true, this);

            if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(make_pair(fw.twin(), tail));
            if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(make_pair(bw.twin(), head_twin));
        }

        for (i = 0; i != v_kmers_size; i++) {

            Kmer& tail = v_kmers[i].first;
            Kmer head_twin = tail.twin();
            Kmer fw, bw;

            const UnitigMap cm(i, 0, 1, k, true, false, true, this);

            if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(make_pair(fw.twin(), tail));
            if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(make_pair(bw.twin(), head_twin));
        }

        for (size_t i = 0; i != v_unitigs_size; i++) {

            const CompressedSequence& seq = v_unitigs[i]->seq;

            Kmer head_twin = seq.getKmer(0).twin();
            Kmer tail = seq.getKmer(seq.size()-k);
            Kmer fw, bw;

            const UnitigMap cm(i, 0, 1, seq.size(), false, false, true, this);

            if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(make_pair(fw.twin(), tail));
            if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(make_pair(bw.twin(), head_twin));
        }
    }
    else {

        Kmer fw;

        for (auto km : *v_joins){

            UnitigMap cm = find(km, true);

            if (!cm.isEmpty){

                if (!cm.isShort && !cm.isAbundant){

                    if ((cm.dist == 0 && cm.strand) || (cm.dist != 0 && !cm.strand)) km = km.twin();
                    if (checkJoin(km, cm, fw)) joins.insert(make_pair(fw.twin(), km));
                }
                else {

                    if (checkJoin(km, cm, fw)) joins.insert(make_pair(fw.twin(), km));

                    km = km.twin();

                    if (checkJoin(km, cm, fw)) joins.insert(make_pair(fw.twin(), km));
                }
            }
        }

        (*v_joins).clear();
    }

    for (KmerHashTable<Kmer>::iterator it = joins.begin(); it != joins.end(); it++) {

        const Kmer head = it->second;
        const Kmer tail = it->first.twin();

        UnitigMap cmHead = find(head, true);
        UnitigMap cmTail = find(tail, true);

        if (!cmHead.isEmpty && !cmTail.isEmpty) {

            Kmer cmHead_head, cmTail_head;

            if (cmHead.isShort) cmHead_head = v_kmers[cmHead.pos_unitig].first;
            else if (cmHead.isAbundant) cmHead_head = h_kmers_ccov.find(cmHead.pos_unitig)->first;
            else cmHead_head = v_unitigs[cmHead.pos_unitig]->seq.getKmer(0);

            if (cmTail.isShort) cmTail_head = v_kmers[cmTail.pos_unitig].first;
            else if (cmTail.isAbundant) cmTail_head = h_kmers_ccov.find(cmTail.pos_unitig)->first;
            else cmTail_head = v_unitigs[cmTail.pos_unitig]->seq.getKmer(0);

            if (cmHead_head != cmTail_head) { // can't join a sequence with itself, either hairPin, loop or mobius loop

                // both kmers are still end-kmers
                bool headDir;
                bool len_k_head = cmHead.isShort || cmHead.isAbundant;

                if (len_k_head && (head == cmHead_head)) headDir = true;
                else if (!len_k_head && (head == v_unitigs[cmHead.pos_unitig]->seq.getKmer(v_unitigs[cmHead.pos_unitig]->numKmers()-1))) headDir = true;
                else if (head.twin() == cmHead_head) headDir = false;
                else continue;

                bool tailDir;
                bool len_k_tail = cmTail.isShort || cmTail.isAbundant;

                if (tail == cmTail_head) tailDir = true;
                else if (len_k_tail){
                    if (tail.twin() == cmTail_head) tailDir = false;
                    else continue;
                }
                else if (tail.twin() == v_unitigs[cmTail.pos_unitig]->seq.getKmer(v_unitigs[cmTail.pos_unitig]->numKmers()-1)) tailDir = false;
                else continue;

                //Compute join sequence
                string joinSeq, tailSeq;

                if (headDir) joinSeq = len_k_head ? cmHead_head.toString() : v_unitigs[cmHead.pos_unitig]->seq.toString();
                else joinSeq = len_k_head ? cmHead_head.twin().toString() : v_unitigs[cmHead.pos_unitig]->seq.rev().toString();

                if (tailDir) tailSeq = len_k_tail ? cmTail_head.toString() : v_unitigs[cmTail.pos_unitig]->seq.toString();
                else tailSeq = len_k_tail ? cmTail_head.twin().toString() : v_unitigs[cmTail.pos_unitig]->seq.rev().toString();

                assert(joinSeq.substr(joinSeq.size()-k+1) == tailSeq.substr(0,k-1));

                joinSeq.append(tailSeq, k-1, string::npos);

                //Compute new coverage
                uint64_t covsum;

                if (len_k_head){

                    CompressedCoverage& ccov = cmHead.isShort ? v_kmers[cmHead.pos_unitig].second : h_kmers_ccov.find(cmHead.pos_unitig)->second;
                    covsum = (ccov.isFull() ? ccov.cov_full : ccov.covAt(0));
                }
                else covsum = v_unitigs[cmHead.pos_unitig]->coveragesum;

                if (len_k_tail){

                    CompressedCoverage& ccov = cmTail.isShort ? v_kmers[cmTail.pos_unitig].second : h_kmers_ccov.find(cmTail.pos_unitig)->second;
                    covsum += (ccov.isFull()? ccov.cov_full : ccov.covAt(0));
                }
                else covsum += v_unitigs[cmTail.pos_unitig]->coveragesum;

                Unitig* unitig;

                if (cmHead.isShort){ //If head is a short unitig, swap and delete it

                    v_kmers_size--;

                    if (cmHead.pos_unitig != v_kmers_size){

                        swapUnitigs(true, cmHead.pos_unitig, v_kmers_size);

                        // If the last unitig of the vector used for the swap was the tail
                        if (cmTail.isShort && (v_kmers_size == cmTail.pos_unitig)) cmTail.pos_unitig = cmHead.pos_unitig;
                    }

                    deleteUnitig(true, false, v_kmers_size);
                }
                else if (cmHead.isAbundant) deleteUnitig(false, true, cmHead.pos_unitig);

                if (cmTail.isShort){ //If they are short

                    v_kmers_size--;

                    if (cmTail.pos_unitig != v_kmers_size){

                        swapUnitigs(true, cmTail.pos_unitig, v_kmers_size);

                        if (cmHead.isShort && (v_kmers_size == cmHead.pos_unitig)) cmHead.pos_unitig = cmTail.pos_unitig;
                    }

                    deleteUnitig(true, false, v_kmers_size);
                }
                else if (cmTail.isAbundant) deleteUnitig(false, true, cmTail.pos_unitig);

                if (len_k_head && len_k_tail){

                    addUnitig(joinSeq, v_unitigs_size);
                    unitig = v_unitigs[v_unitigs_size];
                    v_unitigs_size++;
                }
                else if (len_k_head){

                    deleteUnitig(false, false, cmTail.pos_unitig);
                    addUnitig(joinSeq, cmTail.pos_unitig);
                    unitig = v_unitigs[cmTail.pos_unitig];
                }
                else {

                    if (!len_k_tail){

                        v_unitigs_size--;

                        if (cmTail.pos_unitig != v_unitigs_size){

                            swapUnitigs(false, cmTail.pos_unitig, v_unitigs_size);

                            if (v_unitigs_size == cmHead.pos_unitig) cmHead.pos_unitig = cmTail.pos_unitig;
                        }

                        deleteUnitig(false, false, v_unitigs_size);
                    }

                    deleteUnitig(false, false, cmHead.pos_unitig);
                    addUnitig(joinSeq, cmHead.pos_unitig);
                    unitig = v_unitigs[cmHead.pos_unitig];
                }

                unitig->coveragesum = covsum;
                if (covsum >= unitig->ccov.cov_full * unitig->numKmers()) unitig->ccov.setFull();

                joined++;
            }
        }
    }

    if (v_unitigs_size < v_unitigs.size()) v_unitigs.resize(v_unitigs_size);
    if (v_kmers_size < v_kmers.size()) v_kmers.resize(v_kmers_size);

    return joined;
}

bool CompactedDBG::checkJoin(const Kmer& a, const UnitigMap& cm_a, Kmer& b) {

    size_t k = Kmer::k, g = Minimizer::g;
    size_t fw_count = 0, bw_count = 0;

    Kmer fw_cand;

    UnitigMap cm_cand, cm_cand_tmp;

    //char km_tmp[k + 1];

    //a.forwardBase('A').toString(km_tmp);

    //uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
    //std::pair<uint64_t*, uint64_t*> block = bf.getBlock(it_min_h);

    //RepHash rep_h(k - 1), rep_h_cpy;
    //rep_h.init(km_tmp);

    for (size_t i = 0; i < 4; i++) {

        //rep_h_cpy = rep_h;
        //rep_h_cpy.extendFW(alpha[i]);

        //if (bf.contains(rep_h_cpy.hash(), block)){

            Kmer fw = a.forwardBase(alpha[i]);
            cm_cand_tmp = find(fw, true);

            if (!cm_cand_tmp.isEmpty) {

                fw_count++;
                if (fw_count > 1) break;
                fw_cand = fw;
                cm_cand = cm_cand_tmp;
            }
        //}
    }

    if (fw_count == 1) {

        Kmer cand_head, ac_head;

        if (cm_cand.isShort) cand_head = v_kmers[cm_cand.pos_unitig].first;
        else if (cm_cand.isAbundant) cand_head = h_kmers_ccov.find(cm_cand.pos_unitig)->first;
        else cand_head = v_unitigs[cm_cand.pos_unitig]->seq.getKmer(0);

        if (cm_a.isShort) ac_head = v_kmers[cm_a.pos_unitig].first;
        else if (cm_a.isAbundant) ac_head = h_kmers_ccov.find(cm_a.pos_unitig)->first;
        else ac_head = v_unitigs[cm_a.pos_unitig]->seq.getKmer(0);

        if (cand_head != ac_head) {

            Kmer fw_cpy = fw_cand.twin();

            //fw_cpy.forwardBase('A').toString(km_tmp);

            //it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            //block = bf.getBlock(it_min_h);

            //rep_h.init(km_tmp);

            for (size_t j = 0; j < 4; j++) {

                //rep_h_cpy = rep_h;
                //rep_h_cpy.extendFW(alpha[j]);

                //if (bf.contains(rep_h_cpy.hash(), block)){

                    Kmer fw = fw_cpy.forwardBase(alpha[j]);
                    cm_cand_tmp = find(fw, true);

                    if (!cm_cand_tmp.isEmpty) {

                        bw_count++;
                        if (bw_count > 1) break;
                    }
                //}
            }

            if (bw_count == 1) { // join up

                if (cand_head == fw_cand) {

                    b = fw_cand;
                    return true;
                }

                Kmer candLast;

                if (cm_cand.isShort || cm_cand.isAbundant) candLast = cand_head;
                else candLast = v_unitigs[cm_cand.pos_unitig]->seq.getKmer(v_unitigs[cm_cand.pos_unitig]->seq.size()-k);

                if (candLast.twin() == fw_cand) {

                    b = fw_cand;
                    return true;
                }

                return true;
            }
        }
    }

    return false;
}

void CompactedDBG::check_fp_tips(KmerHashTable<bool>& ignored_km_tips){

    uint64_t nb_real_short_tips = 0;

    size_t k = Kmer::k, g = Minimizer::g;

    size_t nxt_pos_insert_v_unitigs = v_unitigs.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t v_kmers_sz = v_kmers.size();

    char km_tmp[k + 1];

    vector<pair<int,int>> sp;

    for (KmerHashTable<bool>::iterator it = ignored_km_tips.begin(); it != ignored_km_tips.end(); it++) {

        Kmer km = it->first;

        UnitigMap cm = find(km, true); // Check if the (short) tip actually exists

        if (!cm.isEmpty){ // IF the tip exists

            nb_real_short_tips++;

            bool not_found = true;

            km.backwardBase('A').toString(km_tmp);

            //uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            //std::pair<uint64_t*, uint64_t*> block = bf.getBlock(it_min_h);

            //RepHash rep_h(k - 1), rep_h_cpy;
            //rep_h.init(km_tmp + 1);

            for (size_t i = 0; (i < 4) && not_found; ++i) {

                //rep_h_cpy = rep_h;
                //rep_h_cpy.extendBW(alpha[i]);

                //if (bf.contains(rep_h_cpy.hash(), block)){ // Query for all of its predecessors

                    km_tmp[0] = alpha[i];

                    UnitigMap cm_bw = find(Kmer(km_tmp));

                    if (!cm_bw.isEmpty && !cm_bw.isAbundant && !cm_bw.isShort){

                        if (cm_bw.strand) cm_bw.dist++;

                        if ((cm_bw.dist != 0) && (cm_bw.dist != cm_bw.size - k + 1)){

                            sp.push_back(make_pair(0, cm_bw.dist));
                            sp.push_back(make_pair(cm_bw.dist, cm_bw.size - k + 1));

                            splitUnitig(cm_bw.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

                            sp.clear();
                        }

                        not_found = false;
                    }
                //}
            }

            km.forwardBase('A').toString(km_tmp);

            //it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            //block = bf.getBlock(it_min_h);

            //rep_h.init(km_tmp);

            for (size_t i = 0; (i < 4) && not_found; ++i) {

                //rep_h_cpy = rep_h;
                //rep_h_cpy.extendFW(alpha[i]);

                //if (bf.contains(rep_h_cpy.hash(), block)){

                    km_tmp[k - 1] = alpha[i];

                    UnitigMap cm_fw = find(Kmer(km_tmp));

                    if (!cm_fw.isEmpty && !cm_fw.isAbundant && !cm_fw.isShort){

                        if (!cm_fw.strand) cm_fw.dist++;

                        if ((cm_fw.dist != 0) && (cm_fw.dist != cm_fw.size - k + 1)){

                            sp.push_back(make_pair(0, cm_fw.dist));
                            sp.push_back(make_pair(cm_fw.dist, cm_fw.size - k + 1));

                            splitUnitig(cm_fw.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

                            sp.clear();
                        }

                        not_found = false;
                    }
                //}
            }
        }
    }

    if (nxt_pos_insert_v_unitigs < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert_v_unitigs);
    if (v_kmers_sz < v_kmers.size()) v_kmers.resize(v_kmers_sz);
}

void CompactedDBG::writeGFA(string graphfilename) const {

    const size_t k = Kmer::k, g = Minimizer::g;
    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = v_kmers.size();

    size_t i, labelA, labelB, id = v_unitigs_sz + v_kmers_sz + 1;

    //char km_tmp[k + 1];

    //std::pair<uint64_t*, uint64_t*> block;

    UnitigMap cand;

    //RepHash rep_h(k - 1), rep_h_cpy;

    Unitig* unitig = NULL;

    ofstream graphfile;
    ostream graph(0);

    graphfile.open(graphfilename.c_str());
    graph.rdbuf(graphfile.rdbuf());
    assert(!graphfile.fail());

    KmerHashTable<size_t> idmap(h_kmers_ccov.size());

    // gfa header
    graph << "H\tVN:Z:1.0\n";

    for (labelA = 1; labelA <= v_unitigs_sz; labelA++) {

        unitig = v_unitigs[labelA - 1];

        graph << "S\t" << labelA << "\t" << unitig->seq.toString() << "\tLN:i:" <<
        unitig->seq.size() << "\tXC:i:" << unitig->coveragesum << "\n";
    }

    for (labelA = 1; labelA <= v_kmers_sz; labelA++) {

        const pair<Kmer, CompressedCoverage>& p = v_kmers[labelA - 1];

        graph << "S\t" << (labelA + v_unitigs_sz) << "\t" << p.first.toString() << "\tLN:i:" <<
        k << "\tXC:i:" << (p.second.isFull() ? p.second.cov_full : p.second.covAt(0)) << "\n";
    }

    for (h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); it++) {

        id++;
        idmap.insert({it->first, id});

        graph << "S\t" << id << "\t" << it->first.toString() << "\tLN:i:" << k << "\tXC:i:" <<
        (it->second.isFull() ? it->second.cov_full : it->second.covAt(0)) << "\n";
    }

    // We need to deal with the tail of long unitigs
    for (labelA = 1; labelA <= v_unitigs_sz; labelA++) {

        Unitig* unitig = v_unitigs[labelA - 1];

        Kmer head = unitig->seq.getKmer(0);

        //head.backwardBase('A').toString(km_tmp);

        //block = bf.getBlock(minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash());

        //rep_h.init(km_tmp + 1);

        for (i = 0; i < 4; ++i) {

            //rep_h_cpy = rep_h;
            //rep_h_cpy.extendBW(alpha[i]);

            //if (bf.contains(rep_h_cpy.hash(), block)){

                Kmer b = head.backwardBase(alpha[i]);
                cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    graph << "L\t" << labelA << "\t-\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            //}
        }

        Kmer tail = unitig->seq.getKmer(unitig->seq.size() - k);

        //tail.forwardBase('A').toString(km_tmp);

        //block = bf.getBlock(minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash());

        //rep_h.init(km_tmp);

        for (i = 0; i < 4; ++i) {

            //rep_h_cpy = rep_h;
            //rep_h_cpy.extendFW(alpha[i]);

            //if (bf.contains(rep_h_cpy.hash(), block)){

                Kmer b = tail.forwardBase(alpha[i]);
                cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    graph << "L\t" << labelA << "\t+\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            //}
        }
    }

    for (labelA = v_unitigs_sz + 1; labelA <= v_kmers_sz + v_unitigs_sz; labelA++) {

        const pair<Kmer, CompressedCoverage>& p = v_kmers[labelA - v_unitigs_sz - 1];

        //p.first.backwardBase('A').toString(km_tmp);

        //block = bf.getBlock(minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash());

        //rep_h.init(km_tmp + 1);

        for (i = 0; i < 4; ++i) {

            //rep_h_cpy = rep_h;
            //rep_h_cpy.extendBW(alpha[i]);

            //if (bf.contains(rep_h_cpy.hash(), block)){

                Kmer b = p.first.backwardBase(alpha[i]);
                cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    graph << "L\t" << labelA << "\t-\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            //}
        }

        //p.first.forwardBase('A').toString(km_tmp);

        //block = bf.getBlock(minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash());

        //rep_h.init(km_tmp);

        for (i = 0; i < 4; ++i) {

            //rep_h_cpy = rep_h;
            //rep_h_cpy.extendFW(alpha[i]);

            //if (bf.contains(rep_h_cpy.hash(), block)){

                Kmer b = p.first.forwardBase(alpha[i]);
                cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    graph << "L\t" << labelA << "\t+\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            //}
        }
    }

    for (KmerHashTable<size_t>::iterator it = idmap.begin(); it != idmap.end(); it++) {

        labelA = it->second;

        //it->first.backwardBase('A').toString(km_tmp);

        //block = bf.getBlock(minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash());

        //rep_h.init(km_tmp + 1);

        for (i = 0; i < 4; ++i) {

            //rep_h_cpy = rep_h;
            //rep_h_cpy.extendBW(alpha[i]);

            //if (bf.contains(rep_h_cpy.hash(), block)){

                Kmer b = it->first.backwardBase(alpha[i]);
                cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    graph << "L\t" << labelA << "\t-\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            //}
        }

        //it->first.forwardBase('A').toString(km_tmp);

        //block = bf.getBlock(minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash());

        //rep_h.init(km_tmp);

        for (i = 0; i < 4; ++i) {

            //rep_h_cpy = rep_h;
            //rep_h_cpy.extendFW(alpha[i]);

            //if (bf.contains(rep_h_cpy.hash(), block)){

                Kmer b = it->first.forwardBase(alpha[i]);
                cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    graph << "L\t" << labelA << "\t+\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            //}
        }
    }

    graphfile.close();
}

size_t CompactedDBG::removeUnitigs(bool rmIsolated, bool clipTips, vector<Kmer>& v){

    if (!rmIsolated && !clipTips) return 0;

    const bool rm_and_clip = rmIsolated && clipTips;

    const size_t k = Kmer::k, g = Minimizer::g;

    size_t v_unitigs_sz = v_unitigs.size();
    size_t v_kmers_sz = v_kmers.size();
    size_t i;
    size_t removed = 0;

    //uint64_t it_min_h;

    int64_t j;

    //char km_tmp[k + 1];

    const int lim = (clipTips ? 1 : 0);

    int nb_pred, nb_succ;

    //std::pair<uint64_t*, uint64_t*> block;

    Kmer km;

    //RepHash rep_h(k - 1), rep_h_cpy;

    Unitig* unitig = NULL;

    for (j = 0; j < v_unitigs_sz; j++) {

        Unitig* unitig = v_unitigs[j];

        if (unitig->numKmers() < k){

            Kmer head = unitig->seq.getKmer(0);

            //head.backwardBase('A').toString(km_tmp);
            //rep_h.init(km_tmp + 1);

            //it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            //block = bf.getBlock(it_min_h);
            nb_pred = 0;

            for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

                //rep_h_cpy = rep_h;
                //rep_h_cpy.extendBW(alpha[i]);

                if (/*bf.contains(rep_h_cpy.hash(), block) &&*/ !find(head.backwardBase(alpha[i]), true).isEmpty){

                    nb_pred++;
                    if (clipTips) km = head.backwardBase(alpha[i]);
                }
            }

            if (nb_pred <= lim){

                Kmer tail = unitig->seq.getKmer(unitig->seq.size() - k);

                //tail.forwardBase('A').toString(km_tmp);
                //rep_h.init(km_tmp);

                //it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
                //block = bf.getBlock(it_min_h);
                nb_succ = 0;

                for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                    //rep_h_cpy = rep_h;
                    //rep_h_cpy.extendFW(alpha[i]);

                    if (/*bf.contains(rep_h_cpy.hash(), block) && */!find(tail.forwardBase(alpha[i]), true).isEmpty){

                        nb_succ++;
                        if (clipTips) km = tail.forwardBase(alpha[i]);
                    }
                }

                if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))) { //Unitig is isolated

                    removed++;
                    v_unitigs_sz--;

                    if (j != v_unitigs_sz){

                        swapUnitigs(false, j, v_unitigs_sz),
                        j--;
                    }

                    if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
                }
            }
        }
    }

    for (j = 0; j < v_kmers_sz; j++) {

        const pair<Kmer, CompressedCoverage>& p = v_kmers[j];

        //p.first.backwardBase('A').toString(km_tmp);
        //rep_h.init(km_tmp + 1);

        //it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
        //block = bf.getBlock(it_min_h);
        nb_pred = 0;

        for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

            //rep_h_cpy = rep_h;
            //rep_h_cpy.extendBW(alpha[i]);

            if (/*bf.contains(rep_h_cpy.hash(), block) &&*/ !find(p.first.backwardBase(alpha[i]), true).isEmpty){

                nb_pred++;
                if (clipTips) km = p.first.backwardBase(alpha[i]);
            }
        }

        if (nb_pred <= lim){

            //p.first.forwardBase('A').toString(km_tmp);
            //rep_h.init(km_tmp);

            //it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            //block = bf.getBlock(it_min_h);
            nb_succ = 0;

            for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                //rep_h_cpy = rep_h;
                //rep_h_cpy.extendFW(alpha[i]);

                if (/*bf.contains(rep_h_cpy.hash(), block) &&*/ !find(p.first.forwardBase(alpha[i]), true).isEmpty){

                    nb_succ++;
                    if (clipTips) km = p.first.forwardBase(alpha[i]);
                }
            }

            if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))) { //Unitig is isolated

                removed++;
                v_kmers_sz--;

                if (j != v_kmers_sz){

                    swapUnitigs(true, j, v_kmers_sz),
                    j--;
                }

                if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
            }
        }
    }

    for (h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); it++) {

        //it->first.backwardBase('A').toString(km_tmp);
        //rep_h.init(km_tmp + 1);

        //it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
        //block = bf.getBlock(it_min_h);
        nb_pred = 0;

        for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

            //rep_h_cpy = rep_h;
            //rep_h_cpy.extendBW(alpha[i]);

            if (/*bf.contains(rep_h_cpy.hash(), block) &&*/ !find(it->first.backwardBase(alpha[i]), true).isEmpty){

                nb_pred++;
                if (clipTips) km = it->first.backwardBase(alpha[i]);
            }
        }

        if (nb_pred <= lim){

            //it->first.forwardBase('A').toString(km_tmp);
            //rep_h.init(km_tmp);

            //it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            //block = bf.getBlock(it_min_h);
            nb_succ = 0;

            for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                //rep_h_cpy = rep_h;
                //rep_h_cpy.extendFW(alpha[i]);

                if (/*bf.contains(rep_h_cpy.hash(), block) &&*/ !find(it->first.forwardBase(alpha[i]), true).isEmpty){

                    nb_succ++;
                    if (clipTips) km = it->first.forwardBase(alpha[i]);
                }
            }

            if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))){

                removed++;

                it->second = CompressedCoverage();

                if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
            }
        }
    }

    for (j = v_unitigs_sz; j < v_unitigs.size(); j++) deleteUnitig(false, false, j);
    v_unitigs.resize(v_unitigs_sz);

    for (j = v_kmers_sz; j < v_kmers.size(); j++) deleteUnitig(true, false, j);
    v_kmers.resize(v_kmers_sz);

    for (h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); it++){

        if (it->second.size() == 0) deleteUnitig(false, true, it.getHash());
    }

    return removed;
}

size_t CompactedDBG::cstrMatch(const char* a, const char* b) const {

    char* _a = const_cast<char*>(a);
    char* _b = const_cast<char*>(b);

    while ((*_a != '\0') && (*_a == *_b)){ _a++; _b++; }

    return _a - a;
}

CompactedDBG::iterator CompactedDBG::begin() {

    CompactedDBG::iterator it(this);
    it++;
    return it;
}

CompactedDBG::const_iterator CompactedDBG::begin() const {

    CompactedDBG::const_iterator it(this);
    it++;
    return it;
}

CompactedDBG::iterator CompactedDBG::end() { return CompactedDBG::iterator(); }
CompactedDBG::const_iterator CompactedDBG::end() const { return CompactedDBG::const_iterator(); }
