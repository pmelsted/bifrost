static const uint8_t bits[256] = {
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};

/** Compacted de Bruijn graph constructor (set up an empty Compacted dBG).
* @param kmer_length is the length k of k-mers used in the graph (each unitig is of length at least k).
* @param minimizer_length is the length g of minimizers (g < k) used in the graph.
* @return an empty Compacted de Bruijn graph.
*/
template<typename T>
CompactedDBG<T>::CompactedDBG(int kmer_length, int minimizer_length) :  k_(kmer_length), g_(minimizer_length), invalid(false),
                                                                        has_data(!is_void<T>::value) {

    if (kmer_length <= 0){

        cerr << "CompactedDBG::setK(): Length k of k-mers cannot be less than or equal to 0" << endl;
        invalid = true;
    }

    if (kmer_length >= MAX_KMER_SIZE){

        cerr << "CompactedDBG::setK(): Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        invalid = true;
    }

    if (minimizer_length <= 0){

        cerr << "CompactedDBG::setK(): Length g of minimizers cannot be less than or equal to 0" << endl;
        invalid = true;
    }

    if (minimizer_length >= MAX_KMER_SIZE){

        cerr << "CompactedDBG::setK(): Length g of minimizers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        invalid = true;
    }

    if (!invalid){

        Kmer::set_k(k_);
        Minimizer::set_g(g_);
    }
}

/** Compacted de Bruijn graph destructor.
*/
template<typename T>
CompactedDBG<T>::~CompactedDBG() {

    empty();
}

/** Clear the graph: empty the graph and reset its parameters.
*/
template<typename T>
void CompactedDBG<T>::clear(){

    k_ = 0;
    g_ = 0;

    invalid = true;

    empty();
}

/** Empty the graph (does not reset its parameters).
*/
template<typename T>
void CompactedDBG<T>::empty(){

    for (auto unitig : v_unitigs) delete unitig;

    v_unitigs.clear();
    v_kmers.clear();

    hmap_min_unitigs.clear();
    h_kmers_ccov.clear();

    bf.clear();
}

/** Build the Compacted de Bruijn graph.
* @param opt is a structure from which the members are parameters of this function. See CDBG_Build_opt.
* @return boolean indicating if the graph has been built successfully.
*/
template<typename T>
bool CompactedDBG<T>::build(const CDBG_Build_opt& opt){

    size_t max_threads = std::thread::hardware_concurrency();

    bool construct_finished = true;

    if (invalid){

        cerr << "CompactedDBG::build(): Graph is invalid and cannot be built" << endl;
        construct_finished = false;
    }

    if (opt.nb_threads <= 0){

        cerr << "CompactedDBG::build(): Number of threads cannot be less than or equal to 0" << endl;
        construct_finished = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "CompactedDBG::build(): Number of threads cannot exceed " << max_threads << "threads" << endl;
        construct_finished = false;
    }

    if (opt.nb_bits_unique_kmers_bf <= 0){

        cerr << "CompactedDBG::build(): Number of Bloom filter bits per unique k-mer cannot be less than or equal to 0" << endl;
        construct_finished = false;
    }

    if (!opt.reference_mode && (opt.nb_bits_non_unique_kmers_bf <= 0)){

        cerr << "CompactedDBG::build(): Number of Bloom filter bits per non unique k-mer cannot be less than or equal to 0" << endl;
        construct_finished = false;
    }

    if (!opt.reference_mode && (opt.nb_non_unique_kmers > opt.nb_unique_kmers)){

        cerr << "CompactedDBG::build(): The estimated number of non unique k-mers ";
        cerr << "cannot be greater than the estimated number of unique k-mers" << endl;
        construct_finished = false;
    }

    if (opt.read_chunksize <= 0){

        cerr << "CompactedDBG::build(): Chunk size of reads to share among threads cannot be less than or equal to 0" << endl;
        construct_finished = false;
    }

    if (opt.unitig_size <= 0){

        cerr << "CompactedDBG::build(): Maximum unitig size cannot be less than or equal to 0" << endl;
        construct_finished = false;
    }

    if (opt.fastx_filename_in.size() == 0){

        cerr << "CompactedDBG::build(): Number of FASTA/FASTQ files in input cannot be less than or equal to 0" << endl;
        construct_finished = false;
    }

    if (opt.inFilenameBBF.length() != 0){

        FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

        if (fp == NULL) {

            cerr << "CompactedDBG::build(): Could not open input Blocked Bloom filter: " << opt.inFilenameBBF << endl;
            construct_finished = false;
        }
        else fclose(fp);
    }

    if (opt.outFilenameBBF.length() != 0){

        FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

        if (fp == NULL) {

            cerr << "CompactedDBG::build(): Could not open file for writing output Blocked Bloom filter: " << opt.outFilenameBBF << endl;
            construct_finished = false;
        }
        else fclose(fp);
    }

    for (vector<string>::const_iterator it = opt.fastx_filename_in.begin(); it != opt.fastx_filename_in.end(); ++it){

        FILE* fp = fopen(it->c_str(), "r");

        if (fp == NULL) {

            cerr << "CompactedDBG::build(): Could not open input FASTA/FASTQ file " << *it << endl;
            construct_finished = false;
        }
        else fclose(fp);
    }

    empty();

    if (opt.reference_mode) CompressedCoverage::setFullCoverage(1);

    if (opt.inFilenameBBF.length() != 0){

        FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

        if (fp == NULL) {

            cerr << "CompactedDBG::build(): Could not open input Blocked Bloom filter: " << opt.inFilenameBBF << endl;
            construct_finished = false;
        }
        else {

            construct_finished = bf.ReadBloomFilter(fp);

            fclose(fp);
        }
    }
    else construct_finished = filter(opt);

    if (construct_finished){

        if (opt.outFilenameBBF.length() != 0){

            FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

            if (fp == NULL) {

                cerr << "CompactedDBG::build(): Could not open file for writing output Blocked Bloom filter: " << opt.outFilenameBBF << endl;
                construct_finished = false;
            }
            else {

                bf.WriteBloomFilter(fp);

                fclose(fp);
            }
        }

        if (construct_finished){

            if (g_ == DEFAULT_G) hmap_min_unitigs = hmap_min_unitigs_t(opt.nb_non_unique_kmers / 4);

            // Construction step
            construct_finished = construct(opt);
        }

        bf.clear();
    }

    return construct_finished;
}

/** Simplify the Compacted de Bruijn graph: clip short (less than 2k length) tips and/or delete short (less than 2k length) isolated unitigs.
* @param delete_short_isolated_unitigs is a boolean indicating short isolated unitigs must be removed.
* @param clip_short_tips is a boolean indicating short tips must be clipped.
* @param verbose is a boolean indicating if information messages must be printed during the execution of the function.
* @return boolean indicating if the graph has been simplified successfully.
*/
template<typename T>
bool CompactedDBG<T>::simplify(const bool delete_short_isolated_unitigs, const bool clip_short_tips, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::simplify(): Graph is invalid and cannot be simplified" << endl;
        return false;
    }

    if (delete_short_isolated_unitigs || clip_short_tips){

        if (verbose) cout << endl << "CompactedDBG::simplify(): Removing isolated unitigs and/or clipping tips" << endl;

        vector<Kmer> v_joins;
        size_t joined = 0;

        size_t removed = removeUnitigs(delete_short_isolated_unitigs, clip_short_tips, v_joins);

        if (clip_short_tips) joined = joinUnitigs(&v_joins);

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

/** Write the Compacted de Bruijn graph to disk (GFA1 format).
* @param output_filename is a string containing the name of the file in which the graph will be written.
* @param nb_threads is a number indicating how many threads can be used to write the graph to disk.
* @param GFA_output indicates if the graph will be output in GFA format (true) or FASTA format (false).
* @param verbose is a boolean indicating if information messages must be printed during the execution of the function.
* @return boolean indicating if the graph has been written successfully.
*/
template<typename T>
bool CompactedDBG<T>::write(const string output_filename, const size_t nb_threads, const bool GFA_output, const bool verbose) {

    if (invalid){

        cerr << "CompactedDBG::write(): Graph is invalid and cannot be written to disk" << endl;
        return false;
    }

    if (nb_threads > std::thread::hardware_concurrency()){

        cerr << "CompactedDBG::write(): Number of threads cannot exceed " << std::thread::hardware_concurrency() << "threads" << endl;
        return false;
    }

    if (verbose) cout << endl << "CompactedDBG::write(): Writing graph to GFA file" << endl;

    string out = output_filename + (GFA_output ? ".gfa" : ".fasta");

    FILE* fp = fopen(out.c_str(), "w");

    if (fp == NULL) {

        cerr << "CompactedDBG::write(): Could not open file " << out << " for writing graph" << endl;
        return false;
    }
    else fclose(fp);

    GFA_output ? writeGFA(out, nb_threads) : writeFASTA(out);

    return true;
}

/** Find the unitig containing the queried k-mer in the Compacted de Bruijn graph.
* @param km is the queried k-mer (see Kmer class).
* @param extremities_only is a boolean indicating if the k-mer must be searched only in the unitig heads and tails (extremities_only = true).
* By default, the k-mer is searched everywhere (extremities_only = false) but is is slightly slower than looking only in the unitig heads and tails.
* @return UnitigMap<T> object containing the information of the unitig containg the queried k-mer (if present).
*/
template<typename T>
UnitigMap<T> CompactedDBG<T>::find(const Kmer& km, const bool extremities_only) {

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    char km_tmp[k_ + 1];
    km.toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k_, k_, g_, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz = Minimizer(&km_tmp[min_h_res.pos]).rep();
        hmap_min_unitigs_const_iterator it = hmap_min_unitigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            //const tiny_vector<size_t,tiny_vector_sz>& v = it->second;
            const tiny_vector<size_t,tiny_vector_sz>& v = *it;

            const int v_sz = v.size();

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v[i];
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()) return UnitigMap<T>(it_km.getHash(), 0, 1, k_, false, true, km == km_rep, *this);
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

                            if (v_kmers[unitig_id].first == km_rep) return UnitigMap<T>(unitig_id, 0, 1, k_, true, false, true, *this);
                        }
                        else if ((min_h_res.pos == diff - unitig_id_pos) && (v_kmers[unitig_id].first == km_rep)){

                            return UnitigMap<T>(unitig_id, 0, 1, k_, true, false, false, *this);
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - min_h_res.pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km)){

                                return UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, true, *this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km_twin)){

                                return UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, false, *this);
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km)){

                                return UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, true, *this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km_twin)){

                                return UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, false, *this);
                            }
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return UnitigMap<T>();
}

template<typename T>
vector<UnitigMap<T>> CompactedDBG<T>::findPredecessors(const Kmer& km, const bool extremities_only) {

    const Kmer km_pred[4] = {km.backwardBase('A'), km.backwardBase('C'), km.backwardBase('G'), km.backwardBase('T')};
    const Kmer km_rep[4] = {km_pred[0].rep(), km_pred[1].rep(), km_pred[2].rep(), km_pred[3].rep()};

    const Kmer km_twin_a = km_pred[0].twin();

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    char km_tmp[k_ + 1];
    km_pred[0].toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k_, k_, g_, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    vector<UnitigMap<T>> v_um(4);

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz = Minimizer(&km_tmp[min_h_res.pos]).rep();
        hmap_min_unitigs_const_iterator it = hmap_min_unitigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const tiny_vector<size_t,tiny_vector_sz>& v = *it;

            const int v_sz = v.size();

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v[i];
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km;

                        for (size_t j = 0; j != 4; ++j){

                            if ((it_km = h_kmers_ccov.find(km_rep[j])) != h_kmers_ccov.end()){

                                v_um[j] = UnitigMap<T>(it_km.getHash(), 0, 1, k_, false, true, km_pred[j] == km_rep[j], *this);
                            }
                        }
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

                        if ((min_h_res.pos == unitig_id_pos) || (min_h_res.pos == diff - unitig_id_pos)){

                            uint8_t idx = bits[v_kmers[unitig_id].first.getChar(0)];

                            if (v_kmers[unitig_id].first == km_rep[idx]) {

                                v_um[idx] = UnitigMap<T>(unitig_id, 0, 1, k_, true, false, min_h_res.pos == unitig_id_pos, *this);
                            }
                            else {

                                idx = 0x3 - bits[v_kmers[unitig_id].first.getChar(k_ - 1)];

                                if (v_kmers[unitig_id].first == km_rep[idx]) {

                                    v_um[idx] = UnitigMap<T>(unitig_id, 0, 1, k_, true, false, min_h_res.pos == unitig_id_pos, *this);
                                }
                            }
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - min_h_res.pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match + 1, k_ - 1, km)){

                                const uint8_t idx = bits[v_unitigs[unitig_id]->seq.getChar(pos_match)];
                                v_um[idx] = UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, true, *this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_ - 1, km_twin_a)){

                                const uint8_t idx = 0x3 - bits[v_unitigs[unitig_id]->seq.getChar(pos_match + k_ - 1)];
                                v_um[idx] = UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, false, *this);
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match + 1, k_ - 1, km)){

                                const uint8_t idx = bits[v_unitigs[unitig_id]->seq.getChar(pos_match)];
                                v_um[idx] = UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, true, *this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_ - 1, km_twin_a)){

                                const uint8_t idx = 0x3 - bits[v_unitigs[unitig_id]->seq.getChar(pos_match + k_ - 1)];
                                v_um[idx] = UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, false, *this);
                            }
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return v_um;
}

template<typename T>
vector<UnitigMap<T>> CompactedDBG<T>::findSuccessors(const Kmer& km, const size_t limit, const bool extremities_only) {

    if (limit == 0) return vector<UnitigMap<T>>();

    const Kmer km_succ[4] = {km.forwardBase('A'), km.forwardBase('C'), km.forwardBase('G'), km.forwardBase('T')};
    const Kmer km_rep[4] = {km_succ[0].rep(), km_succ[1].rep(), km_succ[2].rep(), km_succ[3].rep()};

    const Kmer km_twin_a = km_succ[0].twin().forwardBase('A');

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    int nb_found = 0;

    const int diff = k_ - g_;

    char km_tmp[k_ + 1];
    km_succ[0].toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k_, k_, g_, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    vector<UnitigMap<T>> v_um(4);

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz = Minimizer(&km_tmp[min_h_res.pos]).rep();
        hmap_min_unitigs_const_iterator it = hmap_min_unitigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            //const tiny_vector<size_t,tiny_vector_sz>& v = it->second;
            const tiny_vector<size_t,tiny_vector_sz>& v = *it;

            const int v_sz = v.size();

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v[i];
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km;

                        for (size_t j = 0; j != 4; ++j){

                            if (v_um[j].isEmpty && ((it_km = h_kmers_ccov.find(km_rep[j])) != h_kmers_ccov.end())){

                                v_um[j] = UnitigMap<T>(it_km.getHash(), 0, 1, k_, false, true, km_succ[j] == km_rep[j], *this);
                                if (++nb_found == limit) return v_um;
                            }
                        }
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

                        if ((min_h_res.pos == unitig_id_pos) || (min_h_res.pos == diff - unitig_id_pos)){

                            uint8_t idx = bits[v_kmers[unitig_id].first.getChar(k_ - 1)];

                            if (v_um[idx].isEmpty && (v_kmers[unitig_id].first == km_rep[idx])){

                                v_um[idx] = UnitigMap<T>(unitig_id, 0, 1, k_, true, false, min_h_res.pos == unitig_id_pos, *this);
                                if (++nb_found == limit) return v_um;
                            }
                            else {

                                idx = 0x3 - bits[v_kmers[unitig_id].first.getChar(0)];

                                if (v_um[idx].isEmpty && (v_kmers[unitig_id].first == km_rep[idx])){

                                    v_um[idx] = UnitigMap<T>(unitig_id, 0, 1, k_, true, false, min_h_res.pos == unitig_id_pos, *this);
                                    if (++nb_found == limit) return v_um;
                                }
                            }
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - min_h_res.pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_ - 1, km_succ[0])){

                                const int idx = bits[v_unitigs[unitig_id]->seq.getChar(pos_match + k_ - 1)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx] = UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, true, *this);
                                    if (++nb_found == limit) return v_um;
                                }
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match + 1, k_ - 1, km_twin_a)){

                                const int idx = 0x3 - bits[v_unitigs[unitig_id]->seq.getChar(pos_match)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx] = UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, false, *this);
                                    if (++nb_found == limit) return v_um;
                                }
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_ - 1, km_succ[0])){

                                const int idx = bits[v_unitigs[unitig_id]->seq.getChar(pos_match + k_ - 1)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx] =UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, true, *this);
                                    if (++nb_found == limit) return v_um;
                                }

                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match + 1, k_ - 1, km_twin_a)){

                                const int idx = 0x3 - bits[v_unitigs[unitig_id]->seq.getChar(pos_match)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx] = UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, false, *this);
                                    if (++nb_found == limit) return v_um;
                                }
                            }
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return v_um;
}

/** Add a sequence to the Compacted de Bruijn graph. Non-{A,C,G,T} characters such as Ns are discarded. The function automatically breaks the sequence into unitigs if required.
* @param seq is a string containing the sequence to insert.
* @param verbose is a boolean indicating if information messages must be printed during the execution of the function.
* @return a boolean indicating if the sequence was successfully inserted in the graph.
*/
template<typename T>
bool CompactedDBG<T>::add(const string& seq, const bool verbose){

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

    for (KmerIterator it_km(seq.c_str()), it_km_end; it_km != it_km_end; ++it_km) { //non-ACGT char. are discarded

        const std::pair<Kmer, int>& p = *it_km;
        const UnitigMap<T> cm = findUnitig(p.first, seq, p.second);

        if (cm.isEmpty){

            for (size_t i = 0; i != 4; ++i){

                UnitigMap<T> cm_bw = find(p.first.backwardBase(alpha[i]));

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

            for (size_t i = 0; i != 4; ++i){

                UnitigMap<T> cm_fw = find(p.first.forwardBase(alpha[i]));

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

                v_kmers[v_kmers_sz].second.ccov.setFull();

                v_kmers_sz++;
                added++;
            }
            else h_kmers_ccov.find(p.first.rep())->second.ccov.setFull();

            v_joins.push_back(p.first);
        }
        else {

            mapRead(cm);
            it_km += cm.len - 1;
        }
    }

    if (nxt_pos_insert_v_unitigs < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert_v_unitigs);
    if (v_kmers_sz < v_kmers.size()) v_kmers.resize(v_kmers_sz);

    size_t joined = joinUnitigs(&v_joins);

    if (verbose){

        cout << "CompactedDBG::add(): Added " << added << " k-mers to the graph" << endl;
        cout << "CompactedDBG::add(): Joined " << joined << " k-mers from the graph" << endl;
    }

    return true;
}

/** Remove a unitig from the Compacted de Bruijn graph.
* @param um is a UnitigMap object containing the information of the unitig to remove from the graph.
* @param verbose is a boolean indicating if information messages must be printed during the execution of the function.
* @return a boolean indicating if the unitig was successfully removed from the graph.
*/
template<typename T>
bool CompactedDBG<T>::remove(const UnitigMap<T>& um, const bool verbose){

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

        const Unitig<T>* unitig = v_unitigs[um.pos_unitig];

        head = unitig->seq.getKmer(0);
        tail = unitig->seq.getKmer(unitig->numKmers() - 1);
    }

    for (size_t i = 0; i != 4; ++i){

        const Kmer bw = head.backwardBase(alpha[i]);

        if (!find(bw, true).isEmpty) v_km.push_back(bw);
    }

    for (size_t i = 0; i != 4; ++i){

        const Kmer fw = tail.forwardBase(alpha[i]);

        if (!find(fw, true).isEmpty) v_km.push_back(fw);
    }

    const size_t swap_position = (um.isShort ? v_kmers.size() : v_unitigs.size()) - 1;

    if (!um.isAbundant && (um.pos_unitig != swap_position)) swapUnitigs(um.isShort, um.pos_unitig, swap_position);

    deleteUnitig(um.isShort, um.isAbundant, um.pos_unitig);

    if (um.isShort) v_kmers.resize(swap_position);
    else if (!um.isAbundant) v_unitigs.resize(swap_position);

    joinUnitigs(&v_km);

    return true;
}

/** Create an iterator to the first unitig of the Compacted de Bruijn graph (unitigs are NOT sorted lexicographically).
* @return an iterator to the first unitig of the graph.
*/
template<typename T>
typename CompactedDBG<T>::iterator CompactedDBG<T>::begin() {

    iterator it(this);
    ++it;
    return it;
}


/** Create an constant iterator to the first unitig of the Compacted de Bruijn graph (unitigs are NOT sorted lexicographically).
* @return a constant iterator to the first unitig of the graph.
*/
template<typename T>
typename CompactedDBG<T>::const_iterator CompactedDBG<T>::begin() const {

    const_iterator it(this);
    ++it;
    return it;
}

/** Create an iterator to the "past-the-last" unitig of the Compacted de Bruijn graph (unitigs are NOT sorted lexicographically).
* @return an iterator to the "past-the-last" unitig of the graph.
*/
template<typename T>
typename CompactedDBG<T>::iterator CompactedDBG<T>::end() { return iterator(); }

/** Create a constant iterator to the "past-the-last" unitig of the Compacted de Bruijn graph (unitigs are NOT sorted lexicographically).
* @return a constant iterator to the "past-the-last" unitig of the graph.
*/
template<typename T>
typename CompactedDBG<T>::const_iterator CompactedDBG<T>::end() const { return const_iterator(); }

template<typename T>
bool CompactedDBG<T>::join(const UnitigMap<T>& um, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::join(): Graph is invalid and no unitigs can be joined in it" << endl;
        return false;
    }

    if (verbose) cout << "CompactedDBG::join(): Start joining unitigs" << endl;

    vector<Kmer> v_km;

    v_km.push_back(um.getHead());
    if (!um.isShort && !um.isAbundant) v_km.push_back(um.getTail());

    size_t joined = joinUnitigs(&v_km);

    if (verbose) cout << "CompactedDBG::join(): " << joined << " unitigs have been joined" << endl;

    return true;
}

template<typename T>
bool CompactedDBG<T>::join(const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::join(): Graph is invalid and no unitigs can be joined in it" << endl;
        return false;
    }

    if (verbose) cout << "CompactedDBG::join(): Start joining all unitigs of the graph" << endl;

    vector<Kmer> v_km;

    size_t joined = joinUnitigs(&v_km);

    if (verbose) cout << "CompactedDBG::join(): " << joined << " unitigs have been joined" << endl;

    return true;
}

template<typename T>
bool CompactedDBG<T>::filter(const CDBG_Build_opt& opt) {

    if (invalid){

        cerr << "CompactedDBG::filter(): Graph is invalid and it cannot be built" << endl;
        return false;
    }

    BlockedBloomFilter bf_tmp;

    if (opt.reference_mode){

        BlockedBloomFilter tmp(opt.nb_unique_kmers, opt.nb_bits_unique_kmers_bf);
        bf.get(tmp);
        tmp.clear();
    }
    else {

        BlockedBloomFilter tmp(opt.nb_unique_kmers, opt.nb_bits_unique_kmers_bf);
        bf_tmp.get(tmp);
        tmp.clear();

        BlockedBloomFilter tmp2(opt.nb_non_unique_kmers, opt.nb_bits_non_unique_kmers_bf);
        bf.get(tmp2);
        tmp2.clear();
    }

    bool done = false;

    const bool multi_threaded = (opt.nb_threads != 1);

    string s;

    size_t file_id = 0;

    uint64_t n_read = 0;

    atomic<uint64_t> num_kmers(0), num_ins(0);

    FastqFile FQ(opt.fastx_filename_in);

    vector<string> readv;

    // Main worker thread
    auto worker_function = [&](vector<string>::iterator a, vector<string>::iterator b) {

        uint64_t l_num_kmers = 0, l_num_ins = 0;

        // for each input
        for (auto x = a; x != b; ++x) {

            const char* str = x->c_str();
            const int len = x->length();

            KmerHashIterator<RepHash> it_kmer_h(str, len, k_), it_kmer_h_end;
            minHashIterator<RepHash> it_min(str, len, k_, g_, RepHash(), true);

            for (int last_pos = -1; it_kmer_h != it_kmer_h_end; ++it_kmer_h, ++it_min, ++l_num_kmers) {

                std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                // If one or more k-mer were jumped because contained non-ACGT char.
                if (p_.second != last_pos + 1)
                    it_min = minHashIterator<RepHash>(&str[p_.second], len - p_.second, k_, g_, RepHash(), true);

                last_pos = p_.second;
                const uint64_t min_hr = it_min.getHash();

                if (opt.reference_mode){

                    bf.insert(p_.first, min_hr, multi_threaded);
                    ++l_num_ins;
                }
                else if (bf_tmp.insert(p_.first, min_hr, multi_threaded)) ++l_num_ins;
                else bf.insert(p_.first, min_hr, multi_threaded);
            }
        }

        // atomic adds
        num_kmers += l_num_kmers;
        num_ins += l_num_ins;
    };

    while (!done) {

        size_t reads_now = 0;

        while (reads_now < opt.read_chunksize) {

            if (FQ.read_next(s, file_id) >= 0) {

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
        size_t batch_size = readv.size() / opt.nb_threads;
        size_t leftover   = readv.size() % opt.nb_threads;

        for (size_t i = 0; i < opt.nb_threads; i++) {

            size_t jump = batch_size + ((i < leftover) ? 1 : 0);

            auto rit_end(rit);
            advance(rit_end, jump);
            workers.push_back(thread(worker_function, rit, rit_end));

            rit = rit_end;
        }

        assert(rit==readv.end());

        for (auto &t : workers) t.join();

        readv.clear();
    }

    FQ.close();

    if (opt.verbose) {

        cout << "CompactedDBG::filter(): Closed all fasta/fastq files" << endl;
        cout << "CompactedDBG::filter(): Processed " << num_kmers << " k-mers in " << n_read  << " reads" << endl;
        cout << "CompactedDBG::filter(): Found " << num_ins << " unique k-mers" << endl;
        cout << "CompactedDBG::filter(): Number of blocks in Bloom filter is " << bf.getNbBlocks() << endl;
    }

    if (opt.useMercyKmers && !opt.reference_mode){

        string mbbf_uniq_filename = opt.prefixFilenameOut + "_uniq";

        FILE* f_mbbf = fopen(mbbf_uniq_filename.c_str(), "wb");

        if (f_mbbf == NULL){

            cerr << "CompactedDBG::filter(): Minimizer Blocked Bloom filter of unique k-mers cannot be written to disk" << endl;
            return false;
        }

        bf_tmp.WriteBloomFilter(f_mbbf);

        fclose(f_mbbf);
    }

    return true;
}

template<typename T>
bool CompactedDBG<T>::construct(const CDBG_Build_opt& opt){

    if (invalid){

        cerr << "CompactedDBG::construct(): Graph is invalid and cannot be built" << endl;
        return false;
    }

    KmerIterator iter, iterend;

    FastqFile FQ(opt.fastx_filename_in);

    string s;

    size_t file_id = 0;

    uint64_t n_read = 0;

    vector<string> readv;

    const bool multi_threaded = (opt.nb_threads != 1);

    tiny_vector<Kmer, 2>* fp_candidate = opt.reference_mode ? nullptr : new tiny_vector<Kmer, 2>[bf.getNbBlocks()];

    KmerHashTable<bool> ignored_km_tips;

    auto worker_function = [&](vector<string>::const_iterator a, vector<string>::const_iterator b,
                                vector<NewUnitig>* smallv,
                                vector<tuple<bool, uint64_t, Kmer>>* v_fp_cand,
                                vector<Kmer>* l_ignored_km_tip) {

        uint64_t it_min_h, last_it_min_h;

        BlockedBloomFilter::BBF_Blocks block_bf;

        Kmer km;
        RepHash rep;

        // for each input
        for (auto x = a; x != b; ++x) {

            const char* s_x = x->c_str();

            KmerHashIterator<RepHash> it_kmer_h(s_x, x->length(), k_), it_kmer_h_end;
            preAllocMinHashIterator<RepHash> it_min(s_x, x->length(), k_, g_, rep, true);

            for (int last_pos_km = -2; it_kmer_h != it_kmer_h_end; ++it_kmer_h, ++it_min) {

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

                size_t r = bf.contains_block(p_.first, it_min_h, block_bf);

                if (r != 0){

                    UnitigMap<T> cm = findUnitig(km, *x, p_.second, it_min);

                    if (cm.isEmpty) { // kmer did not map, push into queue for next unitig generation round

                        bool selfLoop = false;
                        bool isIsolated = false;

                        string newseq;

                        size_t pos_match = findUnitigSequenceBBF(km, newseq, selfLoop, isIsolated, *l_ignored_km_tip); //Build unitig from Bloom filter

                        if (!opt.reference_mode && isIsolated){ // According to the BF, k-mer is isolated in the graph and is a potential false positive

                            //const uint64_t block = ((r == 1 ? block_bf.first : block_bf.second) - bf.getTable_ptr()) / NB_ELEM_BLOCK;
                            const uint64_t block = (r == 1 ? block_bf.first : block_bf.second);

                            Kmer km_rep = km.rep();
                            const tiny_vector<Kmer, 2>& v = fp_candidate[block];
                            tuple<bool, uint64_t, Kmer> t_fp_cand = make_tuple(true, block, km_rep);

                            for (const auto km_tmp : v){ // Go over global list of existing FP candidates

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

    /*const size_t chunk_size = 100;

    vector<vector<NewUnitig>> parray(opt.read_chunksize / chunk_size);
    vector<vector<tuple<bool, uint64_t, Kmer>>> v_fp_cand_threads(opt.read_chunksize / chunk_size);
    vector<vector<Kmer>> v_ignored_km_tip_thread(opt.read_chunksize / chunk_size);

    int round = 0;

    bool done = false;

    while (!done) {

        size_t reads_now = 0;

        while (reads_now < opt.read_chunksize) {

            if (FQ.read_next(s, file_id) >= 0){

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

        if (opt.verbose) cout << "CompactedDBG::construct(): Starting round " << round << endl;

        vector<thread> workers; // need to keep track of threads so we can join them
        queue<function<void()>> tasks; // the task queue

        // synchronization
        mutex queue_mutex;
        condition_variable condition;
        bool stop = false;

        for (size_t i = 0; i < opt.nb_threads; ++i){

            workers.emplace_back(
                [&tasks, &queue_mutex, &condition, &stop]
                {
                    for(;;)
                    {
                        function<void()> task;

                        {
                            unique_lock<mutex> lock(queue_mutex);
                            condition.wait(lock, [&stop, &tasks]{ return stop || !tasks.empty(); });
                            if (stop && tasks.empty()) return;
                            task = move(tasks.front());
                            tasks.pop();
                        }

                        task();
                    }
                }
            );
        }

        auto it_read = readv.begin();
        auto it_read_end = readv.end();

        uint64_t j = 0;

        while (it_read != it_read_end) {

            auto it_read_tmp(it_read);

            if (distance(it_read_tmp, it_read_end) >= chunk_size) advance(it_read_tmp, chunk_size);
            else it_read_tmp = it_read_end;

            std::function<void()> task = bind(worker_function, it_read, it_read_tmp, &parray[j], &v_fp_cand_threads[j], &v_ignored_km_tip_thread[j]);

            unique_lock<mutex> lock(queue_mutex);

            tasks.emplace(task);

            condition.notify_one();

            ++j;

            it_read = it_read_tmp;
        }

        {
            unique_lock<mutex> lock(queue_mutex);
            stop = true;
        }

        condition.notify_all();

        for(thread &worker: workers) worker.join();*/

    vector<vector<NewUnitig>> parray(opt.nb_threads);
    vector<vector<tuple<bool, uint64_t, Kmer>>> v_fp_cand_threads(opt.nb_threads);
    vector<vector<Kmer>> v_ignored_km_tip_thread(opt.nb_threads);

    int round = 0;

    bool done = false;

    while (!done) {

        size_t reads_now = 0;

        while (reads_now < opt.read_chunksize) {

            if (FQ.read_next(s, file_id) >= 0){

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

        if (opt.verbose) cout << "CompactedDBG::construct(): Starting round " << round << endl;

        // run parallel code
        vector<thread> workers;

        auto rit = readv.begin();
        size_t batch_size = readv.size() / opt.nb_threads;
        size_t leftover   = readv.size() % opt.nb_threads;

        for (size_t i = 0; i < opt.nb_threads; ++i) {

            size_t jump = batch_size + (i < leftover ? 1 : 0);
            auto rit_end(rit);

            advance(rit_end, jump);
            workers.push_back(thread(worker_function, rit, rit_end, &parray[i], &v_fp_cand_threads[i], &v_ignored_km_tip_thread[i]));

            rit = rit_end;
        }

        assert(rit == readv.end());

        for (auto& t : workers) t.join();

        for (auto& v : parray) { // for each thread

            for (const auto& x : v) addUnitigSequenceBBF(x.km, x.read, x.pos, x.seq, v_ignored_km_tip_thread[0]); // add each unitig for this thread

            v.clear(); //clear the map
        }

        for (auto& v : v_fp_cand_threads) { // for each thread

            for (const auto& x : v){ //for each FP candidate to add or delete

                Kmer km_fp = std::get<2>(x);

                tiny_vector<Kmer, 2>& v = fp_candidate[std::get<1>(x)];

                size_t i = 0;

                for (; i < v.size(); ++i){ //Search vector for k-mer

                    if (v[i] == km_fp) break;
                }

                if (std::get<0>(x) == false){ // If FP is to delete

                    if (i < v.size()) v.remove(i); // If k-mer has not already been deleted

                    UnitigMap<T> cm = find(km_fp); // Find FP to delete (was inserted into data structure as TP)

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

                    UnitigMap<T> cm = find(km_fp); // Find FP to delete (was inserted into data structure as TP)

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

            for (const auto x : v) ignored_km_tips.insert(x, false);

            v.clear();
        }

        if (opt.read_chunksize > 1 && opt.verbose) {

            cout << "CompactedDBG::construct(): End of round" << endl;
            cout << "CompactedDBG::construct(): Processed " << size() << " unitigs" << endl;
        }

        readv.clear();
    }

    FQ.close();

    bf.clear();

    if (fp_candidate != nullptr) delete[] fp_candidate;

    if (opt.verbose) cout << "CompactedDBG::construct(): Closed all fasta/fastq files" << endl;

    const size_t unitigsBefore = size();

    if (opt.verbose) cout << endl << "CompactedDBG::construct(): Splitting unitigs (1/" << (opt.reference_mode ? "1" : "2" ) << ")" << endl;

    pair<size_t, size_t> unitigSplit = splitAllUnitigs();

    const int unitigsAfter1 = size();

    //if (!opt.reference_mode){

        if (opt.verbose) cout << endl << "CompactedDBG::construct(): Splitting unitigs (2/2)" << endl;

        check_fp_tips(ignored_km_tips);
        ignored_km_tips.clear();
    //}

    const int unitigsAfter2 = size();

    if (opt.verbose) {

        cout << "CompactedDBG::construct(): Before split: " << unitigsBefore << " unitigs" << endl;
        cout << "CompactedDBG::construct(): After split (1/" << (opt.reference_mode ? "1" : "2" ) << "): " << unitigsAfter1 << " unitigs" <<  endl;
        if (!opt.reference_mode) cout << "CompactedDBG::construct(): After split (2/2): " << unitigsAfter2 << " unitigs" <<  endl;
        cout << "CompactedDBG::construct(): Unitigs split: " << unitigSplit.first << endl;
        cout << "CompactedDBG::construct(): Unitigs deleted: " << unitigSplit.second << endl;
    }

    if (opt.verbose) cout << endl << "CompactedDBG::construct(): Joining unitigs" << endl;

    const size_t joined = joinUnitigs(nullptr, opt.nb_threads);

    const int unitigsAfter3 = size();

    if (opt.verbose) {

        cout << "CompactedDBG::construct(): After join: " << unitigsAfter3 << " unitigs" << endl;
        cout << "CompactedDBG::construct(): Joined " << joined << " unitigs" << endl;
    }

    if (opt.useMercyKmers && !opt.reference_mode){

        string filename_mbbf_uniq_km = opt.prefixFilenameOut + "_uniq";

        joinTips(filename_mbbf_uniq_km, opt.nb_threads, opt.verbose);

        if (opt.verbose) cout << "CompactedDBG::construct(): After join tips using mercy k-mers: " << size() << " unitigs" << endl;

        if (std::remove(filename_mbbf_uniq_km.c_str()) != 0) {

            cerr << "CompactedDBG::construct(): Minimizer Blocked Bloom filter file of unique k-mers cannot be removed from disk" << endl;
        }
    }

    return true;
}

// use: b = cm.addUnitig(km,read)
// pre:
// post: either unitig string containsin has been added and b == true
//       or it was present and the coverage information was updated, b == false
//       NOT Threadsafe!
//bool UnitigMap<T>per::addUnitigSequenceBBF(Kmer km, const string& read, size_t pos, const string& seq) {
template<typename T>
bool CompactedDBG<T>::addUnitigSequenceBBF(Kmer km, const string& read, size_t pos, const string& seq, vector<Kmer>& l_ignored_km_tip) {

    string s;
    bool selfLoop = false;
    bool isIsolated = false;

    if (!seq.empty()) s = seq;
    else findUnitigSequenceBBF(km, s, selfLoop, isIsolated, l_ignored_km_tip);

    if (selfLoop) {

        bool foundAny = false;

        for (KmerIterator it(s.c_str()), it_end; it != it_end; ++it) {

            const UnitigMap<T> um = find(it->first);

            if (!um.isEmpty) {

                mapRead(um);
                foundAny = true;
            }
        }

        if (!foundAny) {

            addUnitig(s, s.length() == k_ ? v_kmers.size() : v_unitigs.size());

            for (KmerIterator it(s.c_str()), it_end; it != it_end; ++it) mapRead(find(it->first));
        }

        return true;
    }

    UnitigMap<T> cm = findUnitig(km, read, pos);

    if (cm.isEmpty){

        addUnitig(s, s.length() == k_ ? v_kmers.size() : v_unitigs.size());
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
template<typename T>
size_t CompactedDBG<T>::findUnitigSequenceBBF(Kmer km, string& s, bool& selfLoop, bool& isIsolated, vector<Kmer>& l_ignored_km_tip) {

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

        ++j;

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

            ++j;

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

    s.reserve(k_ + fw_s.size() + bw_s.size());
    s.append(bw_s);

    char tmp[Kmer::MAX_K];
    km.toString(tmp);
    s.append(tmp);

    s.append(fw_s);

    return bw_s.size();
}

template<typename T>
bool CompactedDBG<T>::bwStepBBF(Kmer km, Kmer& front, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand) const {

    size_t i, j = -1, j_tmp;
    size_t nb_neigh = 0;

    char km_tmp[k_ + 1];

    front.backwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash();

    RepHash rep_h(k_ - 1), rep_h_cpy;
    rep_h.init(km_tmp + 1);

    int found_fp_bw = 0;

    Kmer km_fp;

    bool pres_neigh_bw[4] = {false, false, false, false};
    uint64_t hashes_bw[4];

    for (i = 0; i != 4; ++i) {

        rep_h_cpy = rep_h;
        rep_h_cpy.extendBW(alpha[i]);

        hashes_bw[i] = rep_h_cpy.hash();
    }

    nb_neigh = bf.contains(hashes_bw, it_min_h, pres_neigh_bw, check_fp_cand ? 4 : 2);

    for (i = 0; i != 4; ++i) {

        if (pres_neigh_bw[i]){

            j = i;

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

                if (has_no_neighbor_tmp && fwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false)) ++found_fp_bw;
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
        uint64_t hashes_fw[4];

        Kmer bw = front.backwardBase(alpha[j]);

        bw.forwardBase('A').toString(km_tmp);

        it_min_h = minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash();

        rep_h.init(km_tmp);

        for (i = 0; i < 4; ++i) {

            rep_h_cpy = rep_h;
            rep_h_cpy.extendFW(alpha[i]);

            hashes_fw[i] = rep_h_cpy.hash();
        }

        nb_neigh = bf.contains(hashes_fw, it_min_h, pres_neigh_fw, 4);

        if (nb_neigh >= 2){

            for (i = 0; i < 4; ++i) {

                if (pres_neigh_fw[i]){

                    char dummy;
                    bool add = true, has_no_neighbor_tmp = false;

                    km_tmp[k_ - 1] = alpha[i];
                    km_fp = Kmer(km_tmp);

                    fwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                    if (has_no_neighbor_tmp && bwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false)){

                        if (km_fp != km) ++found_fp_fw;
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

                    km_tmp[k_ - 1] = alpha[i];
                    km_fp = Kmer(km_tmp).rep();

                    l_ignored_km_tip.push_back(km_fp);

                    --found_fp_fw;
                }
            }

            front.backwardBase('A').toString(km_tmp);

            for (i = 0; (i < 4) && (found_fp_bw != 0); ++i) {

                if (pres_neigh_bw[i]){

                    km_tmp[0] = alpha[i];
                    km_fp = Kmer(km_tmp).rep();

                    l_ignored_km_tip.push_back(km_fp);

                    --found_fp_bw;
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

template<typename T>
bool CompactedDBG<T>::fwStepBBF(Kmer km, Kmer& end, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand) const {

    size_t i, j = -1, j_tmp;
    size_t nb_neigh = 0;

    char km_tmp[k_ + 1];

    end.forwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash();

    RepHash rep_h(k_ - 1), rep_h_cpy;
    rep_h.init(km_tmp);

    int found_fp_fw = 0;
    Kmer km_fp;

    bool pres_neigh_fw[4] = {false, false, false, false};
    uint64_t hashes_fw[4];

    for (i = 0; i < 4; ++i) {

        rep_h_cpy = rep_h;
        rep_h_cpy.extendFW(alpha[i]);

        hashes_fw[i] = rep_h_cpy.hash();
    }

    nb_neigh = bf.contains(hashes_fw, it_min_h, pres_neigh_fw, check_fp_cand ? 4 : 2);

    for (i = 0; i < 4; ++i){

        if (pres_neigh_fw[i]){

            j = i;

            if (!check_fp_cand && (nb_neigh >= 2)) break;
        }
    }

    if (check_fp_cand && (nb_neigh >= 2)){

        for (i = 0; i < 4; ++i) {

            if (pres_neigh_fw[i]){

                km_tmp[k_ - 1] = alpha[i];
                km_fp = Kmer(km_tmp);

                char dummy;
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
        uint64_t hashes_bw[4];

        Kmer fw = end.forwardBase(alpha[j]);

        // check bw from fw link
        fw.backwardBase('A').toString(km_tmp);

        it_min_h = minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash();

        rep_h.init(km_tmp + 1);

        for (i = 0; i < 4; ++i) {

            rep_h_cpy = rep_h;
            rep_h_cpy.extendBW(alpha[i]);

            hashes_bw[i] = rep_h_cpy.hash();
        }

        nb_neigh = bf.contains(hashes_bw, it_min_h, pres_neigh_bw, 4);

        if (nb_neigh >= 2){

            for (i = 0; i < 4; ++i) {

                if (pres_neigh_bw[i]){

                    char dummy;
                    bool has_no_neighbor_tmp = false;

                    km_tmp[0] = alpha[i];
                    km_fp = Kmer(km_tmp);

                    bwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                    if (has_no_neighbor_tmp && fwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false)){

                        if (km_fp != km) ++found_fp_bw;
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

                    --found_fp_bw;
                }
            }

            end.forwardBase('A').toString(km_tmp);

            for (i = 0; (i < 4) && (found_fp_fw != 0); ++i) {

                if (pres_neigh_fw[i]){

                    km_tmp[k_ - 1] = alpha[i];
                    km_fp = Kmer(km_tmp).rep();

                    l_ignored_km_tip.push_back(km_fp);

                    --found_fp_fw;
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

// use:  cc = cm.findUnitig(km,s,pos)
// pre:  s[pos,pos+k-1] is the kmer km
// post: cc contains either the reference to the unitig position
//       or empty if none found
template<typename T>
UnitigMap<T> CompactedDBG<T>::findUnitig(const Kmer& km, const string& s, size_t pos) {

    // need to check if we find it right away, need to treat this common case
    UnitigMap<T> cc = find(km);

    if (!cc.isEmpty && !cc.isShort && !cc.isAbundant){

        const CompressedSequence& seq = v_unitigs[cc.pos_unitig]->seq;
        size_t km_dist = cc.dist;
        size_t jlen = 0;

        if (cc.strand) jlen = seq.jump(s.c_str(), pos, cc.dist, false) - k_ + 1;
        else {

            jlen = seq.jump(s.c_str(), pos, cc.dist + k_ - 1, true) - k_ + 1; // match s_fw to comp(seq)_bw
            km_dist -= jlen - 1;
        }

        return UnitigMap<T>(cc.pos_unitig, km_dist, jlen, cc.size, false, false, cc.strand, *this);
    }

    return cc;
}

template<typename T>
UnitigMap<T> CompactedDBG<T>::findUnitig(const Kmer& km, const string& s, size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h) {

    // need to check if we find it right away, need to treat this common case
    UnitigMap<T> cc = find(km, it_min_h);

    if (!cc.isEmpty && !cc.isShort && !cc.isAbundant){

        const CompressedSequence& seq = v_unitigs[cc.pos_unitig]->seq;
        size_t km_dist = cc.dist;
        size_t jlen = 0;

        if (cc.strand) jlen = seq.jump(s.c_str(), pos, cc.dist, false) - k_ + 1;
        else {

            jlen = seq.jump(s.c_str(), pos, cc.dist + k_ - 1, true) - k_ + 1; // match s_fw to comp(seq)_bw
            km_dist -= jlen - 1;
        }

        return UnitigMap<T>(cc.pos_unitig, km_dist, jlen, cc.size, false, false, cc.strand, *this);
    }

    return cc;
}

template<typename T>
bool CompactedDBG<T>::addUnitig(const string& str_unitig, const size_t id_unitig){

    int pos;

    size_t len = str_unitig.size();
    size_t pos_id_unitig = id_unitig << 32;

    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    bool isShort = false, isAbundant = false, isForbidden = false;

    char* c_str = const_cast<char*>(str_unitig.c_str());

    char km_tmp[k_ + 1];

    Kmer km_rep;

    if (len == k_){ // Unitig to add is short, maybe abundant as well

        isShort = true;

        pos_id_unitig |= MASK_CONTIG_TYPE;

        km_rep = Kmer(c_str).rep();
        km_rep.toString(km_tmp);
        c_str = km_tmp;
    }

    minHashIterator<RepHash> it_min(c_str, len, k_, g_, RepHash(), true), it_min_end;

    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){

        //If current minimizer was not seen before
        if ((last_pos_min < it_min.getPosition()) || isForbidden){

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&c_str[min_h_res.pos]).rep(); //Get the minimizer to insert
                std::pair<hmap_min_unitigs_iterator, bool> p = hmap_min_unitigs.insert(minz_rep, tiny_vector<size_t,tiny_vector_sz>());
                size_t v_sz = p.first->size();

                pos = min_h_res.pos;
                pos_id_unitig = (pos_id_unitig & mask) | ((size_t) pos);

                if (!isShort){

                    mhr = min_h_res;

                    while ((v_sz >= max_abundance_lim) || ((v_sz > 0) && (((*(p.first))[v_sz-1] & mask) == mask))){

                        mhr_tmp = it_min.getNewMin(mhr);
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            if (((*(p.first))[v_sz-1] & mask) != mask){ // Minimizer was never signaled before as overcrowded

                                // If minimizer bin already contains abundant k-mer, just set flag for unitig overcrowding
                                if (((*(p.first))[v_sz-1] & MASK_CONTIG_ID) == MASK_CONTIG_ID) (*(p.first))[v_sz-1] |= MASK_CONTIG_TYPE;
                                else p.first->push_back(mask);
                            }

                            mhr = mhr_tmp;
                            minz_rep = Minimizer(&c_str[mhr.pos]).rep();
                            p = hmap_min_unitigs.insert(minz_rep, tiny_vector<size_t,tiny_vector_sz>());
                            v_sz = p.first->size();
                        }
                        else break;
                    }
                }

                tiny_vector<size_t, tiny_vector_sz>& v = *(p.first);

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
                ++it_it_min;
            }
        }
    }

    if (isAbundant){

        if (id_unitig == v_kmers.size()) v_kmers.push_back(make_pair(km_rep, CompressedCoverage_t<T>(1)));
        else v_kmers[id_unitig] = make_pair(km_rep, CompressedCoverage_t<T>(1));

        deleteUnitig(true, false, id_unitig);
        if (id_unitig == v_kmers.size() - 1) v_kmers.resize(v_kmers.size() - 1);

        it_min = minHashIterator<RepHash>(c_str, len, k_, g_, RepHash(), true);

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){

            if (last_pos_min < it_min.getPosition()){ //If current minimizer was not seen before

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&c_str[min_h_res.pos]).rep(); //Get the minimizer to insert

                    std::pair<hmap_min_unitigs_iterator, bool> p = hmap_min_unitigs.insert(minz_rep, tiny_vector<size_t,tiny_vector_sz>());
                    tiny_vector<size_t,tiny_vector_sz>& v = *(p.first);

                    size_t v_sz = v.size();

                    if ((v_sz > 0) && ((v[v_sz-1] & MASK_CONTIG_ID) == MASK_CONTIG_ID)) v[v_sz-1]++;
                    else v.push_back(MASK_CONTIG_ID + 1);

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }

        h_kmers_ccov.insert(km_rep, CompressedCoverage_t<T>(1));
    }
    else if (isShort){

        if (id_unitig == v_kmers.size()) v_kmers.push_back(make_pair(km_rep,  CompressedCoverage_t<T>(1)));
        else v_kmers[id_unitig] = make_pair(km_rep,  CompressedCoverage_t<T>(1));
    }
    else if (id_unitig == v_unitigs.size()) v_unitigs.push_back(new Unitig<T>(c_str)); //Push unitig to list of unitigs
    else v_unitigs[id_unitig] = new Unitig<T>(c_str);

    return isAbundant;
}

template<typename T>
void CompactedDBG<T>::deleteUnitig(const bool isShort, const bool isAbundant, const size_t id_unitig){

    if (isAbundant){

        char km_str[k_ + 1];

        Kmer km = h_kmers_ccov.find(id_unitig).getKey();

        km.toString(km_str);

        minHashIterator<RepHash> it_min(km_str, k_, k_, g_, RepHash(), true), it_min_end;

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

            if (last_pos_min < it_min.getPosition()){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&km_str[min_h_res.pos]).rep(); // Get canonical minimizer

                    hmap_min_unitigs_iterator it_h = hmap_min_unitigs.find(minz_rep); // Look for the minimizer in the hash table

                    if (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        tiny_vector<size_t, tiny_vector_sz>& v = *it_h;

                        size_t last_pos_v = v.size() - 1;

                        v[last_pos_v]--;

                        if (((v[last_pos_v] & RESERVED_ID) == 0) && ((v[last_pos_v] & MASK_CONTIG_TYPE) != MASK_CONTIG_TYPE)){

                            if (last_pos_v == 0) hmap_min_unitigs.erase(minz_rep);
                            else v.remove(v.size() - 1);
                        }
                    }

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
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

    minHashIterator<RepHash> it_min(s, len, k_, g_, RepHash(), true), it_min_end;

    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

        if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer hash is found in unitig to delete

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep(); // Get canonical minimizer
                hmap_min_unitigs_iterator it_h = hmap_min_unitigs.find(minz_rep); // Look for the minimizer in the hash table

                mhr = min_h_res;

                while (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                    tiny_vector<size_t, tiny_vector_sz>& v = *it_h;

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
                ++it_it_min;
            }
        }
    }

    // The unitig is deleted but its space in the unitig vector is not because:
    // 1 - It would change indices in the minimizer hash table
    if (isShort) v_kmers[id_unitig].first.set_deleted();
    else {

        delete v_unitigs[id_unitig];
        v_unitigs[id_unitig] = nullptr;
    }
}

template<typename T>
void CompactedDBG<T>::swapUnitigs(const bool isShort, const size_t id_a, const size_t id_b){

    size_t shift_id_unitig_a = id_a << 32;
    size_t shift_id_unitig_b = id_b << 32;

    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    string str;

    size_t h_sz = 0;

    // Swap the unitig pointers in v_unitigs
    if (isShort){

        std::swap(v_kmers[id_a], v_kmers[id_b]);

        shift_id_unitig_a |= MASK_CONTIG_TYPE;
        shift_id_unitig_b |= MASK_CONTIG_TYPE;

        str = v_kmers[id_a].first.toString();
        h_sz = 4;
    }
    else {

        std::swap(v_unitigs[id_a], v_unitigs[id_b]);

        str = v_unitigs[id_a]->seq.toString();
        h_sz = (v_unitigs[id_a]->length() / 4) + (v_unitigs[id_b]->length() / 4);
    }


    MinimizerHashTable<uint8_t> minimizers(h_sz);

    bool isForbidden = false;

    // Swap the unitig IDs in the minimizer hash table
    const char* s = str.c_str();

    size_t len = str.size();

    minHashIterator<RepHash> it_min(s, len, k_, g_, RepHash(), true), it_min_end;
    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig

        if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer is found in unitig

            minHashResultIterator<RepHash> it_it_min(*it_min), it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep();
                hmap_min_unitigs_iterator it_h = hmap_min_unitigs.find(minz_rep);

                if (it_h != hmap_min_unitigs.end()){

                    if (minimizers.find(minz_rep) == minimizers.end()){

                        minimizers.insert(minz_rep, 0);

                        for (auto& unitig_id : *it_h){

                             // Swap the unitig ids but do not change positions;
                            if ((unitig_id & mask) == shift_id_unitig_b) unitig_id = shift_id_unitig_a | (unitig_id & MASK_CONTIG_POS);
                            else if ((unitig_id & mask) == shift_id_unitig_a) unitig_id = shift_id_unitig_b | (unitig_id & MASK_CONTIG_POS);
                        }
                    }

                    size_t v_sz = it_h->size();

                    mhr = min_h_res;

                    while (((*it_h)[v_sz-1] & mask) == mask){

                        mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            minz_rep = Minimizer(&s[mhr_tmp.pos]).rep();
                            it_h = hmap_min_unitigs.find(minz_rep);

                            if (it_h == hmap_min_unitigs.end()) break;

                            if (minimizers.find(minz_rep) == minimizers.end()){

                                minimizers.insert(minz_rep, 0);

                                for (auto& unitig_id : *it_h){

                                     // Swap the unitig ids but do not change positions;
                                    if ((unitig_id & mask) == shift_id_unitig_b) unitig_id = shift_id_unitig_a | (unitig_id & MASK_CONTIG_POS);
                                    else if ((unitig_id & mask) == shift_id_unitig_a) unitig_id = shift_id_unitig_b | (unitig_id & MASK_CONTIG_POS);
                                }
                            }

                            mhr = mhr_tmp;
                            v_sz = it_h->size();
                        }
                        else break;
                    }
                }

                last_pos_min = min_h_res.pos;
                ++it_it_min;
            }
        }
    }

    str = isShort ? v_kmers[id_b].first.toString() : v_unitigs[id_b]->seq.toString();
    s = str.c_str();
    len = str.size();

    isForbidden = false;

    minHashIterator<RepHash> it_min2(s, len, k_, g_, RepHash(), true);

    for (int64_t last_pos_min = -1; it_min2 != it_min_end; ++it_min2){ // Iterate over minimizers of unitig

        if ((last_pos_min < it_min2.getPosition()) || isForbidden){ // If a new minimizer is found in unitig

            minHashResultIterator<RepHash> it_it_min(*it_min2), it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep();
                hmap_min_unitigs_iterator it_h = hmap_min_unitigs.find(minz_rep);

                if (it_h != hmap_min_unitigs.end()){

                    if (minimizers.find(minz_rep) == minimizers.end()){

                        minimizers.insert(minz_rep, 0);

                        for (auto& unitig_id : *it_h){ // Swap unitig ids;

                            if ((unitig_id & mask) == shift_id_unitig_a) unitig_id = shift_id_unitig_b | (unitig_id & MASK_CONTIG_POS);
                        }
                    }

                    size_t v_sz = it_h->size();

                    mhr = min_h_res;

                    while (((*it_h)[v_sz-1] & mask) == mask){

                        mhr_tmp = it_min2.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            minz_rep = Minimizer(&s[mhr_tmp.pos]).rep();
                            it_h = hmap_min_unitigs.find(minz_rep);

                            if (it_h == hmap_min_unitigs.end()) break;

                            if (minimizers.find(minz_rep) == minimizers.end()){

                                minimizers.insert(minz_rep, 0);

                                for (auto& unitig_id : *it_h){ // Swap unitig ids;

                                    if ((unitig_id & mask) == shift_id_unitig_a) unitig_id = shift_id_unitig_b | (unitig_id & MASK_CONTIG_POS);
                                }
                            }

                            mhr = mhr_tmp;
                            v_sz = it_h->size();
                        }
                        else break;
                    }
                }

                last_pos_min = min_h_res.pos;
                ++it_it_min;
            }
        }
    }
}

/*
template<typename T>
bool CompactedDBG<T>::splitUnitig(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs, size_t& v_unitigs_sz, size_t& v_kmers_sz,
                                  const vector<pair<int,int>>& sp){

    bool first_long_unitig = true;
    bool deleted = true;

    if (!sp.empty()){

        Unitig<T>* unitig = v_unitigs[pos_v_unitigs];

        pair<size_t, size_t> lowpair = unitig->ccov.lowCoverageInfo();

        size_t totalcoverage = unitig->coveragesum - lowpair.second;
        size_t ccov_size = unitig->ccov.size();

        const string str = unitig->seq.toString();

        size_t i = 0;

        for (vector<pair<int,int>>::const_iterator sit = sp.begin(); sit != sp.end(); ++sit, ++i) { //Iterate over created split unitigs

            size_t pos = sit->first;
            size_t len = sit->second - pos;

            string split_str = str.substr(pos, len + k_ - 1); // Split unitig sequence
            uint64_t cov_tmp = (totalcoverage * len) / (ccov_size - lowpair.first); // Split unitig coverage

            Unitig<T> data_tmp;

            if (has_data){

                UnitigMap<T> um(pos_v_unitigs, 0, 1, unitig->length(), false, false, true, *this);

                data_tmp = um.splitData(pos, len); //Split the data
            }

            if (split_str.length() == k_){

                if (addUnitig(split_str, v_kmers_sz)){

                    CompressedCoverage_t<T>& cc_t = *h_kmers_ccov.find(Kmer(split_str.c_str()).rep());
                    cc_t.ccov.setFull();

                    if (has_data) cc_t.setData(data_tmp.getData());
                }
                else {

                    v_kmers[v_kmers_sz].second.ccov.setFull(); // We don't care about the coverage per k-mer anymore

                    if (has_data) v_kmers[v_kmers_sz].second.setData(data_tmp.getData());

                    ++v_kmers_sz;
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

                if (has_data) v_unitigs[pos_v_unitigs]->setData(data_tmp.getData());

                first_long_unitig = false;
            }
            else {

                addUnitig(split_str, nxt_pos_insert_v_unitigs);

                v_unitigs[nxt_pos_insert_v_unitigs]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                v_unitigs[nxt_pos_insert_v_unitigs]->coveragesum = cov_tmp;

                if (has_data) v_unitigs[nxt_pos_insert_v_unitigs]->setData(data_tmp.getData());

                ++nxt_pos_insert_v_unitigs;
            }
        }

        deleted = false;
    }

    if (first_long_unitig){

        --nxt_pos_insert_v_unitigs; //Position of the last unitig in the vector which is not NULL

        if (pos_v_unitigs != nxt_pos_insert_v_unitigs){ // Do not proceed to swap if swap positions are the same

            swapUnitigs(false, pos_v_unitigs, nxt_pos_insert_v_unitigs); // Swap unitigs

            // If the swapped unitig, previously in position nxt_pos_insert, was a split unitig
            // created in this method, do not try to split it again
            if (nxt_pos_insert_v_unitigs >= v_unitigs_sz) ++pos_v_unitigs;
            else --v_unitigs_sz;
        }
        else --v_unitigs_sz;

        deleteUnitig(false, false, nxt_pos_insert_v_unitigs);
    }
    else ++pos_v_unitigs;

    return deleted;
}
*/

template<>
inline bool CompactedDBG<void>::splitUnitig(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs, size_t& v_unitigs_sz, size_t& v_kmers_sz,
                                            const vector<pair<int,int>>& sp){
    bool first_long_unitig = true;
    bool deleted = true;

    if (!sp.empty()){

        Unitig<void>* unitig = v_unitigs[pos_v_unitigs];

        pair<size_t, size_t> lowpair = unitig->ccov.lowCoverageInfo();

        size_t totalcoverage = unitig->coveragesum - lowpair.second;
        size_t ccov_size = unitig->ccov.size();

        const string str = unitig->seq.toString();

        size_t i = 0;

        for (vector<pair<int,int>>::const_iterator sit = sp.begin(); sit != sp.end(); ++sit, ++i) { //Iterate over created split unitigs

            size_t pos = sit->first;
            size_t len = sit->second - pos;

            string split_str = str.substr(pos, len + k_ - 1); // Split unitig sequence
            uint64_t cov_tmp = (totalcoverage * len) / (ccov_size - lowpair.first); // Split unitig coverage

            if (split_str.length() == k_){

                if (addUnitig(split_str, v_kmers_sz)){

                    CompressedCoverage_t<void>& cc_t = *h_kmers_ccov.find(Kmer(split_str.c_str()).rep());
                    cc_t.ccov.setFull();
                }
                else {

                    v_kmers[v_kmers_sz].second.ccov.setFull(); // We don't care about the coverage per k-mer anymore

                    ++v_kmers_sz;
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

                ++nxt_pos_insert_v_unitigs;
            }
        }

        deleted = false;
    }

    if (first_long_unitig){

        --nxt_pos_insert_v_unitigs; //Position of the last unitig in the vector which is not NULL

        if (pos_v_unitigs != nxt_pos_insert_v_unitigs){ // Do not proceed to swap if swap positions are the same

            swapUnitigs(false, pos_v_unitigs, nxt_pos_insert_v_unitigs); // Swap unitigs

            // If the swapped unitig, previously in position nxt_pos_insert, was a split unitig
            // created in this method, do not try to split it again
            if (nxt_pos_insert_v_unitigs >= v_unitigs_sz) ++pos_v_unitigs;
            else --v_unitigs_sz;
        }
        else --v_unitigs_sz;

        deleteUnitig(false, false, nxt_pos_insert_v_unitigs);
    }
    else ++pos_v_unitigs;

    return deleted;
}

template<typename T>
bool CompactedDBG<T>::splitUnitig(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs, size_t& v_unitigs_sz, size_t& v_kmers_sz,
                                  const vector<pair<int,int>>& sp){

    bool long_unitig = false;
    bool deleted = true;

    if (!sp.empty()){

        Unitig<T>* unitig = v_unitigs[pos_v_unitigs];

        pair<size_t, size_t> lowpair = unitig->ccov.lowCoverageInfo();

        size_t totalcoverage = unitig->coveragesum - lowpair.second;
        size_t ccov_size = unitig->ccov.size();

        const string str = unitig->seq.toString();

        size_t i = 0;

        for (vector<pair<int,int>>::const_iterator sit = sp.begin(); sit != sp.end(); ++sit, ++i) { //Iterate over created split unitigs

            size_t pos = sit->first;
            size_t len = sit->second - pos;

            string split_str = str.substr(pos, len + k_ - 1); // Split unitig sequence
            uint64_t cov_tmp = (totalcoverage * len) / (ccov_size - lowpair.first); // Split unitig coverage

            UnitigMap<T> um(pos_v_unitigs, 0, 1, unitig->length(), false, false, true, *this);
            Unitig<T> data_tmp = um.splitData(pos, len); //Split the data

            if (split_str.length() == k_){

                if (addUnitig(split_str, v_kmers_sz)){

                    CompressedCoverage_t<T>& cc_t = *h_kmers_ccov.find(Kmer(split_str.c_str()).rep());

                    cc_t.ccov.setFull();
                    cc_t.setData(data_tmp.getData());
                }
                else {

                    v_kmers[v_kmers_sz].second.ccov.setFull(); // We don't care about the coverage per k-mer anymore
                    v_kmers[v_kmers_sz].second.setData(data_tmp.getData());

                    ++v_kmers_sz;
                }
            }
            else {

                addUnitig(split_str, nxt_pos_insert_v_unitigs);

                v_unitigs[nxt_pos_insert_v_unitigs]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                v_unitigs[nxt_pos_insert_v_unitigs]->coveragesum = cov_tmp;
                v_unitigs[nxt_pos_insert_v_unitigs]->setData(data_tmp.getData());

                ++nxt_pos_insert_v_unitigs;

                long_unitig = true;
            }
        }

        deleted = false;
    }

    //if (long_unitig){

        --nxt_pos_insert_v_unitigs; //Position of the last unitig in the vector which is not NULL

        if (pos_v_unitigs != nxt_pos_insert_v_unitigs){ // Do not proceed to swap if swap positions are the same

            swapUnitigs(false, pos_v_unitigs, nxt_pos_insert_v_unitigs); // Swap unitigs

            // If the swapped unitig, previously in position nxt_pos_insert, was a split unitig
            // created in this method, do not try to split it again
            if (nxt_pos_insert_v_unitigs >= v_unitigs_sz) ++pos_v_unitigs;
            else --v_unitigs_sz;
        }
        else --v_unitigs_sz;

        deleteUnitig(false, false, nxt_pos_insert_v_unitigs);
    //}
    //else ++pos_v_unitigs;

    return deleted;
}

template<typename T>
UnitigMap<T> CompactedDBG<T>::find(const Kmer& km, const preAllocMinHashIterator<RepHash>& it_min_h) {

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort;

    size_t unitig_id, unitig_id_pos_tmp, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    preAllocMinHashIterator<RepHash> it_min(it_min_h, k_);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz = Minimizer(&it_min.s[min_h_res.pos]).rep();
        hmap_min_unitigs_const_iterator it = hmap_min_unitigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const tiny_vector<size_t, tiny_vector_sz>& v = *it;

            it = hmap_min_unitigs.end();

            for (const auto unitig_id_pos : v){ // For each unitig associated with the minimizer

                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()){

                            return UnitigMap<T>(it_km.getHash(), 0, 1, k_, false, true, km == km_rep, *this);
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
                    unitig_id_pos_tmp = unitig_id_pos & MASK_CONTIG_POS;

                    if (isShort){

                        if (min_h_res.pos == unitig_id_pos_tmp){

                            if (v_kmers[unitig_id].first == km_rep){

                                return UnitigMap<T>(unitig_id, 0, 1, k_, true, false, true, *this);
                            }
                        }
                        else if ((min_h_res.pos == diff - unitig_id_pos_tmp) && (v_kmers[unitig_id].first == km_rep)){

                            return UnitigMap<T>(unitig_id, 0, 1, k_, true, false, false, *this);
                        }
                    }
                    else {

                        pos_match = unitig_id_pos_tmp - min_h_res.pos;
                        len = v_unitigs[unitig_id]->length() - k_;

                        if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km)){

                            return UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, true, *this);
                        }

                        pos_match = unitig_id_pos_tmp - diff + min_h_res.pos;

                        if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km_twin)){

                            return UnitigMap<T>(unitig_id, pos_match, 1, len + k_, false, false, false, *this);
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return UnitigMap<T>();
}

// use:  split, deleted = mapper.splitAllUnitigs()
// post: All unitigs with 1 coverage somewhere have been split where the coverage is 1
//       split is the number of unitigs splitted
//       deleted is the number of unitigs deleted
//       Now every unitig in mapper has coverage >= 2 everywhere
template<typename T>
pair<size_t, size_t> CompactedDBG<T>::splitAllUnitigs() {

    size_t i;
    size_t split = 0, deleted = 0;
    size_t v_kmers_sz = v_kmers.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t nxt_pos_insert = v_unitigs.size();

    for (typename h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); ++it) {

        if (!it->ccov.isFull()){

            deleteUnitig(false, true, it.getHash());
            deleted++;
        }
    }

    for (i = 0; i < v_kmers_sz;) {

        if (!v_kmers[i].second.ccov.isFull()) {

            --v_kmers_sz;

            if (i != v_kmers_sz) swapUnitigs(true, i, v_kmers_sz);

            deleteUnitig(true, false, v_kmers_sz);

            ++deleted;
        }
        else ++i;
    }

    for (i = 0; i < v_unitigs_sz;) { // Iterate over unitigs created so far

        if (!v_unitigs[i]->ccov.isFull()) { //Coverage not full, unitig must be splitted

            vector<pair<int,int>> sp = v_unitigs[i]->ccov.splittingVector();

            if (splitUnitig(i, nxt_pos_insert, v_unitigs_sz, v_kmers_sz, sp)) deleted++;
            else {

                ++split;
                sp.clear();
            }
        }
        else ++i;
    }

    if (nxt_pos_insert < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert);
    if (v_kmers_sz < v_kmers.size()) v_kmers.resize(v_kmers_sz);

    return make_pair(split, deleted);
}

// use:  joined = mapper.joinUnitigs()
// pre:  no short unitigs exist in sUnitigs.
// post: all unitigs that could be connected have been connected
//       joined is the number of joined unitigs
template<typename T>
size_t CompactedDBG<T>::joinUnitigs(vector<Kmer>* v_joins, const size_t nb_threads) {

    size_t i;
    size_t joined = 0;
    size_t cov_full = CompressedCoverage::getFullCoverage();
    size_t v_unitigs_size = v_unitigs.size();
    size_t v_kmers_size = v_kmers.size();

    const size_t chunk_size = 10000;

    // a and b are candidates for joining
    KmerHashTable<Kmer> joins;

    if (v_joins == nullptr){

        for (typename h_kmers_ccov_t::iterator it_ccov = h_kmers_ccov.begin(); it_ccov != h_kmers_ccov.end(); ++it_ccov) {

            Kmer tail = it_ccov.getKey();
            Kmer head_twin = tail.twin();
            Kmer fw, bw;

            const UnitigMap<T> cm(it_ccov.getHash(), 0, 1, k_, false, true, true, *this);

            if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail);
            if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin);
        }

        if (nb_threads == 1){

            for (i = 0; i != v_kmers_size; ++i) {

                Kmer tail = v_kmers[i].first;
                Kmer head_twin = tail.twin();
                Kmer fw, bw;

                const UnitigMap<T> cm(i, 0, 1, k_, true, false, true, *this);

                if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail);
                if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin);
            }

            for (size_t i = 0; i != v_unitigs_size; ++i) {

                const CompressedSequence& seq = v_unitigs[i]->seq;

                Kmer head_twin = seq.getKmer(0).twin();
                Kmer tail = seq.getKmer(seq.size() - k_);
                Kmer fw, bw;

                const UnitigMap<T> cm(i, 0, 1, seq.size(), false, false, true, *this);

                if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail);
                if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin);
            }
        }
        else {

            auto worker_v_kmers = [&joins, this](typename vector<pair<Kmer, CompressedCoverage_t<T>>>::const_iterator a,
                                                 typename vector<pair<Kmer, CompressedCoverage_t<T>>>::const_iterator b,
                                                 vector<pair<Kmer, Kmer>>* v_out){

                for (size_t i = a - v_kmers.begin(), end = b - v_kmers.begin(); i != end; ++i) {

                    Kmer tail = v_kmers[i].first;
                    Kmer head_twin = tail.twin();
                    Kmer fw, bw;

                    const UnitigMap<T> cm(i, 0, 1, k_, true, false, true, *this);

                    if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) v_out->push_back(make_pair(fw.twin(), tail));
                    if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) v_out->push_back(make_pair(bw.twin(), head_twin));
                }
            };

            auto worker_v_unitigs = [&joins, this](typename vector<Unitig<T>*>::const_iterator a,
                                                   typename vector<Unitig<T>*>::const_iterator b,
                                                   vector<pair<Kmer, Kmer>>* v_out){

                for (size_t i = a - v_unitigs.begin(), end = b - v_unitigs.begin(); i != end; ++i) {

                    const CompressedSequence& seq = v_unitigs[i]->seq;

                    Kmer head_twin = seq.getKmer(0).twin();
                    Kmer tail = seq.getKmer(seq.size() - k_);
                    Kmer fw, bw;

                    const UnitigMap<T> cm(i, 0, 1, seq.size(), false, false, true, *this);

                    if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) v_out->push_back(make_pair(fw.twin(), tail));
                    if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) v_out->push_back(make_pair(bw.twin(), head_twin));
                }
            };

            vector<vector<pair<Kmer, Kmer>>> v_out(nb_threads);

            auto it_kmer = v_kmers.begin();
            auto it_kmer_end = v_kmers.end();

            while (it_kmer != it_kmer_end) {

                vector<thread> workers;

                for (size_t i = 0; i != nb_threads; ++i) {

                    auto it_kmer_tmp(it_kmer);

                    if (distance(it_kmer_tmp, it_kmer_end) >= chunk_size){

                        advance(it_kmer_tmp, chunk_size);
                        workers.push_back(thread(worker_v_kmers, it_kmer, it_kmer_tmp, &v_out[i]));

                        it_kmer = it_kmer_tmp;
                    }
                    else {

                        workers.push_back(thread(worker_v_kmers, it_kmer, it_kmer_end, &v_out[i]));

                        it_kmer = it_kmer_end;

                        break;
                    }
                }

                for (auto &t : workers) t.join();

                for (size_t i = 0; i < nb_threads; ++i) {

                    for (const auto& p : v_out[i]) joins.insert(p.first, p.second);

                    v_out[i].clear();
                }
            }

            auto it_unitig = v_unitigs.begin();
            auto it_unitig_end = v_unitigs.end();

            while (it_unitig != it_unitig_end) {

                vector<thread> workers;

                for (size_t i = 0; i != nb_threads; ++i) {

                    auto it_unitig_tmp(it_unitig);

                    if (distance(it_unitig_tmp, it_unitig_end) >= chunk_size){

                        advance(it_unitig_tmp, chunk_size);
                        workers.push_back(thread(worker_v_unitigs, it_unitig, it_unitig_tmp, &v_out[i]));

                        it_unitig = it_unitig_tmp;
                    }
                    else {

                        workers.push_back(thread(worker_v_unitigs, it_unitig, it_unitig_end, &v_out[i]));

                        it_unitig = it_unitig_end;

                        break;
                    }
                }

                for (auto &t : workers) t.join();

                for (size_t i = 0; i < nb_threads; ++i) {

                    for (const auto& p : v_out[i]) joins.insert(p.first, p.second);

                    v_out[i].clear();
                }
            }
        }
    }
    else {

        Kmer fw;

        for (auto km : *v_joins){

            UnitigMap<T> cm = find(km, true);

            if (!cm.isEmpty){

                if (!cm.isShort && !cm.isAbundant){

                    if ((cm.dist == 0 && cm.strand) || (cm.dist != 0 && !cm.strand)) km = km.twin();
                    if (checkJoin(km, cm, fw)) joins.insert(fw.twin(), km);
                }
                else {

                    if (checkJoin(km, cm, fw)) joins.insert(fw.twin(), km);

                    km = km.twin();

                    if (checkJoin(km, cm, fw)) joins.insert(fw.twin(), km);
                }
            }
        }

        (*v_joins).clear();
    }

    for (KmerHashTable<Kmer>::iterator it = joins.begin(); it != joins.end(); ++it) {

        const Kmer head = *it;
        const Kmer tail = it.getKey().twin();

        UnitigMap<T> cmHead = find(head, true);
        UnitigMap<T> cmTail = find(tail, true);

        if (!cmHead.isEmpty && !cmTail.isEmpty) {

            Kmer cmHead_head, cmTail_head;

            if (cmHead.isShort) cmHead_head = v_kmers[cmHead.pos_unitig].first;
            else if (cmHead.isAbundant) cmHead_head = h_kmers_ccov.find(cmHead.pos_unitig).getKey();
            else cmHead_head = v_unitigs[cmHead.pos_unitig]->seq.getKmer(0);

            if (cmTail.isShort) cmTail_head = v_kmers[cmTail.pos_unitig].first;
            else if (cmTail.isAbundant) cmTail_head = h_kmers_ccov.find(cmTail.pos_unitig).getKey();
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

                assert(joinSeq.substr(joinSeq.size() - k_ + 1) == tailSeq.substr(0, k_ - 1));

                joinSeq.append(tailSeq, k_ - 1, string::npos);

                //Compute new coverage
                uint64_t covsum;

                if (len_k_head){

                    CompressedCoverage& ccov = cmHead.isShort ? v_kmers[cmHead.pos_unitig].second.ccov : h_kmers_ccov.find(cmHead.pos_unitig)->ccov;
                    covsum = ccov.covAt(0);
                }
                else covsum = v_unitigs[cmHead.pos_unitig]->coveragesum;

                if (len_k_tail){

                    CompressedCoverage& ccov = cmTail.isShort ? v_kmers[cmTail.pos_unitig].second.ccov : h_kmers_ccov.find(cmTail.pos_unitig)->ccov;
                    covsum += ccov.covAt(0);
                }
                else covsum += v_unitigs[cmTail.pos_unitig]->coveragesum;

                Unitig<T> data_tmp; //Store temporarily the new merged data
                Unitig<T>* unitig; //New unitig

                cmTail.strand = tailDir;
                cmHead.strand = headDir;

                if (cmHead.isShort || cmHead.isAbundant){

                    if (has_data){

                        cmTail.mergeData(cmHead); //If data, merge them
                        data_tmp.setData(cmTail.getData());
                    }

                    if (cmHead.isShort){ //If head is a short unitig, swap and delete it

                        --v_kmers_size;

                        if (cmHead.pos_unitig != v_kmers_size){

                            swapUnitigs(true, cmHead.pos_unitig, v_kmers_size);

                            // If the last unitig of the vector used for the swap was the tail
                            if (cmTail.isShort && (v_kmers_size == cmTail.pos_unitig)) cmTail.pos_unitig = cmHead.pos_unitig;
                        }

                        deleteUnitig(true, false, v_kmers_size);
                    }
                    else if (cmHead.isAbundant) deleteUnitig(false, true, cmHead.pos_unitig);
                }

                if (cmTail.isShort || cmTail.isAbundant){

                    if (has_data && !cmHead.isShort && !cmHead.isAbundant){

                        cmHead.mergeData(cmTail); //If data, merge them
                        data_tmp.setData(cmHead.getData());
                    }

                    if (cmTail.isShort){ //If tail is a short unitig, swap and delete it

                        --v_kmers_size;

                        if (cmTail.pos_unitig != v_kmers_size){

                            swapUnitigs(true, cmTail.pos_unitig, v_kmers_size);

                            if (cmHead.isShort && (v_kmers_size == cmHead.pos_unitig)) cmHead.pos_unitig = cmTail.pos_unitig;
                        }

                        deleteUnitig(true, false, v_kmers_size);
                    }
                    else if (cmTail.isAbundant) deleteUnitig(false, true, cmTail.pos_unitig);
                }

                if (len_k_head && len_k_tail){

                    addUnitig(joinSeq, v_unitigs_size);
                    unitig = v_unitigs[v_unitigs_size];
                    ++v_unitigs_size;
                }
                else if (len_k_head){

                    deleteUnitig(false, false, cmTail.pos_unitig);
                    addUnitig(joinSeq, cmTail.pos_unitig);
                    unitig = v_unitigs[cmTail.pos_unitig];
                }
                else {

                    if (!len_k_tail){

                        --v_unitigs_size;

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
                if (covsum >= cov_full * unitig->numKmers()) unitig->ccov.setFull();

                if (has_data) unitig->setData(data_tmp.getData()); //Remove the computed data if it exists

                ++joined;
            }
        }
    }

    if (v_unitigs_size < v_unitigs.size()) v_unitigs.resize(v_unitigs_size);
    if (v_kmers_size < v_kmers.size()) v_kmers.resize(v_kmers_size);

    return joined;
}

template<typename T>
bool CompactedDBG<T>::checkJoin(const Kmer& a, const UnitigMap<T>& cm_a, Kmer& b) {

    size_t i, j, count_succ;

    vector<UnitigMap<T>> v_um = findSuccessors(a, 2, true);

    for (i = 0, count_succ = 0; i != 4; ++i){

        if (!v_um[i].isEmpty){ ++count_succ; j = i; }
    }

    if (count_succ == 1) {

        Kmer cand_head, ac_head;
        const Kmer fw_cand = a.forwardBase(alpha[j]);
        const UnitigMap<T> cm_cand = v_um[j];

        if (cm_cand.isShort) cand_head = v_kmers[cm_cand.pos_unitig].first;
        else if (cm_cand.isAbundant) cand_head = h_kmers_ccov.find(cm_cand.pos_unitig).getKey();
        else cand_head = v_unitigs[cm_cand.pos_unitig]->seq.getKmer(0);

        if (cm_a.isShort) ac_head = v_kmers[cm_a.pos_unitig].first;
        else if (cm_a.isAbundant) ac_head = h_kmers_ccov.find(cm_a.pos_unitig).getKey();
        else ac_head = v_unitigs[cm_a.pos_unitig]->seq.getKmer(0);

        if (cand_head != ac_head) {

            v_um.clear();

            v_um = findSuccessors(fw_cand.twin(), 2, true);

            for (i = 0, count_succ = 0; i != 4; ++i) count_succ += !v_um[i].isEmpty;

            if (count_succ == 1) {

                b = fw_cand;
                return true;
            }
        }
    }

    return false;
}

template<typename T>
void CompactedDBG<T>::check_fp_tips(KmerHashTable<bool>& ignored_km_tips){

    uint64_t nb_real_short_tips = 0;

    size_t nxt_pos_insert_v_unitigs = v_unitigs.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t v_kmers_sz = v_kmers.size();

    vector<pair<int,int>> sp;

    for (KmerHashTable<bool>::iterator it = ignored_km_tips.begin(); it != ignored_km_tips.end(); ++it) {

        Kmer km = it.getKey();

        UnitigMap<T> cm = find(km, true); // Check if the (short) tip actually exists

        if (!cm.isEmpty){ // IF the tip exists

            nb_real_short_tips++;

            bool not_found = true;

            for (size_t i = 0; (i < 4) && not_found; ++i) {

                UnitigMap<T> cm_bw = find(km.backwardBase(alpha[i]));

                if (!cm_bw.isEmpty && !cm_bw.isAbundant && !cm_bw.isShort){

                    if (cm_bw.strand) cm_bw.dist++;

                    if ((cm_bw.dist != 0) && (cm_bw.dist != cm_bw.size - k_ + 1)){

                        sp.push_back(make_pair(0, cm_bw.dist));
                        sp.push_back(make_pair(cm_bw.dist, cm_bw.size - k_ + 1));

                        splitUnitig(cm_bw.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

                        sp.clear();
                    }

                    not_found = false;
                }
            }

            for (size_t i = 0; (i < 4) && not_found; ++i) {

                UnitigMap<T> cm_fw = find(km.forwardBase(alpha[i]));

                if (!cm_fw.isEmpty && !cm_fw.isAbundant && !cm_fw.isShort){

                    if (!cm_fw.strand) cm_fw.dist++;

                    if ((cm_fw.dist != 0) && (cm_fw.dist != cm_fw.size - k_ + 1)){

                        sp.push_back(make_pair(0, cm_fw.dist));
                        sp.push_back(make_pair(cm_fw.dist, cm_fw.size - k_ + 1));

                        splitUnitig(cm_fw.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

                        sp.clear();
                    }

                    not_found = false;
                }
            }
        }
    }

    if (nxt_pos_insert_v_unitigs < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert_v_unitigs);
    if (v_kmers_sz < v_kmers.size()) v_kmers.resize(v_kmers_sz);
}

template<typename T>
size_t CompactedDBG<T>::removeUnitigs(bool rmIsolated, bool clipTips, vector<Kmer>& v){

    if (!rmIsolated && !clipTips) return 0;

    const bool rm_and_clip = rmIsolated && clipTips;

    size_t v_unitigs_sz = v_unitigs.size();
    size_t v_kmers_sz = v_kmers.size();
    size_t removed = 0;
    size_t i;

    int64_t j;

    const int lim = (clipTips ? 1 : 0);

    int nb_pred, nb_succ;

    Kmer km;

    Unitig<T>* unitig = nullptr;

    for (j = 0; j < v_unitigs_sz; j++) {

        unitig = v_unitigs[j];

        if (unitig->numKmers() < k_){

            Kmer head = unitig->seq.getKmer(0);

            nb_pred = 0;

            for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

                if (!find(head.backwardBase(alpha[i]), true).isEmpty){

                    nb_pred++;
                    if (clipTips) km = head.backwardBase(alpha[i]);
                }
            }

            if (nb_pred <= lim){

                Kmer tail = unitig->seq.getKmer(unitig->seq.size() - k_);

                nb_succ = 0;

                for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                    if (!find(tail.forwardBase(alpha[i]), true).isEmpty){

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

        const pair<Kmer, CompressedCoverage_t<T>>& p = v_kmers[j];

        nb_pred = 0;

        for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

            if (!find(p.first.backwardBase(alpha[i]), true).isEmpty){

                nb_pred++;
                if (clipTips) km = p.first.backwardBase(alpha[i]);
            }
        }

        if (nb_pred <= lim){

            nb_succ = 0;

            for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                if (!find(p.first.forwardBase(alpha[i]), true).isEmpty){

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

    for (typename h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); ++it) {

        nb_pred = 0;

        for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

            if (!find(it.getKey().backwardBase(alpha[i]), true).isEmpty){

                ++nb_pred;
                if (clipTips) km = it.getKey().backwardBase(alpha[i]);
            }
        }

        if (nb_pred <= lim){

            nb_succ = 0;

            for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                if (!find(it.getKey().forwardBase(alpha[i]), true).isEmpty){

                    ++nb_succ;
                    if (clipTips) km = it.getKey().forwardBase(alpha[i]);
                }
            }

            if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))){

                ++removed;

                *it = CompressedCoverage_t<T>();

                if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
            }
        }
    }

    for (j = v_unitigs_sz; j < v_unitigs.size(); j++) deleteUnitig(false, false, j);
    v_unitigs.resize(v_unitigs_sz);

    for (j = v_kmers_sz; j < v_kmers.size(); j++) deleteUnitig(true, false, j);
    v_kmers.resize(v_kmers_sz);

    for (typename h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); ++it){

        if (it->ccov.size() == 0) deleteUnitig(false, true, it.getHash());
    }

    return removed;
}

template<typename T>
void CompactedDBG<T>::writeFASTA(string graphfilename) {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = v_kmers.size();
    const size_t graph_sz = size();

    size_t i = 0;

    ofstream graphfile;
    ostream graph(0);

    graphfile.open(graphfilename.c_str());
    graph.rdbuf(graphfile.rdbuf());
    assert(!graphfile.fail());

    for (size_t j = 0; j < v_unitigs_sz; ++j, ++i) {

        graph << ">" << i << "\n" << v_unitigs[j]->seq.toString() << (i == graph_sz - 1 ? "\0" : "\n");
    }

    for (size_t j = 0; j < v_kmers_sz; ++j, ++i) {

        graph << ">" << i << "\n" << v_kmers[j].first.toString() << (i == graph_sz - 1 ? "\0" : "\n");
    }

    for (typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); ++it, ++i) {

        graph << ">" << i << "\n" << it.getKey().toString() << (i == graph_sz - 1 ? "\0" : "\n");
    }

    graphfile.close();
}

template<typename T>
void CompactedDBG<T>::writeGFA(string graphfilename, const size_t nb_threads) {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = v_kmers.size();

    size_t i, labelA, labelB, id = v_unitigs_sz + v_kmers_sz + 1;

    UnitigMap<T> cand;

    Unitig<T>* unitig = nullptr;

    ofstream graphfile;
    ostream graph(0);

    graphfile.open(graphfilename.c_str());
    graph.rdbuf(graphfile.rdbuf());
    assert(!graphfile.fail());

    KmerHashTable<size_t> idmap(h_kmers_ccov.size());

    // gfa header
    graph << "H\tVN:Z:1.0\n";

    for (labelA = 1; labelA <= v_unitigs_sz; ++labelA) {

        unitig = v_unitigs[labelA - 1];

        graph << "S\t" << labelA << "\t" << unitig->seq.toString() << "\tLN:i:" <<
        unitig->seq.size() << "\tXC:i:" << unitig->coveragesum << "\n";
    }

    for (labelA = 1; labelA <= v_kmers_sz; ++labelA) {

        const pair<Kmer, CompressedCoverage_t<T>>& p = v_kmers[labelA - 1];

        graph << "S\t" << (labelA + v_unitigs_sz) << "\t" << p.first.toString() <<
        "\tLN:i:" << k_ << "\tXC:i:" << p.second.ccov.covAt(0) << "\n";
    }

    for (typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); ++it) {

        ++id;

        idmap.insert(it.getKey(), id);

        graph << "S\t" << id << "\t" << it.getKey().toString() << "\tLN:i:" << k_ <<
        "\tXC:i:" << it->ccov.covAt(0) << "\n";
    }

    if (nb_threads == 1){

        for (labelA = 1; labelA <= v_unitigs_sz; labelA++) {

            unitig = v_unitigs[labelA - 1];

            Kmer head = unitig->seq.getKmer(0);

            for (i = 0; i < 4; ++i) {

                Kmer b = head.backwardBase(alpha[i]);
                cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    graph << "L\t" << labelA << "\t-\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k_-1) << "M\n";
                }
            }

            Kmer tail = unitig->seq.getKmer(unitig->seq.size() - k_);

            for (i = 0; i < 4; ++i) {

                Kmer b = tail.forwardBase(alpha[i]);
                cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    graph << "L\t" << labelA << "\t+\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k_ - 1) << "M\n";
                }
            }
        }

        for (labelA = v_unitigs_sz + 1; labelA <= v_kmers_sz + v_unitigs_sz; labelA++) {

            const pair<Kmer, CompressedCoverage_t<T>>& p = v_kmers[labelA - v_unitigs_sz - 1];

            for (i = 0; i < 4; ++i) {

                Kmer b = p.first.backwardBase(alpha[i]);
                cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz : 0);

                    graph << "L\t" << labelA << "\t-\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k_ - 1) << "M\n";
                }
            }

            for (i = 0; i < 4; ++i) {

                Kmer b = p.first.forwardBase(alpha[i]);
                cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    graph << "L\t" << labelA << "\t+\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k_ - 1) << "M\n";
                }
            }
        }
    }
    else {

        auto worker_v_unitigs = [v_unitigs_sz, &idmap, this](const size_t labelA_start, const size_t labelA_end,
                                                             vector<pair<pair<size_t, bool>, pair<size_t, bool>>>* v_out){

            // We need to deal with the tail of long unitigs
            for (size_t labelA = labelA_start; labelA < labelA_end; ++labelA) {

                const Unitig<T>* unitig = v_unitigs[labelA - 1];

                const Kmer head = unitig->seq.getKmer(0);
                const Kmer tail = unitig->seq.getKmer(unitig->seq.size() - k_);

                vector<UnitigMap<T>> v_um = findPredecessors(head, true);

                for (const auto& um : v_um) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ? *(idmap.find(um.getHead().rep())) : um.pos_unitig + 1 + (um.isShort ? v_unitigs_sz : 0));
                        v_out->push_back(make_pair(make_pair(labelA, false), make_pair(labelB, um.strand)));
                    }
                }

                v_um.clear();

                v_um = findSuccessors(tail, 4, true);

                for (const auto& um : v_um) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ? *(idmap.find(um.getHead().rep())) : um.pos_unitig + 1 + (um.isShort ? v_unitigs_sz : 0));
                        v_out->push_back(make_pair(make_pair(labelA, true), make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        auto worker_v_kmers = [v_unitigs_sz, &idmap, this](const size_t labelA_start, const size_t labelA_end,
                                                             vector<pair<pair<size_t, bool>, pair<size_t, bool>>>* v_out){

            // We need to deal with the tail of long unitigs
            for (size_t labelA = labelA_start; labelA < labelA_end; ++labelA) {

                const pair<Kmer, CompressedCoverage_t<T>>& p = v_kmers[labelA - v_unitigs_sz - 1];

                vector<UnitigMap<T>> v_um = findPredecessors(p.first, true);

                for (const auto& um : v_um) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ? *(idmap.find(um.getHead().rep())) : um.pos_unitig + 1 + (um.isShort ? v_unitigs_sz : 0));
                        v_out->push_back(make_pair(make_pair(labelA, false), make_pair(labelB, um.strand)));
                    }
                }

                v_um.clear();

                v_um = findSuccessors(p.first, 4, true);

                for (const auto& um : v_um) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ? *(idmap.find(um.getHead().rep())) : um.pos_unitig + 1 + (um.isShort ? v_unitigs_sz : 0));
                        v_out->push_back(make_pair(make_pair(labelA, true), make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        const int chunk_size = 1000;

        vector<vector<pair<pair<size_t, bool>, pair<size_t, bool>>>> v_out(nb_threads);

        labelA = 1;

        while (labelA <= v_unitigs_sz) {

            vector<thread> workers;

            for (size_t i = 0; i != nb_threads; ++i) {

                if (labelA + chunk_size <= v_unitigs_sz){

                    workers.push_back(thread(worker_v_unitigs, labelA, labelA + chunk_size, &v_out[i]));

                    labelA += chunk_size;
                }
                else {

                    workers.push_back(thread(worker_v_unitigs, labelA, v_unitigs_sz + 1, &v_out[i]));

                    labelA += chunk_size;

                    break;
                }
            }

            for (auto &t : workers) t.join();

            for (size_t i = 0; i < nb_threads; ++i) {

                for (const auto& p : v_out[i]){

                    graph << "L\t" << p.first.first << "\t" << (p.first.second ? "+" : "-") << "\t" <<
                    p.second.first << "\t" << (p.second.second ? "+" : "-") << "\t" << (k_-1) << "M\n";
                }

                v_out[i].clear();
            }
        }

        labelA = v_unitigs_sz + 1;

        while (labelA <= v_kmers_sz + v_unitigs_sz) {

            vector<thread> workers;

            for (size_t i = 0; i != nb_threads; ++i) {

                if (labelA + chunk_size <= v_kmers_sz + v_unitigs_sz){

                    workers.push_back(thread(worker_v_kmers, labelA, labelA + chunk_size, &v_out[i]));

                    labelA += chunk_size;
                }
                else {

                    workers.push_back(thread(worker_v_kmers, labelA, v_kmers_sz + v_unitigs_sz + 1, &v_out[i]));

                    labelA += chunk_size;

                    break;
                }
            }

            for (auto &t : workers) t.join();

            for (size_t i = 0; i < nb_threads; ++i) {

                for (const auto& p : v_out[i]){

                    graph << "L\t" << p.first.first << "\t" << (p.first.second ? "+" : "-") << "\t" <<
                    p.second.first << "\t" << (p.second.second ? "+" : "-") << "\t" << (k_-1) << "M\n";
                }

                v_out[i].clear();
            }
        }
    }

    for (KmerHashTable<size_t>::iterator it = idmap.begin(); it != idmap.end(); it++) {

        labelA = *it;

        for (i = 0; i < 4; ++i) {

            Kmer b = it.getKey().backwardBase(alpha[i]);
            cand = find(b, true);

            if (!cand.isEmpty) {

                if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                graph << "L\t" << labelA << "\t-\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k_ - 1) << "M\n";
            }
        }

        for (i = 0; i < 4; ++i) {

            Kmer b = it.getKey().forwardBase(alpha[i]);
            cand = find(b, true);

            if (!cand.isEmpty) {

                if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                graph << "L\t" << labelA << "\t+\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k_ - 1) << "M\n";
            }
        }
    }

    graphfile.close();
}

template<typename T>
void CompactedDBG<T>::mapRead(const UnitigMap<T>& cc) {

    if (cc.isEmpty) return; // nothing maps, move on

    if (cc.isShort) v_kmers[cc.pos_unitig].second.ccov.cover(cc.dist, cc.dist + cc.len - 1);
    else if (cc.isAbundant) h_kmers_ccov.find(cc.pos_unitig)->ccov.cover(cc.dist, cc.dist + cc.len - 1);
    else v_unitigs[cc.pos_unitig]->ccov.cover(cc.dist, cc.dist + cc.len - 1);
}

template<typename T>
vector<Kmer> CompactedDBG<T>::extractMercyKmers(BlockedBloomFilter& bf_uniq_km, const size_t nb_threads, const bool verbose) {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = v_kmers.size();

    size_t i, j;

    char km_tmp[k_ + 1];

    KmerHashTable<uint8_t> tips;

    vector<Kmer> v_out;

    for (typename h_kmers_ccov_t::iterator it_ccov = h_kmers_ccov.begin(); it_ccov != h_kmers_ccov.end(); ++it_ccov) {

        const Kmer km = it_ccov.getKey().rep();

        vector<UnitigMap<T>> v_um = findPredecessors(km, true);

        for (i = 0; (i != 4) && v_um[i].isEmpty; ++i){}

        v_um = findSuccessors(km, true);

        for (j = 0; (j != 4) && v_um[j].isEmpty; ++j){}

        if ((i == 4) && (j == 4)) tips.insert(km, 3);
        else if (j == 4) tips.insert(km, 2);
        else if (i == 4) tips.insert(km, 1);
    }

    for (size_t it_v_km = 0; it_v_km != v_kmers_sz; ++it_v_km) {

        const Kmer km = v_kmers[it_v_km].first.rep();

        vector<UnitigMap<T>> v_um = findPredecessors(km, true);

        for (i = 0; (i != 4) && v_um[i].isEmpty; ++i){}

        v_um = findSuccessors(km, true);

        for (j = 0; (j != 4) && v_um[j].isEmpty; ++j){}

        if ((i == 4) && (j == 4)) tips.insert(km, 3);
        else if (j == 4) tips.insert(km, 2);
        else if (i == 4) tips.insert(km, 1);
    }

    for (size_t it_v_unitig = 0; it_v_unitig != v_unitigs_sz; ++it_v_unitig) {

        const CompressedSequence& seq = v_unitigs[it_v_unitig]->seq;

        const Kmer head = seq.getKmer(0);
        const Kmer tail = seq.getKmer(seq.size() - k_);

        vector<UnitigMap<T>> v_um = findPredecessors(head, true);

        for (i = 0; (i != 4) && v_um[i].isEmpty; ++i){}

        if (i == 4){

            const Kmer head_rep = head.rep();

            if (head == head_rep) tips.insert(head, 1);
            else tips.insert(head_rep, 2);
        }

        v_um = findSuccessors(tail, true);

        for (i = 0; (i != 4) && v_um[i].isEmpty; ++i){}

        if (i == 4){

            const Kmer tail_rep = tail.rep();

            if (tail == tail_rep) tips.insert(tail, 2);
            else tips.insert(tail_rep, 1);
        }
    }

    for (typename KmerHashTable<uint8_t>::iterator it_a = tips.begin(); it_a != tips.end(); ++it_a) {

        if ((*it_a == 1) || (*it_a == 3)){ // Corresponding k-mer has no predecessor in the graph

            const Kmer km_a = it_a.getKey();

            km_a.toString(km_tmp); // Get the k-mer

            RepHash rep_h(k_ - 1), rep_h_cpy; // Prepare its hash
            rep_h.init(km_tmp + 1);

            bool pres_neigh_bw[4] = {false, false, false, false};
            uint64_t hashes_bw[4];

            for (i = 0; i != 4; ++i) {

                rep_h_cpy = rep_h;
                rep_h_cpy.extendBW(alpha[i]);

                hashes_bw[i] = rep_h_cpy.hash(); // Prepare the hash of its predecessor
            }

            // Query the MBBF for all possible predecessors
            bf_uniq_km.contains(hashes_bw, minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash(), pres_neigh_bw, 4);

            for (i = 0; i != 4; ++i) {

                if (pres_neigh_bw[i]){ // A potential mercy k-mer has been found

                    const Kmer km_mercy = km_a.backwardBase(alpha[i]);

                    for (j = 0; j != 4; ++j) {

                        const Kmer km_b = km_mercy.backwardBase(alpha[j]); // Possible predecessor of this mercy k-mer
                        const Kmer km_b_rep = km_b.rep();

                        typename KmerHashTable<uint8_t>::iterator it_b = tips.find(km_b_rep); // Query the tips with predecessor

                        if (it_b != tips.end()){ // Possible predecessor of this mercy k-mer exists

                            if ((km_b == km_b_rep) && ((*it_b == 2) || (*it_b == 3))){

                                // Those tips can't be use anymore in these directions
                                *it_b = (*it_b == 2 ? 0 : 1);
                                *it_a = (*it_a == 1 ? 0 : 2);
                                i = 3; j = 3; // Break the loops

                                v_out.push_back(km_mercy); //Push the mercy k-mer found
                            }
                            else if ((km_b != km_b_rep) && ((*it_b == 1) || (*it_b == 3))){

                                // Those tips can't be use anymore in these directions
                                *it_b = (*it_b == 1 ? 0 : 2);
                                *it_a = (*it_a == 1 ? 0 : 2);
                                i = 3; j = 3; // Break the loops

                                v_out.push_back(km_mercy);  //Push the mercy k-mer found
                            }
                        }
                    }
                }
            }
        }

        if ((*it_a == 2) || (*it_a == 3)){ // Corresponding k-mer has no successor in the graph

            const Kmer km_a = it_a.getKey();

            km_a.toString(km_tmp); // Get the k-mer

            RepHash rep_h(k_ - 1), rep_h_cpy; // Prepare its hash
            rep_h.init(km_tmp + 1);

            bool pres_neigh_fw[4] = {false, false, false, false};
            uint64_t hashes_fw[4];

            for (i = 0; i != 4; ++i) {

                rep_h_cpy = rep_h;
                rep_h_cpy.extendFW(alpha[i]);

                hashes_fw[i] = rep_h_cpy.hash(); // Prepare the hash of its successor
            }

            // Query the MBBF for all possible predecessors
            bf_uniq_km.contains(hashes_fw, minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash(), pres_neigh_fw, 4);

            for (i = 0; i != 4; ++i) {

                if (pres_neigh_fw[i]){ // A potential mercy k-mer has been found

                    const Kmer km_mercy = km_a.forwardBase(alpha[i]);

                    for (j = 0; j != 4; ++j) {

                        const Kmer km_b = km_mercy.forwardBase(alpha[j]); // Possible successor of this mercy k-mer
                        const Kmer km_b_rep = km_b.rep();

                        typename KmerHashTable<uint8_t>::iterator it_b = tips.find(km_b_rep); // Query the tips with successor

                        if (it_b != tips.end()){ // Possible successor of this mercy k-mer exists

                            if ((km_b == km_b_rep) && ((*it_b == 1) || (*it_b == 3))){

                                // Those tips can't be use anymore in these directions
                                *it_b = (*it_b == 1 ? 0 : 2);
                                *it_a = (*it_a == 2 ? 0 : 1);
                                i = 3; j = 3; // Break the loops

                                v_out.push_back(km_mercy); //Push the mercy k-mer found
                            }
                            else if ((km_b != km_b_rep) && ((*it_b == 2) || (*it_b == 3))){

                                // Those tips can't be use anymore in these directions
                                *it_b = (*it_b == 2 ? 0 : 1);
                                *it_a = (*it_a == 2 ? 0 : 1);
                                i = 3; j = 3; // Break the loops

                                v_out.push_back(km_mercy);  //Push the mercy k-mer found
                            }
                        }
                    }
                }
            }
        }
    }

    if (verbose) cout << "CompactedDBG::extractMercyKmers(): " << v_out.size() << " k-mers extracted" << endl;

    return v_out;
}

template<typename T>
size_t CompactedDBG<T>::joinTips(string filename_MBBF_uniq_kmers, const size_t nb_threads, const bool verbose) {

    if (invalid){

        cerr << "CompactedDBG::joinTips(): Graph is invalid and tips cannot be joined" << endl;
        return 0;
    }

    BlockedBloomFilter mbbf;

    FILE* f_mbbf;

    if ((f_mbbf = fopen(filename_MBBF_uniq_kmers.c_str(), "rb")) == NULL){

        cerr << "CompactedDBG::joinTips(): Minimizer Blocked Bloom filter file of unique k-mers cannot be opened" << endl;
        return 0;
    }

    mbbf.ReadBloomFilter(f_mbbf);
    fclose(f_mbbf);

    vector<Kmer> v_mercy_km = extractMercyKmers(mbbf, nb_threads, verbose);

    for (const auto& km_mercy : v_mercy_km) addUnitig(km_mercy.rep().toString(), v_kmers.size());

    size_t nb_join = joinUnitigs(&v_mercy_km, nb_threads);

    if (verbose) cout << "CompactedDBG<T>::joinTips(): " << nb_join << " unitigs have been joined using mercy k-mers" << endl;

    return nb_join;
}

template<typename T>
void CompactedDBG<T>::print() const {

    cout << "CompactedDBG::print(): v_unitigs.size() = " << v_unitigs.size() << endl;
    cout << "CompactedDBG::print(): v_kmers.size() = " << v_kmers.size() << endl;
    cout << "CompactedDBG::print(): h_kmers_ccov.size() = " << h_kmers_ccov.size() << endl;
    cout << "CompactedDBG::print(): hmap_min_unitigs.size() = " << hmap_min_unitigs.size() << endl;
}
