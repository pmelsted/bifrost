#ifndef COMPACTED_DBG_TCC
#define COMPACTED_DBG_TCC

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

template<typename U, typename G>
CompactedDBG<U, G>::CompactedDBG(const int kmer_length, const int minimizer_length) :   k_(kmer_length), g_(minimizer_length), invalid(false) {

    setKmerGmerLength(k_, g_);
}

template<typename U, typename G>
CompactedDBG<U, G>::CompactedDBG(const CompactedDBG& o) :   k_(o.k_), g_(o.g_), invalid(o.invalid),
                                                            bf(o.bf), v_kmers(o.v_kmers), v_unitigs(o.v_unitigs.size()),
                                                            data(o.data), h_kmers_ccov(o.h_kmers_ccov),
                                                            hmap_min_unitigs(o.hmap_min_unitigs){

    for (size_t i = 0; i < o.v_unitigs.size(); ++i){

        v_unitigs[i] = new Unitig<U>;
        *(v_unitigs[i]) = *(o.v_unitigs[i]);
    }
}

template<typename U, typename G>
CompactedDBG<U, G>::CompactedDBG(CompactedDBG&& o) :    k_(o.k_), g_(o.g_), invalid(o.invalid),
                                                        bf(move(o.bf)), v_kmers(move(o.v_kmers)), data(move(o.data)),
                                                        v_unitigs(move(o.v_unitigs)), h_kmers_ccov(move(o.h_kmers_ccov)),
                                                        hmap_min_unitigs(move(o.hmap_min_unitigs)){

    o.clear();
}

template<typename U, typename G>
CompactedDBG<U, G>::~CompactedDBG() {

    clear();
}

template<typename U, typename G>
CompactedDBG<U, G>& CompactedDBG<U, G>::operator=(const CompactedDBG& o){

    clear();

    k_ = o.k_;
    g_ = o.g_;

    invalid = o.invalid;

    v_kmers = o.v_kmers;

    h_kmers_ccov = o.h_kmers_ccov;
    hmap_min_unitigs = o.hmap_min_unitigs;

    bf = o.bf;

    data = o.data;

    v_unitigs.reserve(o.v_unitigs.size());

    for (size_t i = 0; i < o.v_unitigs.size(); ++i){

        v_unitigs[i] = new Unitig<U>;
        *(v_unitigs[i]) = *(o.v_unitigs[i]);
    }

    return *this;
}

template<typename U, typename G>
CompactedDBG<U, G>& CompactedDBG<U, G>::operator=(CompactedDBG&& o){

    if (this != &o) {

        clear();

        k_ = o.k_;
        g_ = o.g_;

        invalid = o.invalid;

        v_kmers = move(o.v_kmers);
        v_unitigs = move(o.v_unitigs);

        h_kmers_ccov = move(o.h_kmers_ccov);
        hmap_min_unitigs = move(o.hmap_min_unitigs);

        bf = move(o.bf);

        data = move(o.data);

        o.clear();
    }

    return *this;
}

template<typename U, typename G>
CompactedDBG<U, G>& CompactedDBG<U, G>::operator+=(const CompactedDBG& o){

    merge(o, false);

    return *this;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::operator==(const CompactedDBG& o) const {

    if ((!invalid && !o.invalid) && (k_ == o.k_) && (size() == o.size())){

        for (const auto& unitig : *this){

            const_UnitigMap<U, G> unitig_o(o.find(unitig.getUnitigHead(), true));

            if (unitig_o.isEmpty) return false;
            else {

                unitig_o.dist = 0;
                unitig_o.len = unitig_o.size - k_ + 1;

                const string unitig_o_str = unitig_o.strand ? unitig_o.referenceUnitigToString() : reverse_complement(unitig_o.referenceUnitigToString());

                if (unitig_o_str != unitig.referenceUnitigToString()) return false;
            }
        }

        return true;
    }

    return false;
}

template<typename U, typename G>
inline bool CompactedDBG<U, G>::operator!=(const CompactedDBG& o) const {

    return !operator==(o);
}

template<typename U, typename G>
void CompactedDBG<U, G>::clear(){

    k_ = 0;
    g_ = 0;

    invalid = true;

    for (auto unitig : v_unitigs) delete unitig;

    v_unitigs.clear();
    v_kmers.clear();

    hmap_min_unitigs.clear();
    h_kmers_ccov.clear();

    bf.clear();
}

template<typename U, typename G>
bool CompactedDBG<U, G>::build(CDBG_Build_opt& opt){

    size_t max_threads = std::thread::hardware_concurrency();

    bool construct_finished = true;

    const int k_cpy = k_;
    const int g_cpy = g_;
    const bool invalid_cpy = invalid;

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

    if (opt.nb_non_unique_kmers > opt.nb_unique_kmers){

        cerr << "CompactedDBG::build(): The estimated number of non unique k-mers ";
        cerr << "cannot be greater than the estimated number of unique k-mers" << endl;
        construct_finished = false;
    }

    if (opt.read_chunksize <= 0){

        cerr << "CompactedDBG::build(): Chunk size of reads to share among threads cannot be less than or equal to 0" << endl;
        construct_finished = false;
    }

    if (opt.outFilenameBBF.length() != 0){

        FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

        if (fp == NULL) {

            cerr << "CompactedDBG::build(): Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
            construct_finished = false;
        }
        else {

            fclose(fp);

            if (std::remove(opt.outFilenameBBF.c_str()) != 0){

                cerr << "CompactedDBG::build(): Could not remove temporary file " << opt.outFilenameBBF << "." << endl;
            }
        }
    }

    if (opt.inFilenameBBF.length() != 0){

        if (check_file_exists(opt.inFilenameBBF)){

            FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

            if (fp == NULL) {

                cerr << "CompactedDBG::build(): Could not read input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                construct_finished = false;
            }
            else fclose(fp);
        }
        else {

            cerr << "CompactedDBG::build(): Input Blocked Bloom filter file " << opt.inFilenameBBF << " does not exist." << endl;
            construct_finished = false;
        }
    }

    if (opt.filename_seq_in.size() + opt.filename_ref_in.size() == 0){

        cerr << "CompactedDBG::build(): Number of FASTA/FASTQ files in input cannot be 0." << endl;
        construct_finished = false;
    }
    else {

        for (const auto& file : opt.filename_seq_in){

            if (check_file_exists(file)){

                FILE* fp = fopen(file.c_str(), "r");

                if (fp == NULL) {

                    cerr << "CompactedDBG::build(): Could not open input FASTA/FASTQ file " << file << endl;
                    construct_finished = false;
                }
                else fclose(fp);
            }
            else {

                cerr << "CompactedDBG::build(): Input file " << opt.inFilenameBBF << " does not exist." << endl;
                construct_finished = false;
            }
        }

        for (const auto& file : opt.filename_ref_in){

            if (check_file_exists(file)){

                FILE* fp = fopen(file.c_str(), "r");

                if (fp == NULL) {

                    cerr << "CompactedDBG::build(): Could not open input FASTA/FASTQ file " << file << endl;
                    construct_finished = false;
                }
                else fclose(fp);
            }
            else {

                cerr << "CompactedDBG::build(): Input file " << opt.inFilenameBBF << " does not exist." << endl;
                construct_finished = false;
            }
        }
    }

    clear();

    k_ = k_cpy;
    g_ = g_cpy;
    invalid = invalid_cpy;

    if (construct_finished){

        if ((opt.filename_seq_in.size() != 0) && (opt.filename_ref_in.size() != 0)){

            CompactedDBG<U, G> graph_seq(k_, g_);
            CompactedDBG<U, G> graph_ref(k_, g_);

            CDBG_Build_opt opt_seq(opt);
            CDBG_Build_opt opt_ref(opt);

            opt_seq.filename_ref_in.clear();
            opt_ref.filename_seq_in.clear();

            construct_finished = graph_seq.build(opt_seq);

            if (construct_finished) construct_finished = graph_ref.build(opt_ref);

            if (construct_finished){

                if (graph_ref.length() < graph_seq.length()){

                    construct_finished = graph_seq.merge(move(graph_ref), opt.nb_threads, opt.verbose);

                    if (construct_finished) *this = move(graph_seq);
                }
                else {

                    construct_finished = graph_ref.merge(move(graph_seq), opt.nb_threads, opt.verbose);

                    if (construct_finished) *this = move(graph_ref);
                }
            }

            CompressedCoverage::setFullCoverage(2);
        }
        else {

            const bool reference_mode = (opt.filename_ref_in.size() != 0);

            const vector<string>& v_files = reference_mode ? opt.filename_ref_in : opt.filename_seq_in;

            if (reference_mode) CompressedCoverage::setFullCoverage(1);
            else CompressedCoverage::setFullCoverage(2);

            if (opt.inFilenameBBF.length() != 0){

                FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "CompactedDBG::build(): Could not open input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                    construct_finished = false;
                }
                else {

                    construct_finished = bf.ReadBloomFilter(fp);

                    fclose(fp);
                }
            }
            else {

                if ((opt.nb_unique_kmers == 0) || (!reference_mode && (opt.nb_non_unique_kmers == 0))){

                    KmerStream_Build_opt kms_opt;

                    kms_opt.threads = opt.nb_threads;
                    kms_opt.verbose = opt.verbose;
                    kms_opt.k = opt.k;
                    kms_opt.q = 0;

                    for (const auto& s : v_files) kms_opt.files.push_back(s);

                    KmerStream kms(kms_opt);

                    opt.nb_unique_kmers = max(1UL, kms.F0());

                    if (!reference_mode) opt.nb_non_unique_kmers = max(1UL, opt.nb_unique_kmers - min(opt.nb_unique_kmers, kms.f1()));

                    if (opt.verbose){

                        cout << "CompactedDBG::build(): Estimated number of k-mers occurring at least once: " << opt.nb_unique_kmers << endl;

                        if (!reference_mode){

                            cout << "CompactedDBG::build(): Estimated number of k-mers occurring twice or more: " << opt.nb_non_unique_kmers << endl;
                        }
                    }
                }

                construct_finished = filter(opt);
            }

            if (construct_finished){

                if (opt.outFilenameBBF.length() != 0){

                    FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

                    if (fp == NULL) {

                        cerr << "CompactedDBG::build(): Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
                        construct_finished = false;
                    }
                    else {

                        bf.WriteBloomFilter(fp);

                        fclose(fp);
                    }
                }

                if (construct_finished) construct_finished = construct(opt); // Construction step

                bf.clear();
            }
        }
    }

    return construct_finished;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::simplify(const bool delete_short_isolated_unitigs, const bool clip_short_tips, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::simplify(): Graph is invalid and cannot be simplified" << endl;
        return false;
    }

    if (delete_short_isolated_unitigs || clip_short_tips){

        if (verbose) cout << endl << "CompactedDBG::simplify(): Removing isolated unitigs and/or clipping tips" << endl;

        vector<Kmer> v_joins;
        size_t joined = 0;

        size_t removed = removeUnitigs(delete_short_isolated_unitigs, clip_short_tips, v_joins);

        if (clip_short_tips) joined = joinUnitigs_<is_void<U>::value>(&v_joins);

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

template<typename U, typename G>
bool CompactedDBG<U, G>::write(const string& output_filename, const size_t nb_threads, const bool GFA_output, const bool verbose) const {

    if (invalid){

        cerr << "CompactedDBG::write(): Graph is invalid and cannot be written to disk" << endl;
        return false;
    }

    if (nb_threads > std::thread::hardware_concurrency()){

        cerr << "CompactedDBG::write(): Number of threads cannot exceed " << std::thread::hardware_concurrency() << "threads" << endl;
        return false;
    }

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

    GFA_output ? writeGFA(out, nb_threads) : writeFASTA(out);

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::read(const string& input_filename, const bool verbose){

    if (verbose) cout << endl << "CompactedDBG::read(): Reading graph from disk" << endl;

    string s_ext = input_filename.substr(input_filename.find_last_of(".") + 1);

    if ((s_ext == "gz")){

        s_ext = s_ext.substr(s_ext.find_last_of(".") + 1);

        if ((s_ext != "fasta") && (s_ext != "fa") && (s_ext != "gfa")) {

            cerr << "CompactedDBG::read(): Input files must be in FASTA (*.fasta, *.fa, *.fasta.gz, *.fa.gz) or " <<
            "GFA (*.gfa) format" << endl;
            cerr << "CompactedDBG::read(): Erroneous file is " << input_filename << endl;

            return false;
        }
    }
    else if ((s_ext != "fasta") && (s_ext != "fa") && (s_ext != "gfa")) {

        cerr << "CompactedDBG::read(): Input files must be in FASTA (*.fasta, *.fa, *.fasta.gz, *.fa.gz) or " <<
        "GFA (*.gfa) format" << endl;
        cerr << "CompactedDBG::read(): Erroneous file is " << input_filename << endl;

        return false;
    }

    if (s_ext == "gfa"){

        FILE* fp = fopen(input_filename.c_str(), "r");

        if (fp == NULL) {

            cerr << "CompactedDBG::read(): Could not open file " << input_filename << " for reading graph" << endl;
            return false;
        }

        fclose(fp);

        char buffer[4096];

        ifstream graphfile_in(input_filename);
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

        setKmerGmerLength(k, g);

        if (!invalid) readGFA(input_filename);

        return !invalid;
    }
    else {

        const int k = k_;
        const int g = g_;

        clear();

        setKmerGmerLength(k, g);

        readFASTA(input_filename);
    }

    CompressedCoverage::setFullCoverage(1);

    for (auto& unitig : *this) unitig.setFullCoverage();

    invalid = false;

    return true;
}

template<typename U, typename G>
const_UnitigMap<U, G> CompactedDBG<U, G>::find(const Kmer& km, const bool extremities_only) const {

    if (invalid){

        cerr << "CompactedDBG::find(): Graph is invalid and cannot be searched" << endl;
        return const_UnitigMap<U, G>();
    }

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    char km_tmp[MAX_KMER_SIZE];
    km.toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k_, k_, g_, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz(Minimizer(&km_tmp[min_h_res.pos]).rep());
        hmap_min_unitigs_t::const_iterator it = hmap_min_unitigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVal1();
            const uint8_t flag_v = it.getVal2();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()) return const_UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km == km_rep, this);
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

                        if (((min_h_res.pos == unitig_id_pos) || (min_h_res.pos == diff - unitig_id_pos)) && (v_kmers[unitig_id].first == km_rep)){

                            return const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km == km_rep, this);
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - min_h_res.pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km_twin)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km_twin)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return const_UnitigMap<U, G>();
}

template<typename U, typename G>
UnitigMap<U, G> CompactedDBG<U, G>::find(const Kmer& km, const bool extremities_only) {

    if (invalid){

        cerr << "CompactedDBG::find(): Graph is invalid and cannot be searched" << endl;
        return UnitigMap<U, G>();
    }

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    char km_tmp[MAX_KMER_SIZE];
    km.toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k_, k_, g_, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz(Minimizer(&km_tmp[min_h_res.pos]).rep());
        hmap_min_unitigs_t::const_iterator it = hmap_min_unitigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVal1();
            const uint8_t flag_v = it.getVal2();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()) return UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km == km_rep, this);
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

                        if (((min_h_res.pos == unitig_id_pos) || (min_h_res.pos == diff - unitig_id_pos)) && (v_kmers[unitig_id].first == km_rep)){

                            return UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km == km_rep, this);
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - min_h_res.pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km_twin)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km_twin)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return UnitigMap<U, G>();
}

template<typename U, typename G>
vector<const_UnitigMap<U, G>> CompactedDBG<U, G>::findPredecessors(const Kmer& km, const bool extremities_only) const {

    const Kmer km_pred[4] = {km.backwardBase('A'), km.backwardBase('C'), km.backwardBase('G'), km.backwardBase('T')};
    const Kmer km_rep[4] = {km_pred[0].rep(), km_pred[1].rep(), km_pred[2].rep(), km_pred[3].rep()};

    const Kmer km_twin_a = km_pred[0].twin();

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    char km_tmp[MAX_KMER_SIZE];
    km_pred[0].toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k_, k_, g_, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    vector<const_UnitigMap<U, G>> v_um(4, const_UnitigMap<U, G>(1, this));

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz(Minimizer(&km_tmp[min_h_res.pos]).rep());
        hmap_min_unitigs_t::const_iterator it = hmap_min_unitigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVal1();
            const uint8_t flag_v = it.getVal2();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km;

                        for (size_t j = 0; j != 4; ++j){

                            if ((it_km = h_kmers_ccov.find(km_rep[j])) != h_kmers_ccov.end()){

                                v_um[j].partialCopy(const_UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km_pred[j] == km_rep[j], this));
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

                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_pred[idx] == km_rep[idx], this));
                            }
                            else {

                                idx = 0x3 - bits[v_kmers[unitig_id].first.getChar(k_ - 1)];

                                if (v_kmers[unitig_id].first == km_rep[idx]) {

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_pred[idx] == km_rep[idx], this));
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
                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this));
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_ - 1, km_twin_a)){

                                const uint8_t idx = 0x3 - bits[v_unitigs[unitig_id]->seq.getChar(pos_match + k_ - 1)];
                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this));
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match + 1, k_ - 1, km)){

                                const uint8_t idx = bits[v_unitigs[unitig_id]->seq.getChar(pos_match)];
                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this));
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_ - 1, km_twin_a)){

                                const uint8_t idx = 0x3 - bits[v_unitigs[unitig_id]->seq.getChar(pos_match + k_ - 1)];
                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this));
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

template<typename U, typename G>
vector<const_UnitigMap<U, G>> CompactedDBG<U, G>::findSuccessors(const Kmer& km, const size_t limit, const bool extremities_only) const {

    vector<const_UnitigMap<U, G>> v_um(4, const_UnitigMap<U, G>(1, this));

    if (limit == 0) return v_um;

    const Kmer km_succ[4] = {km.forwardBase('A'), km.forwardBase('C'), km.forwardBase('G'), km.forwardBase('T')};
    const Kmer km_rep[4] = {km_succ[0].rep(), km_succ[1].rep(), km_succ[2].rep(), km_succ[3].rep()};

    const Kmer km_twin_a = km_succ[0].twin().forwardBase('A');

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    int nb_found = 0;

    const int diff = k_ - g_;

    char km_tmp[MAX_KMER_SIZE];
    km_succ[0].toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k_, k_, g_, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz(Minimizer(&km_tmp[min_h_res.pos]).rep());
        hmap_min_unitigs_t::const_iterator it = hmap_min_unitigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVal1();
            const uint8_t flag_v = it.getVal2();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km;

                        for (size_t j = 0; j != 4; ++j){

                            if (v_um[j].isEmpty && ((it_km = h_kmers_ccov.find(km_rep[j])) != h_kmers_ccov.end())){

                                v_um[j].partialCopy(const_UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km_succ[j] == km_rep[j], this));
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

                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_succ[idx] == km_rep[idx], this));
                                if (++nb_found == limit) return v_um;
                            }
                            else {

                                idx = 0x3 - bits[v_kmers[unitig_id].first.getChar(0)];

                                if (v_um[idx].isEmpty && (v_kmers[unitig_id].first == km_rep[idx])){

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_succ[idx] == km_rep[idx], this));
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

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this));
                                    if (++nb_found == limit) return v_um;
                                }
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->seq.compareKmer(pos_match + 1, k_ - 1, km_twin_a)){

                                const int idx = 0x3 - bits[v_unitigs[unitig_id]->seq.getChar(pos_match)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this));
                                    if (++nb_found == limit) return v_um;
                                }
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_ - 1, km_succ[0])){

                                const int idx = bits[v_unitigs[unitig_id]->seq.getChar(pos_match + k_ - 1)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this));
                                    if (++nb_found == limit) return v_um;
                                }

                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match + 1, k_ - 1, km_twin_a)){

                                const int idx = 0x3 - bits[v_unitigs[unitig_id]->seq.getChar(pos_match)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this));
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

template<typename U, typename G>
bool CompactedDBG<U, G>::add(const string& seq, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::add(): Graph is invalid and no sequence can be added to it" << endl;
        return false;
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::add(): Input sequence length cannot be less than k = " << k_ << endl;
        return false;
    }

    const char* str_seq = seq.c_str();

    string no_dup_km_seq;

    unordered_set<Kmer, KmerHash> km_seen;

    size_t last_start_pos = 0;

    for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end; ++it_km) { //non-ACGT char. are discarded

        const std::pair<Kmer, int>& p = *it_km;

        if (!km_seen.insert(p.first.rep()).second){

            no_dup_km_seq = seq.substr(last_start_pos, p.second - last_start_pos + k_ - 1);
            last_start_pos = p.second;

            std::transform(no_dup_km_seq.begin(), no_dup_km_seq.end(), no_dup_km_seq.begin(), ::toupper);

            mergeUnitig(no_dup_km_seq, verbose);
            km_seen.clear();
        }
    }

    no_dup_km_seq = seq.substr(last_start_pos, seq.length() - last_start_pos + k_ - 1);

    std::transform(no_dup_km_seq.begin(), no_dup_km_seq.end(), no_dup_km_seq.begin(), ::toupper);

    mergeUnitig(no_dup_km_seq, verbose);

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::remove(const const_UnitigMap<U, G>& um, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::remove(): Graph is invalid and no unitig can be removed from it" << endl;
        return false;
    }

    vector<Kmer> v_km;

    const Kmer head(um.getUnitigHead());
    const Kmer tail(um.getUnitigTail());

    for (size_t i = 0; i != 4; ++i){

        const Kmer bw(head.backwardBase(alpha[i]));

        if (!find(bw, true).isEmpty) v_km.push_back(bw);
    }

    for (size_t i = 0; i != 4; ++i){

        const Kmer fw(tail.forwardBase(alpha[i]));

        if (!find(fw, true).isEmpty) v_km.push_back(fw);
    }

    const size_t swap_position = (um.isShort ? v_kmers.size() : v_unitigs.size()) - 1;

    if (!um.isAbundant && (um.pos_unitig != swap_position)) swapUnitigs(um.isShort, um.pos_unitig, swap_position);

    deleteUnitig_<is_void<U>::value>(um.isShort, um.isAbundant, swap_position);

    if (um.isShort) v_kmers.resize(swap_position);
    else if (!um.isAbundant) v_unitigs.resize(swap_position);

    joinUnitigs_<is_void<U>::value>(&v_km);

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::merge(const CompactedDBG& o, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) cerr << "CompactedDBG::merge(): Current graph is invalid." << endl;
         ret = false;
    }

    if (o.invalid){

         if (verbose) cerr << "CompactedDBG::merge(): Graph to merge is invalid." << endl;
         ret = false;
    }

    if (k_ != o.k_){

         if (verbose) cerr << "CompactedDBG::merge(): The graphs to merge do not have the same k-mer length." << endl;
         ret = false;
    }

    if (g_ != o.g_){

         if (verbose) cerr << "CompactedDBG::merge(): The graphs to merge do not have the same minimizer length." << endl;
         ret = false;
    }

    if (this == &o){

         if (verbose) cerr << "CompactedDBG::merge(): Cannot merge graph with itself." << endl;
         ret = false;
    }

    if (ret){

        const size_t sz_before = size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        if (annotateSplitUnitigs(o, nb_threads, verbose)){

            const size_t sz_after = size();
            const pair<size_t, size_t> p = splitAllUnitigs();
            const size_t joined = (p.second != 0) ? joinUnitigs_<is_void<U>::value>() : 0;

            if (verbose){

                cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Split " << p.first << " unitigs into " << p.second << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << endl;
                cout << "CompactedDBG::merge(): " << size() << " unitigs after merging." << endl;
            }

            if (!is_void<U>::value) mergeData(o, nb_threads, verbose);

            return true;
        }
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::merge(const vector<CompactedDBG>& v, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) cerr << "CompactedDBG::merge(): Current graph is invalid." << endl;
         ret = false;
    }

    for (const auto& cdbg : v){

        if (cdbg.invalid){

             if (verbose) cerr << "CompactedDBG::merge(): One of the graph to merge is invalid." << endl;
             ret = false;
        }

        if (k_ != cdbg.k_){

             if (verbose) cerr << "CompactedDBG::merge(): The graphs to merge do not have the same k-mer length." << endl;
             ret = false;
        }

        if (g_ != cdbg.g_){

             if (verbose) cerr << "CompactedDBG::merge(): The graphs to merge do not have the same minimizer length." << endl;
             ret = false;
        }

        if (this == &cdbg){

             if (verbose) cerr << "CompactedDBG::merge(): Cannot merge graph with itself." << endl;
             ret = false;
        }
    }

    if (ret){

        const size_t sz_before = size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        for (const auto& cdbg : v){

            ret = annotateSplitUnitigs(cdbg, nb_threads, verbose);

            if (!ret) break;
        }

        if (ret){

            const size_t sz_after = size();
            const pair<size_t, size_t> p = splitAllUnitigs();
            const size_t joined = (p.second != 0) ? joinUnitigs_<is_void<U>::value>() : 0;

            if (verbose){

                cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Split " << p.first << " unitigs into " << p.second << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << endl;
                cout << "CompactedDBG::merge(): " << size() << " unitigs after merging." << endl;
            }

            if (!is_void<U>::value){

                for (const auto& cdbg : v) mergeData(cdbg, nb_threads, verbose);
            }

            return true;
        }
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::annotateSplitUnitigs(const CompactedDBG<U, G>& o, const size_t nb_threads, const bool verbose){

    if ((this != &o) && !invalid && !o.invalid) { // TODO: Check if k_ and g_ are the same

        if (verbose){

            cout << "CompactedDBG::annotateSplitUnitigs(): Current graph has " << size() << " unitigs." << endl;
            cout << "CompactedDBG::annotateSplitUnitigs(): Graph to merge has " << o.size() << " unitigs." << endl;
            cout << "CompactedDBG::annotateSplitUnitigs(): Start unitigs merging." << endl;
        }

        if (nb_threads == 1){

            for (const auto& unitig : o) annotateSplitUnitig(unitig.referenceUnitigToString(), false);
        }
        else {

            const size_t chunk = 100;

            vector<thread> workers; // need to keep track of threads so we can join them

            typename CompactedDBG<U, G>::const_iterator g_a(o.begin());
            typename CompactedDBG<U, G>::const_iterator g_b(o.end());

            LockGraph lck_g(nb_threads);

            mutex mutex_o_unitig;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        typename CompactedDBG<U, G>::const_iterator l_a, l_b;

                        while (true) {

                            {
                                unique_lock<mutex> lock(mutex_o_unitig);

                                if (g_a == g_b) return;

                                l_a = g_a;
                                l_b = g_a;

                                for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                                g_a = l_b;
                            }

                            for (auto& it_unitig = l_a; it_unitig != l_b; ++it_unitig) {

                                annotateSplitUnitig(it_unitig->referenceUnitigToString(), lck_g, t, false);
                            }
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }

        if (verbose) cout << "CompactedDBG::annotateSplitUnitigs(): Merging unitigs finished." << endl;

        return true;
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::mergeData(const CompactedDBG<U, G>& o, const size_t nb_threads, const bool verbose){

    if ((this != &o) && !invalid && !o.invalid){ // TODO: Check if k_ and g_ are the same

        if (verbose) cout << "CompactedDBG::mergeData(): Merging data started." << endl;

        const size_t nb_locks = nb_threads * 1024;

        std::atomic_flag* locks_unitig = new std::atomic_flag[nb_locks];

        for (size_t i = 0; i < nb_locks; ++i) locks_unitig[i].clear();

        auto worker_function = [&](typename CompactedDBG<U, G>::const_iterator it_a, typename CompactedDBG<U, G>::const_iterator it_b) {

            while (it_a != it_b) {

                const string str(it_a->referenceUnitigToString());
                const char* str_seq = str.c_str();

                for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end;) { //non-ACGT char. are discarded

                    const std::pair<Kmer, int>& p = *it_km;

                    const UnitigMap<U, G> um_dest = findUnitig(p.first, str_seq, p.second);

                    if (um_dest.isEmpty) ++it_km;
                    else {

                        const_UnitigMap<U, G> um_src(*it_a);

                        um_src.dist = p.second;
                        um_src.len = um_dest.len;

                        const uint64_t id_lock = um_dest.pos_unitig % nb_locks;

                        while (locks_unitig[id_lock].test_and_set(std::memory_order_acquire));

                        mergeData_<is_void<U>::value>(um_dest, um_src);

                        locks_unitig[id_lock].clear(std::memory_order_release);

                        it_km += um_dest.len;
                    }
                }

                ++it_a;
            }
        };

        const size_t chunk = 100;

        vector<thread> workers; // need to keep track of threads so we can join them

        typename CompactedDBG<U, G>::const_iterator g_a = o.begin();
        typename CompactedDBG<U, G>::const_iterator g_b = o.end();

        mutex mutex_it;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    typename CompactedDBG<U, G>::const_iterator l_a, l_b;

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

        if (verbose) cout << "CompactedDBG::mergeData(): Merging data finished." << endl;

        delete[] locks_unitig;

        return true;
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::mergeData(CompactedDBG<U, G>&& o, const size_t nb_threads, const bool verbose){

    if ((this != &o) && !invalid && !o.invalid){ // TODO: Check if k_ and g_ are the same

        if (verbose) cout << "CompactedDBG::mergeData(): Merging data started." << endl;

        const size_t nb_locks = nb_threads * 1024;

        std::atomic_flag* locks_unitig = new std::atomic_flag[nb_locks];

        for (size_t i = 0; i < nb_locks; ++i) locks_unitig[i].clear();

        auto worker_function = [&](typename CompactedDBG<U, G>::iterator it_a, typename CompactedDBG<U, G>::iterator it_b) {

            while (it_a != it_b) {

                const string str(it_a->referenceUnitigToString());
                const char* str_seq = str.c_str();

                for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end;) { //non-ACGT char. are discarded

                    const std::pair<Kmer, int>& p = *it_km;

                    const UnitigMap<U, G> um_dest = findUnitig(p.first, str_seq, p.second);

                    if (um_dest.isEmpty) ++it_km;
                    else {

                        const_UnitigMap<U, G> um_src(*it_a);

                        um_src.dist = p.second;
                        um_src.len = um_dest.len;

                        const uint64_t id_lock = um_dest.pos_unitig % nb_locks;

                        while (locks_unitig[id_lock].test_and_set(std::memory_order_acquire));

                        mergeData_<is_void<U>::value>(um_dest, um_src);

                        locks_unitig[id_lock].clear(std::memory_order_release);

                        it_km += um_dest.len;
                    }
                }

                it_a->getData()->clear(*it_a);

                ++it_a;
            }
        };

        const size_t chunk = 100;

        vector<thread> workers; // need to keep track of threads so we can join them

        typename CompactedDBG<U, G>::iterator g_a = o.begin();
        typename CompactedDBG<U, G>::iterator g_b = o.end();

        mutex mutex_it;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    typename CompactedDBG<U, G>::iterator l_a, l_b;

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

        if (verbose) cout << "CompactedDBG::mergeData(): Merging data finished." << endl;

        delete[] locks_unitig;

        return true;
    }

    return false;
}

template<typename U, typename G>
typename CompactedDBG<U, G>::iterator CompactedDBG<U, G>::begin() {

    if (invalid) return iterator();

    iterator it(this);
    ++it;

    return it;
}

template<typename U, typename G>
typename CompactedDBG<U, G>::const_iterator CompactedDBG<U, G>::begin() const {

    if (invalid) return const_iterator();

    const_iterator it(this);
    ++it;

    return it;
}

template<typename U, typename G>
typename CompactedDBG<U, G>::iterator CompactedDBG<U, G>::end() { return iterator(); }

template<typename U, typename G>
typename CompactedDBG<U, G>::const_iterator CompactedDBG<U, G>::end() const { return const_iterator(); }

template<typename U, typename G>
size_t CompactedDBG<U, G>::length() const {

    size_t len = 0;

    if (!invalid){

        len = (v_kmers.size() + h_kmers_ccov.size()) * k_;

        for (const auto& unitig : v_unitigs) len += unitig->seq.size();
    }

    return len;
}

template<typename U, typename G>
size_t CompactedDBG<U, G>::nbKmers() const {

    size_t nb = 0;

    if (!invalid){

        nb = (v_kmers.size() + h_kmers_ccov.size());

        for (const auto& unitig : v_unitigs) nb += unitig->seq.size() - k_ + 1;
    }

    return nb;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::filter(const CDBG_Build_opt& opt) {

    if (invalid){

        cerr << "CompactedDBG::filter(): Graph is invalid and it cannot be built" << endl;
        return false;
    }

    BlockedBloomFilter bf_tmp;

    const bool reference_mode = (opt.filename_ref_in.size() != 0);

    if (reference_mode){

        BlockedBloomFilter tmp(opt.nb_unique_kmers, opt.nb_bits_unique_kmers_bf);
        bf = move(tmp);
    }
    else {

        {
            BlockedBloomFilter tmp(opt.nb_unique_kmers, opt.nb_bits_unique_kmers_bf);
            bf_tmp = move(tmp);
        }

        {
            BlockedBloomFilter tmp(opt.nb_non_unique_kmers, opt.nb_bits_non_unique_kmers_bf);
            bf = move(tmp);
        }
    }

    string s;

    size_t len_read = 0;
    size_t pos_read = 0;
    size_t nb_seq = 0;

    const size_t max_len_seq = 1024;
    const size_t thread_seq_buf_sz = opt.read_chunksize * max_len_seq;

    const bool multi_threaded = (opt.nb_threads != 1);

    atomic<uint64_t> num_kmers(0), num_ins(0);

    const vector<string>& filename_in = reference_mode ? opt.filename_ref_in : opt.filename_seq_in;

    FileParser fp(filename_in);

    // Main worker thread
    auto worker_function = [&](char* seq_buf, const size_t seq_buf_sz) {

        uint64_t l_num_kmers = 0, l_num_ins = 0;

        char* str = seq_buf;
        const char* str_end = &seq_buf[seq_buf_sz];

        while (str < str_end) { // for each input

            const int len = strlen(str);

            for (char* s = str; s != &str[len]; ++s) *s &= 0xDF; // Put characters in upper case

            KmerHashIterator<RepHash> it_kmer_h(str, len, k_), it_kmer_h_end;
            minHashIterator<RepHash> it_min(str, len, k_, g_, RepHash(), true);

            if (reference_mode){

                for (int last_pos = -1; it_kmer_h != it_kmer_h_end; ++it_kmer_h, ++it_min, ++l_num_kmers) {

                    const std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                    // If one or more k-mer were jumped because contained non-ACGT char.
                    if (p_.second != last_pos + 1)
                        it_min = minHashIterator<RepHash>(&str[p_.second], len - p_.second, k_, g_, RepHash(), true);

                    last_pos = p_.second;

                    bf.insert(p_.first, it_min.getHash(), multi_threaded);

                    ++l_num_ins;
                }
            }
            else {

                for (int last_pos = -1; it_kmer_h != it_kmer_h_end; ++it_kmer_h, ++it_min, ++l_num_kmers) {

                    const std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                    // If one or more k-mer were jumped because contained non-ACGT char.
                    if (p_.second != last_pos + 1)
                        it_min = minHashIterator<RepHash>(&str[p_.second], len - p_.second, k_, g_, RepHash(), true);

                    last_pos = p_.second;
                    const uint64_t min_hr = it_min.getHash();

                    if (bf_tmp.insert(p_.first, min_hr, multi_threaded)) ++l_num_ins;
                    else bf.insert(p_.first, min_hr, multi_threaded);
                }
            }

            str += len + 1;
        }

        // atomic adds
        num_kmers += l_num_kmers;
        num_ins += l_num_ins;
    };

    auto reading_function = [&](char* seq_buf, size_t& seq_buf_sz) {

        size_t file_id = 0;

        const size_t sz_buf = thread_seq_buf_sz - k_;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_buf) {

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

                nb_seq += new_reading;

                pos_read = (new_reading ? 0 : pos_read);
                len_read = s.length();

                s_str = s.c_str();

                if (len_read >= k_){

                    if ((thread_seq_buf_sz - seq_buf_sz - 1) < (len_read - pos_read)){

                        strncpy(&seq_buf[seq_buf_sz], &s_str[pos_read], thread_seq_buf_sz - seq_buf_sz - 1);

                        seq_buf[thread_seq_buf_sz - 1] = '\0';

                        pos_read += sz_buf - seq_buf_sz;
                        seq_buf_sz = thread_seq_buf_sz;

                        break;
                    }
                    else {

                        strcpy(&seq_buf[seq_buf_sz], &s_str[pos_read]);

                        seq_buf_sz += (len_read - pos_read) + 1;
                        pos_read = len_read;
                    }
                }
                else pos_read = len_read;
            }
            else return true;
        }

        return false;
    };

    {
        bool stop = false;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_file;

        char* buffer_seq = new char[opt.nb_threads * thread_seq_buf_sz]();
        size_t* buffer_seq_sz = new size_t[opt.nb_threads]();

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_file);

                            if (stop) return;

                            stop = reading_function(&buffer_seq[t * thread_seq_buf_sz], buffer_seq_sz[t]);
                        }

                        worker_function(&buffer_seq[t * thread_seq_buf_sz], buffer_seq_sz[t]);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        delete[] buffer_seq;
        delete[] buffer_seq_sz;
    }

    fp.close();

    if (opt.verbose) {

        cout << "CompactedDBG::filter(): Closed all fasta/fastq files" << endl;
        cout << "CompactedDBG::filter(): Processed " << num_kmers << " k-mers in " << nb_seq  << " reads" << endl;
        cout << "CompactedDBG::filter(): Found " << num_ins << " unique k-mers" << endl;
        cout << "CompactedDBG::filter(): Number of blocks in Bloom filter is " << bf.getNbBlocks() << endl;
    }

    if (opt.useMercyKmers && !reference_mode){

        const string mbbf_uniq_filename(opt.prefixFilenameOut + "_uniq");

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

template<typename U, typename G>
bool CompactedDBG<U, G>::construct(const CDBG_Build_opt& opt){

    if (invalid){

        cerr << "CompactedDBG::construct(): Graph is invalid and cannot be built" << endl;
        return false;
    }

    const bool reference_mode = (opt.filename_ref_in.size() != 0);

    const vector<string>& filename_in = reference_mode ? opt.filename_ref_in : opt.filename_seq_in;

    FileParser fp(filename_in);

    string s;

    size_t len_read = 0;
    size_t pos_read = 0;

    const size_t max_len_seq = 1024;
    const size_t thread_seq_buf_sz = opt.read_chunksize * max_len_seq;

    tiny_vector<Kmer, 2>* fp_candidate = nullptr;

    KmerHashTable<bool> ignored_km_tips;

    const size_t nb_locks = opt.nb_threads * 1024;

    std::atomic_flag lock_ignored_km_tips = ATOMIC_FLAG_INIT;

    vector<std::atomic_flag> locks_fp;
    vector<std::atomic_flag> locks_mapping(nb_locks);
    vector<std::atomic_flag> locks_unitig(opt.nb_threads);

    for (auto& lck : locks_mapping) lck.clear();
    for (auto& lck : locks_unitig) lck.clear();

    if (!reference_mode){

        fp_candidate = new tiny_vector<Kmer, 2>[bf.getNbBlocks()];
        locks_fp = vector<std::atomic_flag>(nb_locks);

        for (auto& lck : locks_fp) lck.clear();
    }

    auto worker_function = [&](char* seq_buf, const size_t seq_buf_sz, const size_t thread_id) {

        vector<Kmer> l_ignored_km_tips;

        uint64_t it_min_h, last_it_min_h;

        BlockedBloomFilter::BBF_Blocks block_bf;

        Kmer km;
        RepHash rep;

        char* str = seq_buf;
        const char* str_end = &seq_buf[seq_buf_sz];

        while (str < str_end) { // for each input

            const int len = strlen(str);

            for (char* s = str; s != &str[len]; ++s) *s &= 0xDF;

            for (size_t i = 0; i < len - k_ + 1; i += max_len_seq - k_ + 1){

                const int curr_len = min(len - i, max_len_seq);
                const char saved_char = str[i + curr_len];
                const char* str_tmp = &str[i];

                str[i + curr_len] = '\0';

                KmerHashIterator<RepHash> it_kmer_h(str_tmp, curr_len, k_), it_kmer_h_end;
                minHashIterator<RepHash> it_min(str_tmp, curr_len, k_, g_, rep, true);

                for (int last_pos_km = -2; it_kmer_h != it_kmer_h_end; ++it_kmer_h, ++it_min) {

                    std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                    if (p_.second != last_pos_km + 1){ // If one or more k-mer were jumped because contained non-ACGT char.

                        km = Kmer(&str_tmp[p_.second]);

                        it_min += (last_pos_km == -2 ? p_.second : (p_.second - last_pos_km) - 1);
                        it_min_h = it_min.getHash();

                        block_bf = bf.getBlock(it_min_h);
                    }
                    else {

                        km.selfForwardBase(str_tmp[p_.second + k_ - 1]);

                        it_min_h = it_min.getHash();
                        if (it_min_h != last_it_min_h) block_bf = bf.getBlock(it_min_h);
                    }

                    last_pos_km = p_.second;
                    last_it_min_h = it_min_h;

                    const size_t r = bf.contains_block(p_.first, it_min_h, block_bf);

                    if (r != 0){

                        while (locks_unitig[thread_id].test_and_set(std::memory_order_acquire));

                        const UnitigMap<U, G> um = findUnitig(km, str_tmp, p_.second);

                        if (um.isEmpty) { // kmer did not map, push into queue for next unitig generation round

                            locks_unitig[thread_id].clear(std::memory_order_release);

                            bool isIsolated = false;

                            string newseq;

                            const size_t pos_match = findUnitigSequenceBBF(km, newseq, isIsolated, l_ignored_km_tips); //Build unitig from Bloom filter

                            if (!reference_mode && isIsolated){ // According to the BF, k-mer is isolated in the graph and is a potential false positive

                                const uint64_t block = (r == 1 ? block_bf.first : block_bf.second);
                                const uint64_t id_lock = block % nb_locks;
                                const Kmer km_rep = km.rep();

                                tiny_vector<Kmer, 2>& v = fp_candidate[block];

                                size_t i = 0;

                                while (locks_fp[id_lock].test_and_set(std::memory_order_acquire));

                                for (; i < v.size(); ++i){ // Search list of fp candidate for k-mer

                                    if (v[i] == km_rep) break;
                                }

                                if (i >= v.size()) v.push_back(km_rep);
                                else {

                                    v.remove(i);

                                    addUnitigSequenceBBF(km, newseq, pos_match, 1, locks_mapping, locks_unitig, thread_id);
                                    addUnitigSequenceBBF(km, newseq, pos_match, 1, locks_mapping, locks_unitig, thread_id);
                                }

                                locks_fp[id_lock].clear(std::memory_order_release);
                            }
                            else {

                                const size_t len_match_km = 1 + cstrMatch(&str_tmp[p_.second + k_], &(newseq.c_str()[pos_match + k_]));

                                addUnitigSequenceBBF(km, newseq, pos_match, len_match_km, locks_mapping, locks_unitig, thread_id);

                                it_kmer_h += len_match_km - 1;
                            }
                        }
                        else {

                            const uint64_t id_lock = um.pos_unitig % nb_locks;

                            while (locks_mapping[id_lock].test_and_set(std::memory_order_acquire));

                            mapRead(um);

                            locks_mapping[id_lock].clear(std::memory_order_release);
                            locks_unitig[thread_id].clear(std::memory_order_release);

                            it_kmer_h += um.len - 1;
                        }
                    }
                }

                str[i + curr_len] = saved_char;
            }

            str += len + 1;
        }

        while (lock_ignored_km_tips.test_and_set(std::memory_order_acquire));

        for (const auto& km_tip : l_ignored_km_tips) ignored_km_tips.insert(km_tip, false);

        lock_ignored_km_tips.clear(std::memory_order_release);
    };

    auto reading_function = [&](char* seq_buf, size_t& seq_buf_sz) {

        size_t file_id = 0;

        const size_t sz_buf = thread_seq_buf_sz - k_;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_buf) {

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

                pos_read = (new_reading ? 0 : pos_read);
                len_read = s.length();

                s_str = s.c_str();

                if (len_read >= k_){

                    if ((thread_seq_buf_sz - seq_buf_sz - 1) < (len_read - pos_read)){

                        strncpy(&seq_buf[seq_buf_sz], &s_str[pos_read], thread_seq_buf_sz - seq_buf_sz - 1);

                        seq_buf[thread_seq_buf_sz - 1] = '\0';

                        pos_read += sz_buf - seq_buf_sz;
                        seq_buf_sz = thread_seq_buf_sz;

                        break;
                    }
                    else {

                        strcpy(&seq_buf[seq_buf_sz], &s_str[pos_read]);

                        seq_buf_sz += (len_read - pos_read) + 1;
                        pos_read = len_read;
                    }
                }
                else pos_read = len_read;
            }
            else return true;
        }

        return false;
    };

    {
        bool stop = false;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_file;

        char* buffer_seq = new char[opt.nb_threads * thread_seq_buf_sz];
        size_t* buffer_seq_sz = new size_t[opt.nb_threads];

        if (opt.verbose) cout << "CompactedDBG::construct(): Extract approximate unitigs" << endl;

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_file);

                            if (stop) return;

                            stop = reading_function(&buffer_seq[t * thread_seq_buf_sz], buffer_seq_sz[t]);
                        }

                        worker_function(&buffer_seq[t * thread_seq_buf_sz], buffer_seq_sz[t], t);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        delete[] buffer_seq;
        delete[] buffer_seq_sz;
    }

    fp.close();

    bf.clear();
    locks_unitig.clear();
    locks_mapping.clear();
    locks_fp.clear();

    if (fp_candidate != nullptr) delete[] fp_candidate;

    if (opt.verbose) cout << "CompactedDBG::construct(): Closed all input files" << endl;

    const size_t unitigsBefore = size();

    if (opt.verbose) cout << endl << "CompactedDBG::construct(): Splitting unitigs (1/2)" << endl;

    pair<size_t, size_t> unitigSplit = extractAllUnitigs();

    const int unitigsAfter1 = size();

    if (opt.verbose) cout << endl << "CompactedDBG::construct(): Splitting unitigs (2/2)" << endl;

    check_fp_tips(ignored_km_tips);
    ignored_km_tips.clear_tables();

    const int unitigsAfter2 = size();

    if (opt.verbose) {

        cout << "CompactedDBG::construct(): Before split: " << unitigsBefore << " unitigs" << endl;
        cout << "CompactedDBG::construct(): After split (1/" << (reference_mode ? "1" : "2" ) << "): " << unitigsAfter1 << " unitigs" <<  endl;
        if (!reference_mode) cout << "CompactedDBG::construct(): After split (2/2): " << unitigsAfter2 << " unitigs" <<  endl;
        cout << "CompactedDBG::construct(): Unitigs split: " << unitigSplit.first << endl;
        cout << "CompactedDBG::construct(): Unitigs deleted: " << unitigSplit.second << endl;

        cout << endl << "CompactedDBG::construct(): Joining unitigs" << endl;
    }

    const size_t joined = joinUnitigs_<is_void<U>::value>(nullptr, opt.nb_threads);

    const int unitigsAfter3 = size();

    if (opt.verbose) {

        cout << "CompactedDBG::construct(): After join: " << unitigsAfter3 << " unitigs" << endl;
        cout << "CompactedDBG::construct(): Joined " << joined << " unitigs" << endl;
    }

    if (opt.useMercyKmers && !reference_mode){

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
template<typename U, typename G>
bool CompactedDBG<U, G>::addUnitigSequenceBBF(Kmer km, const string& seq, const size_t pos_match_km, const size_t len_match_km,
                                           vector<std::atomic_flag>& locks_mapping, vector<std::atomic_flag>& locks_unitig,
                                           const size_t thread_id) {

    for (auto& lck : locks_unitig){ //Acquire all the locks for insertion

        while (lck.test_and_set(std::memory_order_acquire));
    }

    UnitigMap<U, G> um = find(km); // Look if unitig was already inserted

    if (um.isEmpty){ // If it wasn't already inserted

        addUnitig(seq, seq.length() == k_ ? v_kmers.size() : v_unitigs.size()); // Add the unitig

        for (size_t t = 0; t < locks_unitig.size(); ++t){
            // Release all locks except one
            if (t != thread_id) locks_unitig[t].clear(std::memory_order_release);
        }

        um = find(km); // Get location of the inserted unitig
    }
    else {

        for (size_t t = 0; t < locks_unitig.size(); ++t){
            // Release all locks except one
            if (t != thread_id) locks_unitig[t].clear(std::memory_order_release);
        }
    }

    // Prepare read mapping
    um.len = len_match_km;

    if (!um.isShort && !um.isAbundant && !um.strand) um.dist -= um.len - 1;

    if (um.dist + um.len > um.size - k_){ // This is a self loop

        KmerIterator it(&(seq.c_str()[pos_match_km]));

        for (size_t i = 0; i != len_match_km; ++it, ++i){

            um = find(it->first);

            const uint64_t id_lock = um.pos_unitig % locks_mapping.size();

            while (locks_mapping[id_lock].test_and_set(std::memory_order_acquire));

            mapRead(find(it->first)); // Map k-mers one by one

            locks_mapping[id_lock].clear(std::memory_order_release);
        }
    }
    else {

        const uint64_t id_lock = um.pos_unitig % locks_mapping.size();

        while (locks_mapping[id_lock].test_and_set(std::memory_order_acquire));

        mapRead(um); // Map read

        locks_mapping[id_lock].clear(std::memory_order_release);
    }

    locks_unitig[thread_id].clear(std::memory_order_release);

    return !um.isEmpty;
}

// use:  cm.findUnitigSequenceBBF(km, s, selfLoop)
// pre:  km is in the bloom filter
// post: s is the unitig containing the kmer km
//       and the first k-mer in s is smaller (wrt. < operator)
//       than the last kmer
//       selfLoop is true of the unitig is a loop or hairpin
template<typename U, typename G>
size_t CompactedDBG<U, G>::findUnitigSequenceBBF(Kmer km, string& s, bool& isIsolated, vector<Kmer>& l_ignored_km_tip) {

    string fw_s;

    Kmer end = km;
    Kmer last = end;
    Kmer twin = km.twin();

    char c;

    size_t j = 0;

    bool has_no_neighbor, selfLoop = false;

    isIsolated = false;

    while (fwStepBBF(end, end, c, has_no_neighbor, l_ignored_km_tip)) {

        ++j;

        if (end == km) {
            selfLoop = true;
            break;
        }
        else if ((end == twin) || (end == last.twin())) break;

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

            if ((front == km) || (front == twin) || (front == first.twin())) break;

            bw_s.push_back(c);
            first = front;
        }

        if (isIsolated) isIsolated = (j == 0) && has_no_neighbor;

        reverse(bw_s.begin(), bw_s.end());
    }

    char tmp[MAX_KMER_SIZE];

    km.toString(tmp);

    s.reserve(k_ + fw_s.size() + bw_s.size());

    s.append(bw_s);
    s.append(tmp);
    s.append(fw_s);

    return bw_s.size();
}

template<typename U, typename G>
bool CompactedDBG<U, G>::bwStepBBF(Kmer km, Kmer& front, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand) const {

    size_t i, j = -1, j_tmp;
    size_t nb_neigh = 0;

    char km_tmp[MAX_KMER_SIZE];

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

template<typename U, typename G>
bool CompactedDBG<U, G>::fwStepBBF(Kmer km, Kmer& end, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand) const {

    size_t i, j = -1, j_tmp;
    size_t nb_neigh = 0;

    char km_tmp[MAX_KMER_SIZE];

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
template<typename U, typename G>
UnitigMap<U, G> CompactedDBG<U, G>::findUnitig(const Kmer& km, const char* s, size_t pos) {

    // need to check if we find it right away, need to treat this common case
    const UnitigMap<U, G> cc = find(km);

    if (!cc.isEmpty && !cc.isShort && !cc.isAbundant){

        const CompressedSequence& seq = v_unitigs[cc.pos_unitig]->seq;
        size_t km_dist = cc.dist;
        size_t jlen = 0;

        if (cc.strand) jlen = seq.jump(s, pos, cc.dist, false) - k_ + 1;
        else {

            jlen = seq.jump(s, pos, cc.dist + k_ - 1, true) - k_ + 1; // match s_fw to comp(seq)_bw
            km_dist -= jlen - 1;
        }

        return UnitigMap<U, G>(cc.pos_unitig, km_dist, jlen, cc.size, false, false, cc.strand, this);
    }

    return cc;
}

template<typename U, typename G>
UnitigMap<U, G> CompactedDBG<U, G>::findUnitig(const Kmer& km, const char* s, size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h) {

    // need to check if we find it right away, need to treat this common case
    const UnitigMap<U, G> cc = find(km, it_min_h);

    if (!cc.isEmpty && !cc.isShort && !cc.isAbundant){

        const CompressedSequence& seq = v_unitigs[cc.pos_unitig]->seq;
        size_t km_dist = cc.dist;
        size_t jlen = 0;

        if (cc.strand) jlen = seq.jump(s, pos, cc.dist, false) - k_ + 1;
        else {

            jlen = seq.jump(s, pos, cc.dist + k_ - 1, true) - k_ + 1; // match s_fw to comp(seq)_bw
            km_dist -= jlen - 1;
        }

        return UnitigMap<U, G>(cc.pos_unitig, km_dist, jlen, cc.size, false, false, cc.strand, this);
    }

    return cc;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::addUnitig(const string& str_unitig, const size_t id_unitig){

    int pos;

    const size_t len = str_unitig.size();
    size_t pos_id_unitig = id_unitig << 32;

    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    bool isShort = false, isAbundant = false, isForbidden = false;

    const char* c_str = str_unitig.c_str();

    char km_tmp[MAX_KMER_SIZE];

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

                std::pair<hmap_min_unitigs_t::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                size_t flag = p.first.getVal2();
                size_t v_sz = p.first->size(flag);

                pos = min_h_res.pos;
                pos_id_unitig = (pos_id_unitig & mask) | ((size_t) pos);

                if (!isShort){

                    mhr = min_h_res;

                    while ((v_sz >= max_abundance_lim) || ((v_sz > 0) && (((*(p.first))(v_sz-1, flag) & mask) == mask))){

                        mhr_tmp = it_min.getNewMin(mhr);
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            if (((*(p.first))(v_sz-1, flag) & mask) != mask){ // Minimizer was never signaled before as overcrowded

                                // If minimizer bin already contains abundant k-mer, just set flag for unitig overcrowding
                                if (((*(p.first))(v_sz-1, flag) & MASK_CONTIG_ID) == MASK_CONTIG_ID) (*(p.first))(v_sz-1, flag) |= MASK_CONTIG_TYPE;
                                else p.first.getVal2() = (flag = p.first->push_back(mask, flag));
                            }

                            mhr = mhr_tmp;
                            minz_rep = Minimizer(&c_str[mhr.pos]).rep();

                            p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                            flag = p.first.getVal2();
                            v_sz = p.first->size(flag);
                        }
                        else break;
                    }
                }

                packed_tiny_vector& v = p.first.getVal1();
                uint8_t& flag_v = p.first.getVal2();

                if (v_sz == 0) flag_v = v.push_back(pos_id_unitig, flag_v); //Newly created vector, just push unitig ID
                else if (isShort && (v_sz >= min_abundance_lim)){ //The minimizer (is/might be) too abundant

                    isShort = false;
                    isAbundant = true;

                    it_min = it_min_end;

                    break;
                }
                else if ((v(v_sz-1, flag_v) & MASK_CONTIG_ID) == MASK_CONTIG_ID){ //If minimizer is abundant or crowded

                    if ((v_sz == 1) || (v(v_sz-2, flag_v) != pos_id_unitig)) flag_v = v.insert(pos_id_unitig, v_sz-1, flag_v);
                }
                else if (v(v_sz-1, flag_v) != pos_id_unitig) flag_v = v.push_back(pos_id_unitig, flag_v);

                last_pos_min = min_h_res.pos;
                ++it_it_min;
            }
        }
    }

    if (isAbundant){

        if (id_unitig == v_kmers.size()) v_kmers.push_back(make_pair(km_rep, CompressedCoverage_t<U>(1)));
        else v_kmers[id_unitig] = make_pair(km_rep, CompressedCoverage_t<U>(1));

        deleteUnitig_<is_void<U>::value>(true, false, id_unitig, false);
        if (id_unitig == v_kmers.size() - 1) v_kmers.resize(v_kmers.size() - 1);

        it_min = minHashIterator<RepHash>(c_str, len, k_, g_, RepHash(), true);

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){

            if (last_pos_min < it_min.getPosition()){ //If current minimizer was not seen before

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&c_str[min_h_res.pos]).rep(); //Get the minimizer to insert

                    std::pair<hmap_min_unitigs_t::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                    packed_tiny_vector& v = p.first.getVal1();
                    uint8_t& flag_v = p.first.getVal2();
                    const size_t v_sz = v.size(flag_v);

                    if ((v_sz > 0) && ((v(v_sz-1, flag_v) & MASK_CONTIG_ID) == MASK_CONTIG_ID)) v(v_sz-1, flag_v)++;
                    else flag_v = v.push_back(MASK_CONTIG_ID + 1, flag_v);

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }

        h_kmers_ccov.insert(km_rep, CompressedCoverage_t<U>(1));
    }
    else if (isShort){

        if (id_unitig == v_kmers.size()) v_kmers.push_back(make_pair(km_rep,  CompressedCoverage_t<U>(1)));
        else v_kmers[id_unitig] = make_pair(km_rep,  CompressedCoverage_t<U>(1));
    }
    else if (id_unitig == v_unitigs.size()) v_unitigs.push_back(new Unitig<U>(c_str)); //Push unitig to list of unitigs
    else v_unitigs[id_unitig] = new Unitig<U>(c_str);

    return isAbundant;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::addUnitig(const string& str_unitig, const size_t id_unitig, const size_t id_unitig_r, const size_t is_short_r){

    int pos;

    const size_t pos_id_unitig_r = (id_unitig_r << 32) | (is_short_r ? MASK_CONTIG_TYPE : 0);

    const size_t len = str_unitig.size();
    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    size_t pos_id_unitig = id_unitig << 32;
    size_t i = 0;

    bool isShort = false, isAbundant = false, isForbidden = false;

    const char* c_str = str_unitig.c_str();

    char km_tmp[MAX_KMER_SIZE];

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

                std::pair<hmap_min_unitigs_t::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                size_t flag = p.first.getVal2();
                size_t v_sz = p.first->size(flag);

                pos = min_h_res.pos;
                pos_id_unitig = (pos_id_unitig & mask) | ((size_t) pos);

                for (i = 0; i < v_sz; ++i){

                    if (((*(p.first))(i, flag) & mask) == pos_id_unitig_r){

                        (*(p.first))(i, flag) = pos_id_unitig;
                        break;
                    }
                }

                if (!isShort){

                    mhr = min_h_res;

                    while ((i == v_sz) && ((v_sz >= max_abundance_lim) || ((v_sz > 0) && (((*(p.first))(v_sz-1, flag) & mask) == mask)))){

                        mhr_tmp = it_min.getNewMin(mhr);
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            if (((*(p.first))(v_sz-1, flag) & mask) != mask){ // Minimizer was never signaled before as overcrowded

                                // If minimizer bin already contains abundant k-mer, just set flag for unitig overcrowding
                                if (((*(p.first))(v_sz-1, flag) & MASK_CONTIG_ID) == MASK_CONTIG_ID) (*(p.first))(v_sz-1, flag) |= MASK_CONTIG_TYPE;
                                else p.first.getVal2() = (flag = p.first->push_back(mask, flag));
                            }

                            mhr = mhr_tmp;
                            minz_rep = Minimizer(&c_str[mhr.pos]).rep();

                            p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                            flag = p.first.getVal2();
                            v_sz = p.first->size(flag);

                            for (i = 0; i < v_sz; ++i){

                                if (((*(p.first))(i, flag) & mask) == pos_id_unitig_r){

                                    (*(p.first))(i, flag) = pos_id_unitig;
                                    break;
                                }
                            }
                        }
                        else {

                            i = v_sz;
                            break;
                        }
                    }
                }

                if (i == v_sz){

                    packed_tiny_vector& v = p.first.getVal1();
                    uint8_t& flag_v = p.first.getVal2();

                    if (v_sz == 0) flag_v = v.push_back(pos_id_unitig, flag_v); //Newly created vector, just push unitig ID
                    else if (isShort && (v_sz >= min_abundance_lim)){ //The minimizer (is/might be) too abundant

                        isShort = false;
                        isAbundant = true;

                        it_min = it_min_end;

                        break;
                    }
                    else if ((v(v_sz-1, flag_v) & MASK_CONTIG_ID) == MASK_CONTIG_ID){ //If minimizer is abundant or crowded

                        if ((v_sz == 1) || (v(v_sz-2, flag_v) != pos_id_unitig)) flag_v = v.insert(pos_id_unitig, v_sz-1, flag_v);
                    }
                    else if (v(v_sz-1, flag_v) != pos_id_unitig) flag_v = v.push_back(pos_id_unitig, flag_v);
                }

                last_pos_min = min_h_res.pos;
                ++it_it_min;
            }
        }
    }

    if (isAbundant){

        if (id_unitig == v_kmers.size()) v_kmers.push_back(make_pair(km_rep, CompressedCoverage_t<U>(1)));
        else v_kmers[id_unitig] = make_pair(km_rep, CompressedCoverage_t<U>(1));

        deleteUnitig_<is_void<U>::value>(true, false, id_unitig, false);
        if (id_unitig == v_kmers.size() - 1) v_kmers.resize(v_kmers.size() - 1);

        it_min = minHashIterator<RepHash>(c_str, len, k_, g_, RepHash(), true);

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){

            if (last_pos_min < it_min.getPosition()){ //If current minimizer was not seen before

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&c_str[min_h_res.pos]).rep(); //Get the minimizer to insert

                    std::pair<hmap_min_unitigs_t::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                    packed_tiny_vector& v = p.first.getVal1();
                    uint8_t& flag_v = p.first.getVal2();
                    const size_t v_sz = v.size(flag_v);

                    if ((v_sz > 0) && ((v(v_sz-1, flag_v) & MASK_CONTIG_ID) == MASK_CONTIG_ID)) v(v_sz-1, flag_v)++;
                    else flag_v = v.push_back(MASK_CONTIG_ID + 1, flag_v);

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }

        h_kmers_ccov.insert(km_rep, CompressedCoverage_t<U>(1));
    }
    else if (isShort){

        if (id_unitig == v_kmers.size()) v_kmers.push_back(make_pair(km_rep,  CompressedCoverage_t<U>(1)));
        else v_kmers[id_unitig] = make_pair(km_rep, CompressedCoverage_t<U>(1));
    }
    else if (id_unitig == v_unitigs.size()) v_unitigs.push_back(new Unitig<U>(c_str)); //Push unitig to list of unitigs
    else v_unitigs[id_unitig] = new Unitig<U>(c_str);

    return isAbundant;
}

template<typename U, typename G>
void CompactedDBG<U, G>::swapUnitigs(const bool isShort, const size_t id_a, const size_t id_b){

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
                hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(minz_rep);

                if (it_h != hmap_min_unitigs.end()){

                    if (minimizers.find(minz_rep) == minimizers.end()){

                        minimizers.insert(minz_rep, 0);

                        packed_tiny_vector& v = it_h.getVal1();
                        const uint8_t flag_v = it_h.getVal2();
                        const size_t v_sz = v.size(flag_v);

                        for (size_t i = 0; i < v_sz; ++i){

                             // Swap the unitig ids but do not change positions;
                            if ((v(i, flag_v) & mask) == shift_id_unitig_b) v(i, flag_v) = shift_id_unitig_a | (v(i, flag_v) & MASK_CONTIG_POS);
                            else if ((v(i, flag_v) & mask) == shift_id_unitig_a) v(i, flag_v) = shift_id_unitig_b | (v(i, flag_v) & MASK_CONTIG_POS);
                        }
                    }

                    //size_t v_sz = it_h->size();
                    packed_tiny_vector* v = &(it_h.getVal1());
                    uint8_t flag_v = it_h.getVal2();
                    size_t v_sz = v->size(flag_v);

                    mhr = min_h_res;

                    //while (((*it_h)[v_sz-1] & mask) == mask){
                    while (((*v)(v_sz-1, flag_v) & mask) == mask){

                        mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            minz_rep = Minimizer(&s[mhr_tmp.pos]).rep();
                            it_h = hmap_min_unitigs.find(minz_rep);

                            if (it_h == hmap_min_unitigs.end()) break;

                            if (minimizers.find(minz_rep) == minimizers.end()){

                                minimizers.insert(minz_rep, 0);

                                packed_tiny_vector& v_tmp = it_h.getVal1();
                                const uint8_t flag_v_tmp = it_h.getVal2();
                                const size_t v_sz_tmp = v_tmp.size(flag_v_tmp);

                                for (size_t i = 0; i < v_sz_tmp; ++i){

                                     // Swap the unitig ids but do not change positions;
                                    if ((v_tmp(i, flag_v_tmp) & mask) == shift_id_unitig_b){

                                        v_tmp(i, flag_v_tmp) = shift_id_unitig_a | (v_tmp(i, flag_v_tmp) & MASK_CONTIG_POS);
                                    }
                                    else if ((v_tmp(i, flag_v_tmp) & mask) == shift_id_unitig_a){

                                        v_tmp(i, flag_v_tmp) = shift_id_unitig_b | (v_tmp(i, flag_v_tmp) & MASK_CONTIG_POS);
                                    }
                                }
                            }

                            mhr = mhr_tmp;

                            v = &(it_h.getVal1());
                            flag_v = it_h.getVal2();
                            v_sz = v->size(flag_v);
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
                hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(minz_rep);

                if (it_h != hmap_min_unitigs.end()){

                    if (minimizers.find(minz_rep) == minimizers.end()){

                        minimizers.insert(minz_rep, 0);

                        packed_tiny_vector& v = it_h.getVal1();
                        const uint8_t flag_v = it_h.getVal2();
                        const size_t v_sz = v.size(flag_v);

                        for (size_t i = 0; i < v_sz; ++i){

                            if ((v(i, flag_v) & mask) == shift_id_unitig_a) v(i, flag_v) = shift_id_unitig_b | (v(i, flag_v) & MASK_CONTIG_POS);
                        }
                    }

                    packed_tiny_vector* v = &(it_h.getVal1());
                    uint8_t flag_v = it_h.getVal2();
                    size_t v_sz = v->size(flag_v);

                    mhr = min_h_res;

                    while (((*v)(v_sz-1, flag_v) & mask) == mask){

                        mhr_tmp = it_min2.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            minz_rep = Minimizer(&s[mhr_tmp.pos]).rep();
                            it_h = hmap_min_unitigs.find(minz_rep);

                            if (it_h == hmap_min_unitigs.end()) break;

                            if (minimizers.find(minz_rep) == minimizers.end()){

                                minimizers.insert(minz_rep, 0);

                                packed_tiny_vector& v_tmp = it_h.getVal1();
                                const uint8_t flag_v_tmp = it_h.getVal2();
                                const size_t v_sz_tmp = v_tmp.size(flag_v_tmp);

                                for (size_t i = 0; i < v_sz_tmp; ++i){

                                    if ((v_tmp(i, flag_v_tmp) & mask) == shift_id_unitig_a) {

                                        v_tmp(i, flag_v_tmp) = shift_id_unitig_b | (v_tmp(i, flag_v_tmp) & MASK_CONTIG_POS);
                                    }
                                }
                            }

                            mhr = mhr_tmp;
                            v = &(it_h.getVal1());
                            flag_v = it_h.getVal2();
                            v_sz = v->size(flag_v);
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

// Input sequence must not contain duplicated k-mers
// If it does, use CompactedDBG<U, G>::add().
template<typename U, typename G>
bool CompactedDBG<U, G>::mergeUnitig(const string& seq, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::mergeUnitig(): Graph is invalid and no sequence can be added to it" << endl;
        return false;
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::mergeUnitig(): Input sequence length cannot be less than k = " << k_ << endl;
        return false;
    }

    size_t nxt_pos_insert_v_unitigs = v_unitigs.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t v_kmers_sz = v_kmers.size();

    size_t nb_curr_pred = 0;
    size_t nb_curr_succ = 0;
    size_t nb_prev_pred = 0xFFFFFFFFFFFFFFFF;
    size_t nb_prev_succ = 0xFFFFFFFFFFFFFFFF;

    size_t added = 0;
    size_t split_before = 0, split_after;

    bool prev_found = true;

    string curr_unitig;

    vector<Kmer> v_joins;

    KmerHashTable<vector<size_t>> kht;

    const char* str_seq = seq.c_str();

    const size_t str_seq_len = seq.length();

    auto add_graph_function = [&](const string& unitig){

        const char* str_unitig = unitig.c_str();
        const size_t len_unitig = unitig.length();

        if (len_unitig == k_){

            if (!addUnitig(str_unitig, v_kmers_sz)){

                v_kmers[v_kmers_sz].second.ccov.setFull();

                ++v_kmers_sz;
                ++added;
            }
            else h_kmers_ccov.find(Kmer(str_unitig).rep())->ccov.setFull();
        }
        else {

            addUnitig(str_unitig, nxt_pos_insert_v_unitigs);

            v_unitigs[nxt_pos_insert_v_unitigs]->ccov.setFull();

            ++nxt_pos_insert_v_unitigs;
            ++added;

            v_joins.push_back(Kmer(&str_unitig[len_unitig - k_]));
        }

        v_joins.push_back(Kmer(str_unitig));
    };

    auto add_split_function = [&](){

        KmerHashTable<vector<size_t>>::iterator it(kht.begin());
        const KmerHashTable<vector<size_t>>::iterator it_end(kht.end());

        vector<pair<int,int>> sp;

        while (it != it_end){

            ++split_before;
            ++split_after;

            const_UnitigMap<U, G> um(find(it.getKey(), true));

            vector<size_t>& split_v = *it;

            sort(split_v.begin(), split_v.end());

            size_t prev_split_pos = 0;

            for (const auto pos : split_v){

                if (pos != prev_split_pos){

                    sp.push_back({prev_split_pos, pos});

                    v_joins.push_back(um.getUnitigKmer(prev_split_pos));
                    if (prev_split_pos != pos - 1) v_joins.push_back(um.getUnitigKmer(pos - 1));

                    prev_split_pos = pos;

                    ++split_after;
                }
            }

            sp.push_back({prev_split_pos, um.size - k_ + 1});

            v_joins.push_back(um.getUnitigKmer(prev_split_pos));
            if (prev_split_pos != um.size - k_) v_joins.push_back(um.getUnitigKmer(um.size - k_));

            extractUnitig_<is_void<U>::value>(um.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

            sp.clear();
            split_v.clear();

            ++it;
        }

        kht.clear();
    };

    for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end;) { //non-ACGT char. are discarded

        const std::pair<Kmer, int>& p = *it_km;

        UnitigMap<U, G> cm = findUnitig(p.first, str_seq, p.second);

        if (cm.isEmpty){

            vector<const_UnitigMap<U, G>> um_bw(findPredecessors(p.first));
            vector<const_UnitigMap<U, G>> um_fw(findSuccessors(p.first));

            nb_curr_pred = 0;
            nb_curr_succ = 0;

            for (auto& um : um_bw){

                if (!um.isEmpty){

                    ++nb_curr_pred;

                    if (!um.isAbundant && !um.isShort){

                        if (um.strand) ++(um.dist);
                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)) kht.insert(um.getUnitigHead(), vector<size_t>()).first->push_back(um.dist);
                    }
                }
            }

            for (auto& um : um_fw){

                if (!um.isEmpty){

                    ++nb_curr_succ;

                    if (!um.isAbundant && !um.isShort){

                        if (!um.strand) ++(um.dist);
                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)) kht.insert(um.getUnitigHead(), vector<size_t>()).first->push_back(um.dist);
                    }
                }
            }

            if (prev_found){ // Previous k-mer was found in the graph => current k-mer (not found in the graph) starts a new unitig

                prev_found = false;
                curr_unitig = p.first.toString();
            }
            else if ((nb_prev_succ == 0) && (nb_curr_pred == 0)) curr_unitig.push_back(str_seq[p.second + k_ - 1]);
            else {

                add_graph_function(curr_unitig);

                curr_unitig = p.first.toString();
            }

            nb_prev_succ = nb_curr_succ;
            nb_prev_pred = nb_curr_pred;

            ++it_km;
        }
        else {

            if ((p.second == 0) && !cm.isAbundant && !cm.isShort){

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)) kht.insert(cm.getUnitigHead(), vector<size_t>()).first->push_back(cm.dist);
            }
            else if ((p.second + cm.len == str_seq_len - k_ + 1) && !cm.isAbundant && !cm.isShort){

                cm.dist += cm.len;

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)) kht.insert(cm.getUnitigHead(), vector<size_t>()).first->push_back(cm.dist);
            }

            it_km += cm.len;

            if (!prev_found){

                add_graph_function(curr_unitig);

                prev_found = true;
            }
        }
    }

    add_split_function();

    if (!prev_found) add_graph_function(curr_unitig);

    if (nxt_pos_insert_v_unitigs < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert_v_unitigs);
    if (v_kmers_sz < v_kmers.size()) v_kmers.resize(v_kmers_sz);

    const size_t joined = joinUnitigs_<is_void<U>::value>(&v_joins);

    if (verbose){

        cout << "CompactedDBG::mergeUnitig(): Added " << added << " new unitigs to the graph." << endl;
        cout << "CompactedDBG::mergeUnitig(): Split " << split_before << " unitigs into " << split_after << " new unitigs." << endl;
        cout << "CompactedDBG::mergeUnitig(): Joined " << joined << " unitigs from the graph." << endl;
    }

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::annotateSplitUnitig(const string& seq, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::annotateSplitUnitig(): Graph is invalid and no sequence can be added to it" << endl;
        return false;
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::annotateSplitUnitig(): Input sequence length cannot be less than k = " << k_ << endl;
        return false;
    }

    size_t nxt_pos_insert_v_unitigs = v_unitigs.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t v_kmers_sz = v_kmers.size();

    size_t nb_curr_pred = 0;
    size_t nb_curr_succ = 0;
    size_t nb_prev_pred = 0xFFFFFFFFFFFFFFFF;
    size_t nb_prev_succ = 0xFFFFFFFFFFFFFFFF;

    bool prev_found = true;

    string curr_unitig;

    const char* str_seq = seq.c_str();

    const size_t str_seq_len = seq.length();

    auto add_graph_function = [&](const string& unitig){

        const char* str_unitig = unitig.c_str();
        const size_t len_unitig = unitig.length();

        if (len_unitig == k_){

            if (!addUnitig(str_unitig, v_kmers_sz)) v_kmers[v_kmers_sz++].second.ccov.setFull();
            else h_kmers_ccov.find(Kmer(str_unitig).rep())->ccov.setFull();
        }
        else {

            addUnitig(str_unitig, nxt_pos_insert_v_unitigs);

            v_unitigs[nxt_pos_insert_v_unitigs++]->ccov.setFull();
        }
    };

    for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end;) { //non-ACGT char. are discarded

        const std::pair<Kmer, int>& p = *it_km;

        UnitigMap<U, G> cm = findUnitig(p.first, str_seq, p.second);

        if (cm.isEmpty){

            vector<const_UnitigMap<U, G>> um_bw(findPredecessors(p.first));
            vector<const_UnitigMap<U, G>> um_fw(findSuccessors(p.first));

            nb_curr_pred = 0;
            nb_curr_succ = 0;

            for (auto& um : um_bw){

                if (!um.isEmpty){

                    ++nb_curr_pred;

                    if (!um.isAbundant && !um.isShort){

                        if (um.strand) ++(um.dist);
                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)) unmapRead(um);
                    }
                }
            }

            for (auto& um : um_fw){

                if (!um.isEmpty){

                    ++nb_curr_succ;

                    if (!um.isAbundant && !um.isShort){

                        if (!um.strand) ++(um.dist);
                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)) unmapRead(um);
                    }
                }
            }

            if (prev_found){ // Previous k-mer was found in the graph => current k-mer (not found in the graph) starts a new unitig

                prev_found = false;
                curr_unitig = p.first.toString();
            }
            else if (((nb_prev_succ == 0) && (nb_curr_pred == 0))) curr_unitig.push_back(str_seq[p.second + k_ - 1]);
            else {

                add_graph_function(curr_unitig);

                curr_unitig = p.first.toString();
            }

            nb_prev_succ = nb_curr_succ;
            nb_prev_pred = nb_curr_pred;

            ++it_km;
        }
        else {

            if ((p.second == 0) && !cm.isAbundant && !cm.isShort){

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)){

                    const size_t len = cm.len;

                    cm.len = 1;
                    unmapRead(cm);
                    cm.len = len;
                }
            }
            else if ((p.second + cm.len == str_seq_len - k_ + 1) && !cm.isAbundant && !cm.isShort){

                cm.dist += cm.len;

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)){

                    const size_t len = cm.len;

                    cm.len = 1;
                    unmapRead(cm);
                    cm.len = len;
                }
            }

            it_km += cm.len;

            if (!prev_found){

                add_graph_function(curr_unitig);

                prev_found = true;
            }
        }
    }

    if (!prev_found) add_graph_function(curr_unitig);

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::annotateSplitUnitig(const string& seq, LockGraph& lck_g, const size_t thread_id, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::annotateSplitUnitig(): Graph is invalid and no sequence can be added to it" << endl;
        return false;
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::annotateSplitUnitig(): Input sequence length cannot be less than k = " << k_ << endl;
        return false;
    }

    size_t nb_curr_pred = 0;
    size_t nb_curr_succ = 0;
    size_t nb_prev_pred = 0xFFFFFFFFFFFFFFFF;
    size_t nb_prev_succ = 0xFFFFFFFFFFFFFFFF;

    bool prev_found = true;

    string curr_unitig;

    char km_tmp[MAX_KMER_SIZE];

    const char* str_seq = seq.c_str();

    const size_t str_seq_len = seq.length();

    auto add_graph_function = [&](const string& unitig){

        const char* str_unitig = unitig.c_str();
        const size_t len_unitig = unitig.length();

        lck_g.lock_all_threads();

        if (len_unitig == k_){

            if (!addUnitig(str_unitig, v_kmers.size())) v_kmers[v_kmers.size() - 1].second.ccov.setFull();
            else h_kmers_ccov.find(Kmer(str_unitig).rep())->ccov.setFull();
        }
        else {

            addUnitig(str_unitig, v_unitigs.size());

            v_unitigs[v_unitigs.size() - 1]->ccov.setFull();
        }

        lck_g.unlock_all_threads();
    };

    for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end;) { //non-ACGT char. are discarded

        const std::pair<Kmer, int>& p = *it_km;

        lck_g.lock_thread(thread_id);

        UnitigMap<U, G> cm = findUnitig(p.first, str_seq, p.second);

        if (cm.isEmpty){

            vector<const_UnitigMap<U, G>> um_bw(findPredecessors(p.first));
            vector<const_UnitigMap<U, G>> um_fw(findSuccessors(p.first));

            nb_curr_pred = 0;
            nb_curr_succ = 0;

            for (auto& um : um_bw){

                if (!um.isEmpty){

                    ++nb_curr_pred;

                    if (!um.isAbundant && !um.isShort){

                        if (um.strand) ++(um.dist);
                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)){

                            lck_g.lock_unitig(um.pos_unitig);
                            unmapRead(um);
                            lck_g.unlock_unitig(um.pos_unitig);
                        }
                    }
                }
            }

            for (auto& um : um_fw){

                if (!um.isEmpty){

                    ++nb_curr_succ;

                    if (!um.isAbundant && !um.isShort){

                        if (!um.strand) ++(um.dist);
                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)){

                            lck_g.lock_unitig(um.pos_unitig);
                            unmapRead(um);
                            lck_g.unlock_unitig(um.pos_unitig);
                        }
                    }
                }
            }

            lck_g.unlock_thread(thread_id);

            if (prev_found){ // Previous k-mer was found in the graph => current k-mer (not found in the graph) starts a new unitig

                prev_found = false;
                curr_unitig = p.first.toString();
            }
            else if (((nb_prev_succ == 0) && (nb_curr_pred == 0))) curr_unitig.push_back(str_seq[p.second + k_ - 1]);
            else {

                add_graph_function(curr_unitig);

                curr_unitig = p.first.toString();
            }

            nb_prev_succ = nb_curr_succ;
            nb_prev_pred = nb_curr_pred;

            ++it_km;
        }
        else {

            if ((p.second == 0) && !cm.isAbundant && !cm.isShort){

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)){

                    const size_t len = cm.len;

                    cm.len = 1;

                    lck_g.lock_unitig(cm.pos_unitig);
                    unmapRead(cm);
                    lck_g.unlock_unitig(cm.pos_unitig);

                    cm.len = len;
                }
            }
            else if ((p.second + cm.len == str_seq_len - k_ + 1) && !cm.isAbundant && !cm.isShort){

                cm.dist += cm.len;

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)){

                    const size_t len = cm.len;

                    cm.len = 1;

                    lck_g.lock_unitig(cm.pos_unitig);
                    unmapRead(cm);
                    lck_g.unlock_unitig(cm.pos_unitig);

                    cm.len = len;
                }
            }

            lck_g.unlock_thread(thread_id);

            it_km += cm.len;

            if (!prev_found){

                add_graph_function(curr_unitig);

                prev_found = true;
            }
        }
    }

    if (!prev_found) add_graph_function(curr_unitig);

    return true;
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<!is_void, void>::type CompactedDBG<U, G>::deleteUnitig_(const bool isShort, const bool isAbundant,
                                                                                const size_t id_unitig, const bool delete_data){

    if (isAbundant){

        char km_str[MAX_KMER_SIZE];

        typename h_kmers_ccov_t::iterator it_h = h_kmers_ccov.find(id_unitig);

        const Kmer km = it_h.getKey();

        km.toString(km_str);

        minHashIterator<RepHash> it_min(km_str, k_, k_, g_, RepHash(), true), it_min_end;

        if (delete_data) it_h->getData()->clear(UnitigMap<U, G>(id_unitig, 0, 1, k_, isShort, isAbundant, true, this));

        it_h->ccov.clear();
        h_kmers_ccov.erase(it_h);

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

            if (last_pos_min < it_min.getPosition()){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&km_str[min_h_res.pos]).rep(); // Get canonical minimizer

                    hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(minz_rep); // Look for the minimizer in the hash table

                    if (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        packed_tiny_vector& v = it_h.getVal1();
                        uint8_t& flag_v = it_h.getVal2();

                        const size_t last_pos_v = v.size(flag_v) - 1;

                        v(last_pos_v, flag_v)--;

                        if (((v(last_pos_v, flag_v) & RESERVED_ID) == 0) && ((v(last_pos_v, flag_v) & MASK_CONTIG_TYPE) != MASK_CONTIG_TYPE)){

                            if (last_pos_v == 0) hmap_min_unitigs.erase(minz_rep);
                            else flag_v = v.remove(v.size(flag_v) - 1, flag_v);
                        }
                    }

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }
    }
    else {

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

        const size_t len = str.size();

        minHashIterator<RepHash> it_min(s, len, k_, g_, RepHash(), true), it_min_end;

        minHashResult mhr, mhr_tmp;

        // The unitig is deleted but its space in the unitig vector is not because:
        // 1 - It would change indices in the minimizer hash table
        if (isShort){

            if (delete_data) v_kmers[id_unitig].second.getData()->clear(UnitigMap<U, G>(id_unitig, 0, 1, len, isShort, isAbundant, true, this));
            v_kmers[id_unitig].second.ccov.clear();
            v_kmers[id_unitig].first.set_deleted();
        }
        else {

            if (delete_data) v_unitigs[id_unitig]->getData()->clear(UnitigMap<U, G>(id_unitig, 0, len - k_ + 1, len, isShort, isAbundant, true, this));

            delete v_unitigs[id_unitig];
            v_unitigs[id_unitig] = nullptr;
        }

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

            if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
                isForbidden = false;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep(); // Get canonical minimizer
                    hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(minz_rep); // Look for the minimizer in the hash table

                    mhr = min_h_res;

                    while (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        packed_tiny_vector& v = it_h.getVal1();
                        uint8_t& flag_v = it_h.getVal2();
                        const size_t v_sz = v.size(flag_v);

                        for (size_t i = 0; i < v_sz; i++){

                            if ((v(i, flag_v) & mask) == pos_id_unitig){

                                flag_v = v.remove(i, flag_v);
                                break;
                            }
                        }

                        it_h = hmap_min_unitigs.end();

                        if (v.size(flag_v) == 0) hmap_min_unitigs.erase(minz_rep);
                        else if (!isShort && ((v(v_sz-1, flag_v) & mask) == mask)){ //Minimizer bin is overcrowded


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
    }
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<is_void, void>::type CompactedDBG<U, G>::deleteUnitig_( const bool isShort, const bool isAbundant,
                                                                                const size_t id_unitig, const bool delete_data){

    if (isAbundant){

        char km_str[MAX_KMER_SIZE];

        typename h_kmers_ccov_t::iterator it_h = h_kmers_ccov.find(id_unitig);

        const Kmer km = it_h.getKey();

        km.toString(km_str);

        minHashIterator<RepHash> it_min(km_str, k_, k_, g_, RepHash(), true), it_min_end;

        it_h->ccov.clear();

        h_kmers_ccov.erase(it_h);

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

            if (last_pos_min < it_min.getPosition()){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&km_str[min_h_res.pos]).rep(); // Get canonical minimizer

                    hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(minz_rep); // Look for the minimizer in the hash table

                    if (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        packed_tiny_vector& v = it_h.getVal1();
                        uint8_t& flag_v = it_h.getVal2();

                        const size_t last_pos_v = v.size(flag_v) - 1;

                        v(last_pos_v, flag_v)--;

                        if (((v(last_pos_v, flag_v) & RESERVED_ID) == 0) && ((v(last_pos_v, flag_v) & MASK_CONTIG_TYPE) != MASK_CONTIG_TYPE)){

                            if (last_pos_v == 0) hmap_min_unitigs.erase(minz_rep);
                            else flag_v = v.remove(v.size(flag_v) - 1, flag_v);
                        }
                    }

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }
    }
    else {

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

        const size_t len = str.size();

        minHashIterator<RepHash> it_min(s, len, k_, g_, RepHash(), true), it_min_end;

        minHashResult mhr, mhr_tmp;

        // The unitig is deleted but its space in the unitig vector is not because:
        // 1 - It would change indices in the minimizer hash table
        if (isShort){

            v_kmers[id_unitig].second.ccov.clear();
            v_kmers[id_unitig].first.set_deleted();
        }
        else {

            delete v_unitigs[id_unitig];
            v_unitigs[id_unitig] = nullptr;
        }

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

            if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
                isForbidden = false;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep(); // Get canonical minimizer
                    hmap_min_unitigs_t::iterator it_h = hmap_min_unitigs.find(minz_rep); // Look for the minimizer in the hash table

                    mhr = min_h_res;

                    while (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        packed_tiny_vector& v = it_h.getVal1();
                        uint8_t& flag_v = it_h.getVal2();
                        const size_t v_sz = v.size(flag_v);

                        for (size_t i = 0; i < v_sz; i++){

                            if ((v(i, flag_v) & mask) == pos_id_unitig){

                                flag_v = v.remove(i, flag_v);
                                break;
                            }
                        }

                        it_h = hmap_min_unitigs.end();

                        if (v.size(flag_v) == 0) hmap_min_unitigs.erase(minz_rep);
                        else if (!isShort && ((v(v_sz-1, flag_v) & mask) == mask)){ //Minimizer bin is overcrowded


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
    }
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<!is_void, bool>::type CompactedDBG<U, G>::extractUnitig_(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs,
                                                                               size_t& v_unitigs_sz, size_t& v_kmers_sz, const vector<pair<int,int>>& sp){

    bool deleted = true;

    if (!sp.empty()){

        const Unitig<U>* unitig = v_unitigs[pos_v_unitigs];

        const pair<size_t, size_t> lowpair = unitig->ccov.lowCoverageInfo();

        const size_t totalcoverage = unitig->coveragesum - lowpair.second;
        const size_t ccov_size = unitig->ccov.size();

        const string str = unitig->seq.toString();

        size_t i = 0;

        UnitigMap<U, G> um(pos_v_unitigs, 0, 0, unitig->length(), false, false, true, this);

        vector<Unitig<U>> v_data(sp.size());

        deleted = false;

        for (vector<pair<int,int>>::const_iterator sit = sp.begin(); sit != sp.end(); ++sit, ++i) { //Iterate over created split unitigs

            um.dist = sit->first;
            um.len = sit->second - um.dist;

            if (um.len == 1){

                const string split_str = str.substr(um.dist, um.len + k_ - 1);

                um.strand = (split_str <= reverse_complement(split_str));
            }
            else um.strand = true;

            v_data[i] = std::move(um.splitData(sit+1 == sp.end()));
        }

        i = 0;

        for (vector<pair<int,int>>::const_iterator sit = sp.begin(); sit != sp.end(); ++sit, ++i) { //Iterate over created split unitigs

            const size_t len = sit->second - sit->first;
            //const uint64_t cov_tmp = (totalcoverage * len) / (ccov_size - lowpair.first); // Split unitig coverage
            const uint64_t cov_tmp = len * CompressedCoverage::getFullCoverage();

            string split_str = str.substr(sit->first, len + k_ - 1); // Split unitig sequence

            if (len == 1){

                const string split_str_rev = reverse_complement(split_str);

                if (split_str > split_str_rev) split_str = split_str_rev;

                if (addUnitig(split_str, v_kmers_sz)){

                    CompressedCoverage_t<U>& cc_t = *h_kmers_ccov.find(Kmer(split_str.c_str()).rep());

                    cc_t.ccov.setFull();

                    *(cc_t.getData()) = move(*(v_data[i].getData()));
                }
                else {

                    v_kmers[v_kmers_sz].second.ccov.setFull(); // We don't care about the coverage per k-mer anymore

                    *(v_kmers[v_kmers_sz].second.getData()) = move(*(v_data[i].getData()));

                    ++v_kmers_sz;
                }
            }
            else {

                addUnitig(split_str, nxt_pos_insert_v_unitigs);

                v_unitigs[nxt_pos_insert_v_unitigs]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                v_unitigs[nxt_pos_insert_v_unitigs]->coveragesum = cov_tmp;

                *(v_unitigs[nxt_pos_insert_v_unitigs]->getData()) = move(*(v_data[i].getData()));

                ++nxt_pos_insert_v_unitigs;
            }
        }
    }

    --nxt_pos_insert_v_unitigs; //Position of the last unitig in the vector which is not NULL

    if (pos_v_unitigs != nxt_pos_insert_v_unitigs){ // Do not proceed to swap if swap positions are the same

        swapUnitigs(false, pos_v_unitigs, nxt_pos_insert_v_unitigs); // Swap unitigs

        // If the swapped unitig, previously in position nxt_pos_insert, was a split unitig
        // created in this method, do not try to split it again
        if (nxt_pos_insert_v_unitigs >= v_unitigs_sz) ++pos_v_unitigs;
        else --v_unitigs_sz;
    }
    else --v_unitigs_sz;

    deleteUnitig_<false>(false, false, nxt_pos_insert_v_unitigs);

    return deleted;
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<is_void, bool>::type CompactedDBG<U, G>::extractUnitig_(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs,
                                                                              size_t& v_unitigs_sz, size_t& v_kmers_sz, const vector<pair<int,int>>& sp){

    bool deleted = true;

    if (!sp.empty()){

        const Unitig<U>* unitig = v_unitigs[pos_v_unitigs];

        const pair<size_t, size_t> lowpair = unitig->ccov.lowCoverageInfo();

        const size_t totalcoverage = unitig->coveragesum - lowpair.second;
        const size_t ccov_size = unitig->ccov.size();

        const string str = unitig->seq.toString();

        deleted = false;

        for (vector<pair<int,int>>::const_iterator sit = sp.begin(); sit != sp.end(); ++sit) { //Iterate over created split unitigs

            const size_t len = sit->second - sit->first;
            //const uint64_t cov_tmp = (totalcoverage * len) / (ccov_size - lowpair.first); // Split unitig coverage
            const uint64_t cov_tmp = len * CompressedCoverage::getFullCoverage();
            const string split_str = str.substr(sit->first, len + k_ - 1); // Split unitig sequence

            if (len == 1){

                if (addUnitig(split_str, v_kmers_sz, pos_v_unitigs, false)) h_kmers_ccov.find(Kmer(split_str.c_str()).rep())->ccov.setFull();
                else v_kmers[v_kmers_sz++].second.ccov.setFull(); // We don't care about the coverage per k-mer anymore}
            }
            else {

                addUnitig(split_str, nxt_pos_insert_v_unitigs, pos_v_unitigs, false);

                v_unitigs[nxt_pos_insert_v_unitigs]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                v_unitigs[nxt_pos_insert_v_unitigs]->coveragesum = cov_tmp;

                ++nxt_pos_insert_v_unitigs;
            }
        }
    }

    --nxt_pos_insert_v_unitigs; //Position of the last unitig in the vector which is not NULL

    if (pos_v_unitigs != nxt_pos_insert_v_unitigs){ // Do not proceed to swap if swap positions are the same

        swapUnitigs(false, pos_v_unitigs, nxt_pos_insert_v_unitigs); // Swap unitigs

        // If the swapped unitig, previously in position nxt_pos_insert, was a split unitig
        // created in this method, do not try to split it again
        if (nxt_pos_insert_v_unitigs >= v_unitigs_sz) ++pos_v_unitigs;
        else --v_unitigs_sz;
    }
    else --v_unitigs_sz;

    deleteUnitig_<true>(false, false, nxt_pos_insert_v_unitigs);

    return deleted;
}

template<typename U, typename G>
UnitigMap<U, G> CompactedDBG<U, G>::find(const Kmer& km, const preAllocMinHashIterator<RepHash>& it_min_h) {

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
        hmap_min_unitigs_t::const_iterator it = hmap_min_unitigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVal1();
            const uint8_t flag_v = it.getVal2();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id = v(i, flag_v) >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((v(i, flag_v) & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()){

                            return UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km == km_rep, this);
                        }
                    }

                    if ((v(i, flag_v) & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(&it_min.s[mhr.pos]).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (v(i, flag_v) & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos_tmp = v(i, flag_v) & MASK_CONTIG_POS;

                    if (isShort){

                        if (min_h_res.pos == unitig_id_pos_tmp){

                            if (v_kmers[unitig_id].first == km_rep){

                                return UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, true, this);
                            }
                        }
                        else if ((min_h_res.pos == diff - unitig_id_pos_tmp) && (v_kmers[unitig_id].first == km_rep)){

                            return UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, false, this);
                        }
                    }
                    else {

                        pos_match = unitig_id_pos_tmp - min_h_res.pos;
                        len = v_unitigs[unitig_id]->length() - k_;

                        if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km)){

                            return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                        }

                        pos_match = unitig_id_pos_tmp - diff + min_h_res.pos;

                        if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->seq.compareKmer(pos_match, k_, km_twin)){

                            return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return UnitigMap<U, G>();
}

// pre: Some k-mers in the unitigs might have a coverage which is less than CompressedCoverage::getFullCoverage()
//      Those k-mers must be deleted from the unitigs.
// post: All unitigs have a per k-mer coverage of CompressedCoverage::getFullCoverage(). The graph is not
//       necessarily compacted after calling this function.
template<typename U, typename G>
pair<size_t, size_t> CompactedDBG<U, G>::extractAllUnitigs() {

    size_t i;
    size_t split = 0, deleted = 0;
    size_t v_kmers_sz = v_kmers.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t nxt_pos_insert = v_unitigs.size();

    for (typename h_kmers_ccov_t::iterator it(h_kmers_ccov.begin()); it != h_kmers_ccov.end(); ++it) {

        if (!it->ccov.isFull()){

            deleteUnitig_<is_void<U>::value>(false, true, it.getHash());
            ++deleted;
        }
    }

    for (i = 0; i < v_kmers_sz;) {

        if (!v_kmers[i].second.ccov.isFull()) {

            --v_kmers_sz;

            if (i != v_kmers_sz) swapUnitigs(true, i, v_kmers_sz);

            deleteUnitig_<is_void<U>::value>(true, false, v_kmers_sz);

            ++deleted;
        }
        else ++i;
    }

    for (i = 0; i < v_unitigs_sz;) { // Iterate over unitigs created so far

        if (!v_unitigs[i]->ccov.isFull()) { //Coverage not full, unitig must be splitted

            vector<pair<int,int>> sp = v_unitigs[i]->ccov.splittingVector();

            if (extractUnitig_<is_void<U>::value>(i, nxt_pos_insert, v_unitigs_sz, v_kmers_sz, sp)) deleted++;
            else {

                ++split;
                sp.clear();
            }
        }
        else ++i;
    }

    if (nxt_pos_insert < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert);
    if (v_kmers_sz < v_kmers.size()) v_kmers.resize(v_kmers_sz);

    return {split, deleted};
}

// pre: Some k-mers in the unitigs might have a coverage which is less than CompressedCoverage::getFullCoverage().
//      Unitigs must be split at the positions matching those low-coverage k-mers (not deleted).
// post: All unitigs have a per k-mer coverage of CompressedCoverage::getFullCoverage(). The graph is not
//       necessarily compacted after calling this function.
template<typename U, typename G>
pair<size_t, size_t> CompactedDBG<U, G>::splitAllUnitigs() {

    pair<size_t, size_t> p = {0, 0};

    size_t v_kmers_sz = v_kmers.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t nxt_pos_insert = v_unitigs.size();

    const size_t cov_full = CompressedCoverage::getFullCoverage();

    for (size_t i = 0; i < v_unitigs_sz;) { // Iterate over unitigs created so far

        const CompressedCoverage& ccov = v_unitigs[i]->ccov;

        if (!ccov.isFull()) { //Coverage not full, unitig must be splitted

            size_t prev_split_pos = 0;

            vector<pair<int,int>> sp;

            for (size_t pos = 0; pos < ccov.size(); ++pos){

                if ((ccov.covAt(pos) != cov_full) && (pos != prev_split_pos)){

                    sp.push_back({prev_split_pos, pos});
                    ++(p.second);

                    prev_split_pos = pos;
                }
            }

            sp.push_back({prev_split_pos, ccov.size()});

            ++(p.second);
            ++(p.first);

            extractUnitig_<is_void<U>::value>(i, nxt_pos_insert, v_unitigs_sz, v_kmers_sz, sp);
        }
        else ++i;
    }

    if (nxt_pos_insert < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert);
    if (v_kmers_sz < v_kmers.size()) v_kmers.resize(v_kmers_sz);

    return p;
}

template<typename U, typename G>
pair<size_t, size_t> CompactedDBG<U, G>::getSplitInfoAllUnitigs() const {

    pair<size_t, size_t> p = {0, 0};

    const size_t cov_full = CompressedCoverage::getFullCoverage();

    for (size_t i = 0; i < v_unitigs.size(); ++i) { // Iterate over unitigs created so far

        const CompressedCoverage& ccov = v_unitigs[i]->ccov;

        if (!ccov.isFull()) { //Coverage not full, unitig must be splitted

            size_t prev_split_pos = 0;

            for (size_t pos = 0; pos < ccov.size(); ++pos){

                if ((ccov.covAt(pos) != cov_full) && (pos != prev_split_pos)){

                    ++(p.second);

                    prev_split_pos = pos;
                }
            }

            ++(p.first);
            ++(p.second);
        }
    }

    return p;
}

template<typename U, typename G>
void CompactedDBG<U, G>::createJoinHT(vector<Kmer>* v_joins, KmerHashTable<Kmer>& joins, const size_t nb_threads) const {

    const size_t v_unitigs_size = v_unitigs.size();
    const size_t v_kmers_size = v_kmers.size();

    const size_t chunk_size = 10000;

    if (v_joins == nullptr){

        for (typename h_kmers_ccov_t::const_iterator it_ccov(h_kmers_ccov.begin()); it_ccov != h_kmers_ccov.end(); ++it_ccov) {

            const Kmer tail(it_ccov.getKey());
            const Kmer head_twin(tail.twin());

            Kmer fw, bw;

            const const_UnitigMap<U, G> cm(it_ccov.getHash(), 0, 1, k_, false, true, true, this);

            if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail);
            if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin);
        }

        if (nb_threads == 1){

            for (size_t i = 0; i != v_kmers_size; ++i) {

                const Kmer tail(v_kmers[i].first);
                const Kmer head_twin(tail.twin());

                Kmer fw, bw;

                const const_UnitigMap<U, G> cm(i, 0, 1, k_, true, false, true, this);

                if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail);
                if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin);
            }

            for (size_t i = 0; i != v_unitigs_size; ++i) {

                const CompressedSequence& seq = v_unitigs[i]->seq;

                const Kmer head_twin(seq.getKmer(0).twin());
                const Kmer tail(seq.getKmer(seq.size() - k_));

                Kmer fw, bw;

                const const_UnitigMap<U, G> cm(i, 0, 1, seq.size(), false, false, true, this);

                if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail);
                if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin);
            }
        }
        else {

            vector<vector<pair<Kmer, Kmer>>> t_v_out(nb_threads);
            vector<std::atomic_flag> lcks(nb_threads);

            for (auto& lck : lcks) lck.clear();

            auto worker_v_kmers = [&joins, &lcks, &t_v_out, this](typename vector<pair<Kmer, CompressedCoverage_t<U>>>::const_iterator a,
                                                                  typename vector<pair<Kmer, CompressedCoverage_t<U>>>::const_iterator b,
                                                                  const size_t thread_id){

                while (lcks[thread_id].test_and_set(std::memory_order_acquire));

                for (size_t i = a - v_kmers.begin(), end = b - v_kmers.begin(); i != end; ++i) {

                    const Kmer tail(v_kmers[i].first);
                    const Kmer head_twin(tail.twin());

                    Kmer fw, bw;

                    const const_UnitigMap<U, G> cm(i, 0, 1, k_, true, false, true, this);

                    if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) t_v_out[thread_id].push_back(make_pair(fw.twin(), tail));
                    if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) t_v_out[thread_id].push_back(make_pair(bw.twin(), head_twin));
                }

                lcks[thread_id].clear(std::memory_order_release);

                if (t_v_out[thread_id].size() >= 1000){

                    for (auto& lck : lcks){ //Acquire all the locks for insertion

                        while (lck.test_and_set(std::memory_order_acquire));
                    }

                    for (const auto& p : t_v_out[thread_id]) joins.insert(p.first, p.second);

                    for (auto& lck : lcks) lck.clear(std::memory_order_release); //Acquire all the locks for insertion

                    t_v_out[thread_id].clear();
                }
            };

            auto worker_v_unitigs = [&joins, &lcks, &t_v_out, this](typename vector<Unitig<U>*>::const_iterator a,
                                                                    typename vector<Unitig<U>*>::const_iterator b,
                                                                    const size_t thread_id){

                while (lcks[thread_id].test_and_set(std::memory_order_acquire));

                for (size_t i = a - v_unitigs.begin(), end = b - v_unitigs.begin(); i != end; ++i) {

                    const CompressedSequence& seq = v_unitigs[i]->seq;

                    const Kmer head_twin(seq.getKmer(0).twin());
                    const Kmer tail(seq.getKmer(seq.size() - k_));

                    Kmer fw, bw;

                    const const_UnitigMap<U, G> cm(i, 0, 1, seq.size(), false, false, true, this);

                    if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) t_v_out[thread_id].push_back(make_pair(fw.twin(), tail));
                    if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) t_v_out[thread_id].push_back(make_pair(bw.twin(), head_twin));
                }

                lcks[thread_id].clear(std::memory_order_release);

                if (t_v_out[thread_id].size() >= 1000){

                    for (auto& lck : lcks){ //Acquire all the locks for insertion

                        while (lck.test_and_set(std::memory_order_acquire));
                    }

                    for (const auto& p : t_v_out[thread_id]) joins.insert(p.first, p.second);

                    for (auto& lck : lcks) lck.clear(std::memory_order_release); //Acquire all the locks for insertion

                    t_v_out[thread_id].clear();
                }
            };

            {

                auto it_kmer = v_kmers.begin();
                auto it_kmer_end = v_kmers.end();

                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mutex_it_km;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&, t]{

                            auto l_it_kmer = v_kmers.begin();
                            auto l_it_kmer_end = v_kmers.end();

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_it_km);

                                    if (it_kmer == it_kmer_end) return;

                                    l_it_kmer = it_kmer;
                                    l_it_kmer_end = it_kmer;

                                    if (distance(l_it_kmer, it_kmer_end) >= 1000){

                                        advance(l_it_kmer_end, 1000);

                                        it_kmer = l_it_kmer_end;
                                    }
                                    else {

                                        it_kmer = it_kmer_end;
                                        l_it_kmer_end = it_kmer_end;
                                    }
                                }

                                worker_v_kmers(l_it_kmer, l_it_kmer_end, t);
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();

                for (size_t t = 0; t < nb_threads; ++t){

                    for (const auto& p : t_v_out[t]) joins.insert(p.first, p.second);

                    t_v_out[t].clear();
                }
            }

            {
                auto it_unitig = v_unitigs.begin();
                auto it_unitig_end = v_unitigs.end();

                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mutex_it_unitig;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&, t]{

                            auto l_it_unitig = v_unitigs.begin();
                            auto l_it_unitig_end = v_unitigs.end();

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_it_unitig);

                                    if (it_unitig == it_unitig_end) return;

                                    l_it_unitig = it_unitig;
                                    l_it_unitig_end = it_unitig;

                                    if (distance(l_it_unitig, it_unitig_end) >= 1000){

                                        advance(l_it_unitig_end, 1000);

                                        it_unitig = l_it_unitig_end;
                                    }
                                    else {

                                        it_unitig = it_unitig_end;
                                        l_it_unitig_end = it_unitig_end;
                                    }
                                }

                                worker_v_unitigs(l_it_unitig, l_it_unitig_end, t);
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();

                for (size_t t = 0; t < nb_threads; ++t){

                    for (const auto& p : t_v_out[t]) joins.insert(p.first, p.second);

                    t_v_out[t].clear();
                }
            }
        }
    }
    else {

        Kmer fw;

        for (auto km : *v_joins){

            const const_UnitigMap<U, G> cm = find(km, true);

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
    }
}

// use:  joined = mapper.joinUnitigs()
// pre:  no short unitigs exist in sUnitigs.
// post: all unitigs that could be connected have been connected
//       joined is the number of joined unitigs
template<typename U, typename G>
template<bool is_void>
typename std::enable_if<!is_void, size_t>::type CompactedDBG<U, G>::joinUnitigs_(vector<Kmer>* v_joins, const size_t nb_threads) {

    size_t i;
    size_t joined = 0;
    size_t cov_full = CompressedCoverage::getFullCoverage();
    size_t v_unitigs_size = v_unitigs.size();
    size_t v_kmers_size = v_kmers.size();

    const size_t chunk_size = 10000;

    // a and b are candidates for joining
    KmerHashTable<Kmer> joins;

    createJoinHT(v_joins, joins, nb_threads);

    if (v_joins != nullptr) v_joins->clear();

    for (KmerHashTable<Kmer>::iterator it(joins.begin()); it != joins.end(); ++it) {

        const Kmer head(*it);
        const Kmer tail(it.getKey().twin());

        UnitigMap<U, G> cmHead(find(head, true));
        UnitigMap<U, G> cmTail(find(tail, true));

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
                string joinSeq;

                joinSeq.reserve((len_k_head ? k_ : cmHead.size) + (len_k_tail ? k_ : cmTail.size) - k_ + 1);

                if (headDir) joinSeq = len_k_head ? cmHead_head.toString() : v_unitigs[cmHead.pos_unitig]->seq.toString();
                else joinSeq = len_k_head ? cmHead_head.twin().toString() : v_unitigs[cmHead.pos_unitig]->seq.rev().toString();

                if (tailDir) joinSeq.append(len_k_tail ? cmTail_head.toString() : v_unitigs[cmTail.pos_unitig]->seq.toString(), k_ - 1, string::npos);
                else joinSeq.append(len_k_tail ? cmTail_head.twin().toString() : v_unitigs[cmTail.pos_unitig]->seq.rev().toString(), k_ - 1, string::npos);

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

                Unitig<U> data_tmp; //Store temporarily the new merged data
                Unitig<U>* unitig; //New unitig

                cmHead.strand = headDir;
                cmHead.dist = 0;
                cmHead.len = cmHead.size - k_ + 1;

                cmTail.strand = tailDir;
                cmTail.dist = 0;
                cmTail.len = cmTail.size - k_ + 1;

                data_tmp.getData()->concat(cmHead, cmTail);

                cmHead.getData()->clear(cmHead);
                cmTail.getData()->clear(cmTail);

                if (cmHead.isShort || cmHead.isAbundant){

                    if (cmHead.isShort){ //If head is a short unitig, swap and delete it

                        --v_kmers_size;

                        if (cmHead.pos_unitig != v_kmers_size){

                            swapUnitigs(true, cmHead.pos_unitig, v_kmers_size);

                            // If the last unitig of the vector used for the swap was the tail
                            if (cmTail.isShort && (v_kmers_size == cmTail.pos_unitig)) cmTail.pos_unitig = cmHead.pos_unitig;
                        }

                        deleteUnitig_<false>(true, false, v_kmers_size);
                    }
                    else if (cmHead.isAbundant) deleteUnitig_<false>(false, true, cmHead.pos_unitig);
                }

                if (cmTail.isShort || cmTail.isAbundant){

                    if (cmTail.isShort){ //If tail is a short unitig, swap and delete it

                        --v_kmers_size;

                        if (cmTail.pos_unitig != v_kmers_size){

                            swapUnitigs(true, cmTail.pos_unitig, v_kmers_size);

                            if (cmHead.isShort && (v_kmers_size == cmHead.pos_unitig)) cmHead.pos_unitig = cmTail.pos_unitig;
                        }

                        deleteUnitig_<false>(true, false, v_kmers_size);
                    }
                    else if (cmTail.isAbundant) deleteUnitig_<false>(false, true, cmTail.pos_unitig);
                }

                if (len_k_head && len_k_tail){

                    addUnitig(joinSeq, v_unitigs_size);
                    unitig = v_unitigs[v_unitigs_size];
                    ++v_unitigs_size;
                }
                else if (len_k_head){

                    deleteUnitig_<false>(false, false, cmTail.pos_unitig);
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

                        deleteUnitig_<false>(false, false, v_unitigs_size);
                    }

                    deleteUnitig_<false>(false, false, cmHead.pos_unitig);
                    addUnitig(joinSeq, cmHead.pos_unitig);
                    unitig = v_unitigs[cmHead.pos_unitig];
                }

                unitig->coveragesum = covsum;
                if (covsum >= cov_full * unitig->numKmers()) unitig->ccov.setFull();

                *(unitig->getData()) = std::move(*(data_tmp.getData()));

                ++joined;
            }
        }
    }

    if (v_unitigs_size < v_unitigs.size()) v_unitigs.resize(v_unitigs_size);
    if (v_kmers_size < v_kmers.size()) v_kmers.resize(v_kmers_size);

    return joined;
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<is_void, size_t>::type CompactedDBG<U, G>::joinUnitigs_(vector<Kmer>* v_joins, const size_t nb_threads) {

    size_t i;
    size_t joined = 0;
    size_t cov_full = CompressedCoverage::getFullCoverage();
    size_t v_unitigs_size = v_unitigs.size();
    size_t v_kmers_size = v_kmers.size();

    const size_t chunk_size = 10000;

    // a and b are candidates for joining
    KmerHashTable<Kmer> joins;

    createJoinHT(v_joins, joins, nb_threads);

    if (v_joins != nullptr) v_joins->clear();

    for (KmerHashTable<Kmer>::iterator it = joins.begin(); it != joins.end(); ++it) {

        const Kmer head = *it;
        const Kmer tail = it.getKey().twin();

        UnitigMap<U, G> cmHead = find(head, true);
        UnitigMap<U, G> cmTail = find(tail, true);

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
                string joinSeq;

                joinSeq.reserve((len_k_head ? k_ : cmHead.size) + (len_k_tail ? k_ : cmTail.size) - k_ + 1);

                if (headDir) joinSeq = len_k_head ? cmHead_head.toString() : v_unitigs[cmHead.pos_unitig]->seq.toString();
                else joinSeq = len_k_head ? cmHead_head.twin().toString() : v_unitigs[cmHead.pos_unitig]->seq.rev().toString();

                if (tailDir) joinSeq.append(len_k_tail ? cmTail_head.toString() : v_unitigs[cmTail.pos_unitig]->seq.toString(), k_ - 1, string::npos);
                else joinSeq.append(len_k_tail ? cmTail_head.twin().toString() : v_unitigs[cmTail.pos_unitig]->seq.rev().toString(), k_ - 1, string::npos);

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

                Unitig<U>* unitig; //New unitig

                cmTail.strand = tailDir;
                cmHead.strand = headDir;

                if (cmHead.isShort || cmHead.isAbundant){

                    if (cmHead.isShort){ //If head is a short unitig, swap and delete it

                        --v_kmers_size;

                        if (cmHead.pos_unitig != v_kmers_size){

                            swapUnitigs(true, cmHead.pos_unitig, v_kmers_size);

                            // If the last unitig of the vector used for the swap was the tail
                            if (cmTail.isShort && (v_kmers_size == cmTail.pos_unitig)) cmTail.pos_unitig = cmHead.pos_unitig;
                        }

                        deleteUnitig_<true>(true, false, v_kmers_size);
                    }
                    else if (cmHead.isAbundant) deleteUnitig_<true>(false, true, cmHead.pos_unitig);
                }

                if (cmTail.isShort || cmTail.isAbundant){

                    if (cmTail.isShort){ //If tail is a short unitig, swap and delete it

                        --v_kmers_size;

                        if (cmTail.pos_unitig != v_kmers_size){

                            swapUnitigs(true, cmTail.pos_unitig, v_kmers_size);

                            if (cmHead.isShort && (v_kmers_size == cmHead.pos_unitig)) cmHead.pos_unitig = cmTail.pos_unitig;
                        }

                        deleteUnitig_<true>(true, false, v_kmers_size);
                    }
                    else if (cmTail.isAbundant) deleteUnitig_<true>(false, true, cmTail.pos_unitig);
                }

                if (len_k_head && len_k_tail){

                    addUnitig(joinSeq, v_unitigs_size);
                    unitig = v_unitigs[v_unitigs_size];
                    ++v_unitigs_size;
                }
                else if (len_k_head){

                    deleteUnitig_<true>(false, false, cmTail.pos_unitig);
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

                        deleteUnitig_<true>(false, false, v_unitigs_size);
                    }

                    deleteUnitig_<true>(false, false, cmHead.pos_unitig);
                    addUnitig(joinSeq, cmHead.pos_unitig);
                    unitig = v_unitigs[cmHead.pos_unitig];
                }

                unitig->coveragesum = covsum;
                if (covsum >= cov_full * unitig->numKmers()) unitig->ccov.setFull();

                ++joined;
            }
        }
    }

    if (v_unitigs_size < v_unitigs.size()) v_unitigs.resize(v_unitigs_size);
    if (v_kmers_size < v_kmers.size()) v_kmers.resize(v_kmers_size);

    return joined;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::checkJoin(const Kmer& a, const const_UnitigMap<U, G>& cm_a, Kmer& b) const {

    size_t i, j, count_succ;

    vector<const_UnitigMap<U, G>> v_um(findSuccessors(a, 2, true));

    for (i = 0, count_succ = 0; i != 4; ++i){

        if (!v_um[i].isEmpty){ ++count_succ; j = i; }
    }

    if (count_succ == 1) {

        Kmer cand_head, ac_head;

        const Kmer fw_cand(a.forwardBase(alpha[j]));

        const const_UnitigMap<U, G> cm_cand(v_um[j]);

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

template<typename U, typename G>
void CompactedDBG<U, G>::check_fp_tips(KmerHashTable<bool>& ignored_km_tips){

    uint64_t nb_real_short_tips = 0;

    size_t nxt_pos_insert_v_unitigs = v_unitigs.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t v_kmers_sz = v_kmers.size();

    vector<pair<int,int>> sp;

    for (KmerHashTable<bool>::iterator it(ignored_km_tips.begin()); it != ignored_km_tips.end(); ++it) {

        Kmer km(it.getKey());

        UnitigMap<U, G> cm(find(km, true)); // Check if the (short) tip actually exists

        if (!cm.isEmpty){ // IF the tip exists

            ++nb_real_short_tips;

            bool not_found = true;

            for (size_t i = 0; (i < 4) && not_found; ++i) {

                UnitigMap<U, G> cm_bw(find(km.backwardBase(alpha[i])));

                if (!cm_bw.isEmpty && !cm_bw.isAbundant && !cm_bw.isShort){

                    if (cm_bw.strand) ++cm_bw.dist;

                    if ((cm_bw.dist != 0) && (cm_bw.dist != cm_bw.size - k_ + 1)){

                        sp.push_back(make_pair(0, cm_bw.dist));
                        sp.push_back(make_pair(cm_bw.dist, cm_bw.size - k_ + 1));

                        extractUnitig_<is_void<U>::value>(cm_bw.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

                        sp.clear();
                    }

                    not_found = false;
                }
            }

            for (size_t i = 0; (i < 4) && not_found; ++i) {

                UnitigMap<U, G> cm_fw(find(km.forwardBase(alpha[i])));

                if (!cm_fw.isEmpty && !cm_fw.isAbundant && !cm_fw.isShort){

                    if (!cm_fw.strand) ++cm_fw.dist;

                    if ((cm_fw.dist != 0) && (cm_fw.dist != cm_fw.size - k_ + 1)){

                        sp.push_back(make_pair(0, cm_fw.dist));
                        sp.push_back(make_pair(cm_fw.dist, cm_fw.size - k_ + 1));

                        extractUnitig_<is_void<U>::value>(cm_fw.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

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

template<typename U, typename G>
size_t CompactedDBG<U, G>::removeUnitigs(bool rmIsolated, bool clipTips, vector<Kmer>& v){

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

    Unitig<U>* unitig = nullptr;

    for (j = 0; j < v_unitigs_sz; ++j) {

        unitig = v_unitigs[j];

        if (unitig->numKmers() < k_){

            const Kmer head(unitig->seq.getKmer(0));

            nb_pred = 0;

            for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

                if (!find(head.backwardBase(alpha[i]), true).isEmpty){

                    ++nb_pred;
                    if (clipTips) km = head.backwardBase(alpha[i]);
                }
            }

            if (nb_pred <= lim){

                const Kmer tail(unitig->seq.getKmer(unitig->seq.size() - k_));

                nb_succ = 0;

                for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                    if (!find(tail.forwardBase(alpha[i]), true).isEmpty){

                        ++nb_succ;
                        if (clipTips) km = tail.forwardBase(alpha[i]);
                    }
                }

                if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))) { //Unitig is isolated

                    ++removed;
                    --v_unitigs_sz;

                    if (j != v_unitigs_sz){

                        swapUnitigs(false, j, v_unitigs_sz),
                        --j;
                    }

                    if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
                }
            }
        }
    }

    for (j = 0; j < v_kmers_sz; ++j) {

        const pair<Kmer, CompressedCoverage_t<U>>& p = v_kmers[j];

        nb_pred = 0;

        for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

            if (!find(p.first.backwardBase(alpha[i]), true).isEmpty){

                ++nb_pred;
                if (clipTips) km = p.first.backwardBase(alpha[i]);
            }
        }

        if (nb_pred <= lim){

            nb_succ = 0;

            for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                if (!find(p.first.forwardBase(alpha[i]), true).isEmpty){

                    ++nb_succ;
                    if (clipTips) km = p.first.forwardBase(alpha[i]);
                }
            }

            if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))) { //Unitig is isolated

                ++removed;
                --v_kmers_sz;

                if (j != v_kmers_sz){

                    swapUnitigs(true, j, v_kmers_sz),
                    --j;
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

                *it = CompressedCoverage_t<U>();

                if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
            }
        }
    }

    for (j = v_unitigs_sz; j < v_unitigs.size(); ++j) deleteUnitig_<is_void<U>::value>(false, false, j);
    v_unitigs.resize(v_unitigs_sz);

    for (j = v_kmers_sz; j < v_kmers.size(); ++j) deleteUnitig_<is_void<U>::value>(true, false, j);
    v_kmers.resize(v_kmers_sz);

    for (typename h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); ++it){

        if (it->ccov.size() == 0) deleteUnitig_<is_void<U>::value>(false, true, it.getHash());
    }

    return removed;
}

template<typename U, typename G>
void CompactedDBG<U, G>::writeFASTA(const string& graphfilename) const {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = v_kmers.size();
    const size_t graph_sz = size();

    size_t i = 0;

    ofstream graphfile;
    ostream graph(0);

    graphfile.open(graphfilename.c_str());
    graph.rdbuf(graphfile.rdbuf());
    graph.sync_with_stdio(false);

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

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<!is_void, void>::type CompactedDBG<U, G>::writeGFA_sequence_(GFA_Parser& graph, KmerHashTable<size_t>& idmap) const {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = v_kmers.size();

    size_t i, labelA, labelB, id = v_unitigs_sz + v_kmers_sz + 1;

    for (labelA = 1; labelA <= v_unitigs_sz; ++labelA) {

        const Unitig<U>* unitig = v_unitigs[labelA - 1];
        const string slabelA = std::to_string(labelA);
        const string data = unitig->getData()->serialize();

        graph.write_sequence(slabelA, unitig->seq.size(), unitig->seq.toString(), data == "" ? data : string("DA:Z:" + data));
    }

    for (labelA = 1; labelA <= v_kmers_sz; ++labelA) {

        const pair<Kmer, CompressedCoverage_t<U>>& p = v_kmers[labelA - 1];
        const string slabelA = std::to_string(labelA + v_unitigs_sz);
        const string data = p.second.getData()->serialize();

        graph.write_sequence(slabelA, k_, p.first.toString(), data == "" ? data : string("DA:Z:" + data));
    }

    for (typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); ++it) {

        const string slabelA = std::to_string(id);
        const string data = it->getData()->serialize();

        idmap.insert(it.getKey(), id);

        graph.write_sequence(slabelA, k_, it.getKey().toString(), data == "" ? data : string("DA:Z:" + data));

        ++id;
    }
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<is_void, void>::type CompactedDBG<U, G>::writeGFA_sequence_(GFA_Parser& graph, KmerHashTable<size_t>& idmap) const {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = v_kmers.size();

    size_t i, labelA, labelB, id = v_unitigs_sz + v_kmers_sz + 1;

    for (labelA = 1; labelA <= v_unitigs_sz; ++labelA) {

        const Unitig<U>* unitig = v_unitigs[labelA - 1];
        const string slabelA = std::to_string(labelA);

        graph.write_sequence(slabelA, unitig->seq.size(), unitig->seq.toString(), "");
    }

    for (labelA = 1; labelA <= v_kmers_sz; ++labelA) {

        const pair<Kmer, CompressedCoverage_t<U>>& p = v_kmers[labelA - 1];
        const string slabelA = std::to_string(labelA + v_unitigs_sz);

        graph.write_sequence(slabelA, k_, p.first.toString(), "");
    }

    for (typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); ++it) {

        const string slabelA = std::to_string(id);

        idmap.insert(it.getKey(), id);

        graph.write_sequence(slabelA, k_, it.getKey().toString(), "");

        ++id;
    }
}

template<typename U, typename G>
void CompactedDBG<U, G>::writeGFA(const string& graphfilename, const size_t nb_threads) const {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = v_kmers.size();

    size_t i, labelA, labelB, id = v_unitigs_sz + v_kmers_sz + 1;

    const string header_tag("BV:Z:" + string(BFG_VERSION) + "\t" + "KL:Z:" + to_string(k_) + "\t" + "ML:Z:" + to_string(g_));

    KmerHashTable<size_t> idmap(h_kmers_ccov.size());

    GFA_Parser graph(graphfilename);

    graph.open_write(1, header_tag);

    writeGFA_sequence_<is_void<U>::value>(graph, idmap);

    if (nb_threads == 1){

        for (labelA = 1; labelA <= v_unitigs_sz; labelA++) {

            const Unitig<U>* unitig = v_unitigs[labelA - 1];
            const Kmer head = unitig->seq.getKmer(0);

            for (i = 0; i < 4; ++i) {

                const Kmer b = head.backwardBase(alpha[i]);
                const const_UnitigMap<U, G> cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    graph.write_edge(slabelA, 0, k_-1, false,
                                     slabelB, cand.strand ? 0 : cand.size - k_ + 1, cand.strand ? k_-1 : cand.size, !cand.strand);
                }
            }

            const Kmer tail = unitig->seq.getKmer(unitig->seq.size() - k_);

            for (i = 0; i < 4; ++i) {

                const Kmer b = tail.forwardBase(alpha[i]);
                const const_UnitigMap<U, G> cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    graph.write_edge(slabelA, unitig->seq.size() - k_ + 1, unitig->seq.size(), true,
                                     slabelB, cand.strand ? 0 : cand.size - k_ + 1, cand.strand ? k_-1 : cand.size, cand.strand);
                }
            }
        }

        for (labelA = v_unitigs_sz + 1; labelA <= v_kmers_sz + v_unitigs_sz; labelA++) {

            const pair<Kmer, CompressedCoverage_t<U>>& p = v_kmers[labelA - v_unitigs_sz - 1];

            for (i = 0; i < 4; ++i) {

                const Kmer b = p.first.backwardBase(alpha[i]);
                const const_UnitigMap<U, G> cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz : 0);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    graph.write_edge(slabelA, 0, k_-1, false,
                                     slabelB, cand.strand ? 0 : cand.size - k_ + 1, cand.strand ? k_-1 : cand.size, !cand.strand);
                }
            }

            for (i = 0; i < 4; ++i) {

                const Kmer b = p.first.forwardBase(alpha[i]);
                const const_UnitigMap<U, G> cand = find(b, true);

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                    else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    graph.write_edge(slabelA, 0, k_-1, true,
                                     slabelB, cand.strand ? 0 : cand.size - k_ + 1, cand.strand ? k_-1 : cand.size, cand.strand);
                }
            }
        }
    }
    else {

        auto worker_v_unitigs = [v_unitigs_sz, &idmap, this](const size_t labelA_start, const size_t labelA_end,
                                                             vector<pair<pair<size_t, bool>, pair<size_t, bool>>>* v_out){

            // We need to deal with the tail of long unitigs
            for (size_t labelA = labelA_start; labelA < labelA_end; ++labelA) {

                const Unitig<U>* unitig = v_unitigs[labelA - 1];

                const Kmer head = unitig->seq.getKmer(0);
                const Kmer tail = unitig->seq.getKmer(unitig->seq.size() - k_);

                vector<const_UnitigMap<U, G>> v_um = findPredecessors(head, true);

                for (const auto& um : v_um) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ? *(idmap.find(um.getUnitigHead().rep())) : um.pos_unitig + 1 + (um.isShort ? v_unitigs_sz : 0));
                        v_out->push_back(make_pair(make_pair(labelA, false), make_pair(labelB, !um.strand)));
                    }
                }

                v_um.clear();

                v_um = findSuccessors(tail, 4, true);

                for (const auto& um : v_um) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ? *(idmap.find(um.getUnitigHead().rep())) : um.pos_unitig + 1 + (um.isShort ? v_unitigs_sz : 0));
                        v_out->push_back(make_pair(make_pair(labelA, true), make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        auto worker_v_kmers = [v_unitigs_sz, &idmap, this](const size_t labelA_start, const size_t labelA_end,
                                                             vector<pair<pair<size_t, bool>, pair<size_t, bool>>>* v_out){

            // We need to deal with the tail of long unitigs
            for (size_t labelA = labelA_start; labelA < labelA_end; ++labelA) {

                const pair<Kmer, CompressedCoverage_t<U>>& p = v_kmers[labelA - v_unitigs_sz - 1];

                vector<const_UnitigMap<U, G>> v_um = findPredecessors(p.first, true);

                for (const auto& um : v_um) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ? *(idmap.find(um.getUnitigHead().rep())) : um.pos_unitig + 1 + (um.isShort ? v_unitigs_sz : 0));
                        v_out->push_back(make_pair(make_pair(labelA, false), make_pair(labelB, !um.strand)));
                    }
                }

                v_um.clear();

                v_um = findSuccessors(p.first, 4, true);

                for (const auto& um : v_um) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ? *(idmap.find(um.getUnitigHead().rep())) : um.pos_unitig + 1 + (um.isShort ? v_unitigs_sz : 0));
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

                    const string slabelA = std::to_string(p.first.first);
                    const string slabelB = std::to_string(p.second.first);

                    graph.write_edge(slabelA, 0, k_-1, p.first.second, slabelB, 0, k_-1, p.second.second);
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

                    const string slabelA = std::to_string(p.first.first);
                    const string slabelB = std::to_string(p.second.first);

                    graph.write_edge(slabelA, 0, k_-1, p.first.second, slabelB, 0, k_-1, p.second.second);
                }

                v_out[i].clear();
            }
        }
    }

    for (KmerHashTable<size_t>::iterator it = idmap.begin(); it != idmap.end(); it++) {

        labelA = *it;

        for (i = 0; i < 4; ++i) {

            const Kmer b = it.getKey().backwardBase(alpha[i]);
            const const_UnitigMap<U, G> cand = find(b, true);

            if (!cand.isEmpty) {

                if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                const string slabelA = std::to_string(labelA);
                const string slabelB = std::to_string(labelB);

                graph.write_edge(slabelA, 0, k_-1, false,
                                 slabelB, cand.strand ? 0 : cand.size - k_ + 1, cand.strand ? k_-1 : cand.size, !cand.strand);
            }
        }

        for (i = 0; i < 4; ++i) {

            const Kmer b = it.getKey().forwardBase(alpha[i]);
            const const_UnitigMap<U, G> cand = find(b, true);

            if (!cand.isEmpty) {

                if (cand.isAbundant) labelB = *(idmap.find(b.rep()));
                else labelB = cand.pos_unitig + 1 + (cand.isShort ? v_unitigs_sz: 0);

                const string slabelA = std::to_string(labelA);
                const string slabelB = std::to_string(labelB);

                graph.write_edge(slabelA, 0, k_-1, true,
                                 slabelB, cand.strand ? 0 : cand.size - k_ + 1, cand.strand ? k_-1 : cand.size, cand.strand);
            }
        }
    }

    graph.close();
}

template<typename U, typename G>
void CompactedDBG<U, G>::readGFA(const string& graphfilename) {

    size_t graph_file_id = 0;

    bool new_file_opened = false;

    GFA_Parser graph(graphfilename);

    graph.open_read();

    GFA_Parser::GFA_line r = graph.read(graph_file_id, new_file_opened, true);

    while ((r.first != nullptr) || (r.second != nullptr)){

        if (r.first != nullptr) addUnitig(r.first->seq, (r.first->seq.length() == k_) ? v_kmers.size() : v_unitigs.size());

        r = graph.read(graph_file_id, new_file_opened, true);
    }
}

template<typename U, typename G>
void CompactedDBG<U, G>::readFASTA(const string& graphfilename) {

    size_t graph_file_id = 0;

    string seq;

    FastqFile ff(vector<string>(1, graphfilename));

    while (ff.read_next(seq, graph_file_id) != -1) addUnitig(seq, (seq.length() == k_) ? v_kmers.size() : v_unitigs.size());
}

template<typename U, typename G>
void CompactedDBG<U, G>::mapRead(const const_UnitigMap<U, G>& um) {

    if (um.isEmpty) return; // nothing maps, move on

    if (um.isShort) v_kmers[um.pos_unitig].second.ccov.cover(um.dist, um.dist + um.len - 1);
    else if (um.isAbundant) h_kmers_ccov.find(um.pos_unitig)->ccov.cover(um.dist, um.dist + um.len - 1);
    else v_unitigs[um.pos_unitig]->ccov.cover(um.dist, um.dist + um.len - 1);
}

template<typename U, typename G>
void CompactedDBG<U, G>::unmapRead(const const_UnitigMap<U, G>& um) {

    if (um.isEmpty) return; // nothing maps, move on

    if (um.isShort) v_kmers[um.pos_unitig].second.ccov.uncover(um.dist, um.dist + um.len - 1);
    else if (um.isAbundant) h_kmers_ccov.find(um.pos_unitig)->ccov.uncover(um.dist, um.dist + um.len - 1);
    else v_unitigs[um.pos_unitig]->ccov.uncover(um.dist, um.dist + um.len - 1);
}

template<typename U, typename G>
vector<Kmer> CompactedDBG<U, G>::extractMercyKmers(BlockedBloomFilter& bf_uniq_km, const size_t nb_threads, const bool verbose) {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = v_kmers.size();

    size_t i, j;

    char km_tmp[MAX_KMER_SIZE];

    KmerHashTable<uint8_t> tips;

    vector<Kmer> v_out;

    for (typename h_kmers_ccov_t::iterator it_ccov = h_kmers_ccov.begin(); it_ccov != h_kmers_ccov.end(); ++it_ccov) {

        const Kmer km = it_ccov.getKey().rep();

        vector<const_UnitigMap<U, G>> v_um = findPredecessors(km, true);

        for (i = 0; (i != 4) && v_um[i].isEmpty; ++i){}

        v_um = findSuccessors(km, true);

        for (j = 0; (j != 4) && v_um[j].isEmpty; ++j){}

        if ((i == 4) && (j == 4)) tips.insert(km, 3);
        else if (j == 4) tips.insert(km, 2);
        else if (i == 4) tips.insert(km, 1);
    }

    for (size_t it_v_km = 0; it_v_km != v_kmers_sz; ++it_v_km) {

        const Kmer km = v_kmers[it_v_km].first.rep();

        vector<const_UnitigMap<U, G>> v_um = findPredecessors(km, true);

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

        vector<const_UnitigMap<U, G>> v_um = findPredecessors(head, true);

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

template<typename U, typename G>
size_t CompactedDBG<U, G>::joinTips(string filename_MBBF_uniq_kmers, const size_t nb_threads, const bool verbose) {

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

    size_t nb_join = joinUnitigs_<is_void<U>::value>(&v_mercy_km, nb_threads);

    if (verbose) cout << "CompactedDBG<U, G>::joinTips(): " << nb_join << " unitigs have been joined using mercy k-mers" << endl;

    return nb_join;
}

template<typename U, typename G>
void CompactedDBG<U, G>::setKmerGmerLength(const int kmer_length, const int minimizer_length){

    invalid = false;

    if (kmer_length <= 0){

        cerr << "CompactedDBG::CompactedDBG(): Length k of k-mers cannot be less than or equal to 0" << endl;
        invalid = true;
    }

    if (kmer_length >= MAX_KMER_SIZE){

        cerr << "CompactedDBG::CompactedDBG(): Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        invalid = true;
    }

    if (minimizer_length <= 0){

        cerr << "CompactedDBG::CompactedDBG(): Length g of minimizers cannot be less than or equal to 0" << endl;
        invalid = true;
    }

    if (minimizer_length >= MAX_KMER_SIZE){

        cerr << "CompactedDBG::CompactedDBG(): Length g of minimizers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        invalid = true;
    }

    if (minimizer_length > kmer_length - 2){

        cerr << "CompactedDBG::CompactedDBG(): Length g of minimizers cannot exceed k - 2" << endl;
        invalid = true;
    }

    if (!invalid){

        Kmer::set_k(kmer_length);
        Minimizer::set_g(minimizer_length);

        k_ = kmer_length;
        g_ = minimizer_length;
    }
}

template<typename U, typename G>
void CompactedDBG<U, G>::print() const {

    cout << "CompactedDBG::print(): v_unitigs.size() = " << v_unitigs.size() << endl;
    cout << "CompactedDBG::print(): v_kmers.size() = " << v_kmers.size() << endl;
    cout << "CompactedDBG::print(): h_kmers_ccov.size() = " << h_kmers_ccov.size() << endl;
    cout << "CompactedDBG::print(): hmap_min_unitigs.size() = " << hmap_min_unitigs.size() << endl;
}

#endif
