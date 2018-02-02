#include "ColoredCDBG.hpp"

ColoredCDBG::ColoredCDBG(int kmer_length, int minimizer_length) :   CompactedDBG(kmer_length, minimizer_length),
                                                                    color_sets(nullptr), nb_color_sets(0), nb_seeds(0) {

    std::random_device rd; //Seed
    std::default_random_engine generator(rd()); //Random number generator
    std::uniform_int_distribution<long long unsigned> distribution(0,0xFFFFFFFFFFFFFFFF); //Distribution on which to apply the generator

    //Initialize the hash function seeds for
    for (int i = 0; i < 256; ++i) seeds[i] = distribution(generator);
}

ColoredCDBG::~ColoredCDBG() {

    empty();
}

void ColoredCDBG::clear(){

    CompactedDBG::clear();

    nb_color_sets = 0;
    nb_seeds = 0;

    empty();
}

void ColoredCDBG::empty(){

    if (color_sets != nullptr){

        delete[] color_sets;
        color_sets = nullptr;
    }

    CompactedDBG::empty();
}

bool ColoredCDBG::build(const CCDBG_Build_opt& opt){

    CDBG_Build_opt opt_ = opt.getCDBG_Build_opt();

    return CompactedDBG::build(opt_);
}

bool ColoredCDBG::mapColors(const CCDBG_Build_opt& opt){

    if (opt.filename_seq_in.size() == 0){

        initColorSets(opt);
        buildColorSets(opt);
    }
    else readColorSets(opt);

    return true;
}

void ColoredCDBG::initColorSets(const CCDBG_Build_opt& opt, const size_t max_nb_hash){

    const size_t nb_locks = opt.nb_threads * 256;

    mutex mutex_km_overflow;

    std::atomic_flag* cs_locks = new std::atomic_flag[nb_locks];

    for (size_t i = 0; i < nb_locks; ++i) cs_locks[i].clear();

    nb_seeds = max_nb_hash;
    nb_color_sets = size();

    color_sets = new ColorSet[nb_color_sets];

    auto worker_function = [&](ColoredCDBG::iterator a, ColoredCDBG::iterator b){

        int i;

        uint64_t h_v, id_lock;

        for (auto& unitig = a; unitig != b; ++unitig){

            const Kmer head = unitig->getHead();

            for (i = 0; i < max_nb_hash; ++i){

                h_v = head.hash(seeds[i]) % nb_color_sets; // Hash to which we can possibly put our colorset for current kmer
                id_lock = h_v % nb_locks; // Lock ID for this hash

                cs_locks[id_lock].test_and_set(std::memory_order_acquire); // Set the corresponding lock

                if (color_sets[h_v].isUnoccupied()) break; // If color set is unoccupied, we want to use it, maintain lock

                cs_locks[id_lock].clear(std::memory_order_release); // Else, lock is released
            }

            if (i == max_nb_hash){ // IF we couldn't find a hash matching an unoccupied color set for current k-mer

                unique_lock<mutex> lock(mutex_km_overflow); // Set lock for insertion into hash table of k-mer overflow

                km_overflow.insert(head, ColorSet()); // Insertion
            }
            else {

                color_sets[h_v].setOccupied(); // Set color set to occupied
                cs_locks[id_lock].clear(std::memory_order_release); // Release lock
            }

            const HashID hid(static_cast<uint8_t>(i == max_nb_hash ? 0 : i + 1));
            unitig->setData(&hid); // Set the hash ID for this unitig
        }
    };

    {
        const size_t chunk = 1000;

        vector<thread> workers; // need to keep track of threads so we can join them

        ColoredCDBG::iterator g_a = begin();
        const ColoredCDBG::iterator g_b = end();

        mutex mutex_it;

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    ColoredCDBG::iterator l_a, l_b;

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

    delete[] cs_locks;

    cout << "Number of unitigs not hashed is " << km_overflow.size() << " on " << nb_color_sets << " unitigs." << endl;
}

void ColoredCDBG::buildColorSets(const CCDBG_Build_opt& opt){

    const size_t nb_locks = opt.nb_threads * 256;

    const int k_ = getK();
    const int chunk = 1000;

    size_t prev_file_id = 0;

    bool next_file = true;

    FileParser fp(opt.filename_seq_in);

    mutex mutex_km_overflow;

    std::atomic_flag* cs_locks = new std::atomic_flag[nb_locks];

    for (size_t i = 0; i < nb_locks; ++i) cs_locks[i].clear();

    // Main worker thread
    auto worker_function = [&](const vector<pair<string, size_t>>& v_read_color) {

        // for each input
        for (const auto& read_color : v_read_color) {

            for (KmerIterator it_km(read_color.first.c_str()), it_km_end; it_km != it_km_end; ++it_km) {

                UnitigMap<HashID> um = find(it_km->first);

                if (!um.isEmpty) {

                    if (um.strand || (um.dist != 0)){

                        um.len += um.lcp(read_color.first.c_str(), it_km->second + k_, um.strand ? um.dist + k_ : um.dist - 1, um.strand);
                        it_km += um.len - 1;
                    }

                    HashID* hid = um.getData();

                    if (hid->get() == 0){

                        unique_lock<mutex> lock(mutex_km_overflow);

                        setColor(um, read_color.second);
                    }
                    else {

                        const uint64_t id_lock = getHash(um) % nb_locks;

                        cs_locks[id_lock].test_and_set(std::memory_order_acquire); // Set the corresponding lock

                        setColor(um, read_color.second);

                        cs_locks[id_lock].clear(std::memory_order_release);
                    }
                }
            }
        }
    };

    auto reading_function = [&](vector<pair<string, size_t>>& v_read_color) {

        string s;

        size_t file_id = prev_file_id;
        size_t reads_now = 0;

        while (reads_now < chunk) {

            if (fp.read(s, file_id)) {

                v_read_color.emplace_back(make_pair(s, file_id));

                ++reads_now;
            }
            else {

                next_file = false;
                return true;
            }
        }

        next_file = true;

        if (file_id != prev_file_id){

            prev_file_id = file_id;
            return true;
        }

        prev_file_id = file_id;
        return false;
    };

    {
        bool stop = false;

        vector<thread> workers; // need to keep track of threads so we can join them
        vector<vector<pair<string, size_t>>> reads_colors(opt.nb_threads);

        mutex mutex_file;

        while (next_file){

            stop = false;

            for (size_t t = 0; t < opt.nb_threads; ++t){

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

            for (size_t t = 0; t < opt.nb_threads; ++t) reads_colors[t].clear();

            for (size_t i = 0; i < nb_color_sets; ++i) color_sets[i].optimize();

            for (KmerHashTable<ColorSet>::iterator it = km_overflow.begin(), it_end = km_overflow.end(); it != it_end; ++it) it->optimize();
        }
    }

    for (size_t i = 0; i < nb_color_sets; ++i) color_sets[i].optimize();

    for (KmerHashTable<ColorSet>::iterator it = km_overflow.begin(), it_end = km_overflow.end(); it != it_end; ++it) it->optimize();

    fp.close();

    delete[] cs_locks;
}

bool ColoredCDBG::setColor(const UnitigMap<HashID>& um, const size_t color_id) {

    if (!um.isEmpty && (color_sets != nullptr)){

        ColorSet* color_set = getColorSet(um);

        if (color_set != nullptr){

            color_set->add(um, color_id);

            return true;
        }
    }

    return false;
}

bool ColoredCDBG::joinColors(const UnitigMap<HashID>& um_dest, const UnitigMap<HashID>& um_src) {

    if (!um_dest.isEmpty && !um_src.isEmpty && (color_sets != nullptr)){

        ColorSet* color_set_dest = getColorSet(um_dest);
        ColorSet* color_set_src = getColorSet(um_src);

        if ((color_set_dest != nullptr) && (color_set_src != nullptr)){

            ColorSet new_cs, csd_rev, css_rev;

            size_t prev_color_id = 0xffffffffffffffff;
            size_t prev_km_dist = 0xffffffffffffffff;

            const size_t um_dest_km_sz = um_dest.size - getK() + 1;
            const size_t um_src_km_sz = um_src.size - getK() + 1;

            UnitigMap<HashID> new_um_dest(um_dest.pos_unitig, 0, 0, um_dest.size + um_src.size,
                                          um_dest.isShort, um_dest.isAbundant, um_dest.strand, *(um_dest.cdbg));

            if (!um_dest.strand){

                csd_rev = color_set_dest->reverse(um_dest.size);
                color_set_dest = &csd_rev;
            }

            ColorSet::const_iterator it = color_set_dest->begin(), it_end = color_set_dest->end();

            if (it != it_end){

                prev_color_id = *it / um_dest_km_sz;
                prev_km_dist = *it - (prev_color_id * um_dest_km_sz);

                new_um_dest.dist = prev_km_dist;
                new_um_dest.len = 1;

                ++it;
            }

            // Insert colors layer by layer
            for (; it != it_end; ++it){

                const size_t color_id = *it / um_dest_km_sz;
                const size_t km_dist = *it - (color_id * um_dest_km_sz);

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

            UnitigMap<HashID> new_um_src(um_src.pos_unitig, 0, 0, um_dest.size + um_src.size,
                                         um_src.isShort, um_src.isAbundant, um_src.strand, *(um_src.cdbg));

            if (!um_src.strand){

                css_rev = color_set_src->reverse(um_src.size);
                color_set_src = &css_rev;
            }

            it = color_set_src->begin(), it_end = color_set_src->end();

            if (it != it_end){

                prev_color_id = *it / um_src_km_sz;
                prev_km_dist = *it - (prev_color_id * um_src_km_sz);

                new_um_src.dist = prev_km_dist + um_dest.size;
                new_um_src.len = 1;

                ++it;
            }

            // Insert colors layer by layer
            for (; it != it_end; ++it){

                const size_t color_id = *it / um_src_km_sz;
                const size_t km_dist = *it - (color_id * um_src_km_sz);

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

ColorSet ColoredCDBG::extractColors(const UnitigMap<HashID>& um) const {

    ColorSet new_cs;

    if (!um.isEmpty && (color_sets != nullptr)){

        const ColorSet* cs = getColorSet(um);

        if (cs != nullptr){

            const size_t end = um.dist + um.len;
            const size_t um_km_sz = um.size - getK() + 1;

            UnitigMap<HashID> fake_um(0, 0, 1, um.len, false, false, um.strand, *(um.cdbg));

            ColorSet::const_iterator it = cs->begin(), it_end = cs->end();

            for (ColorSet::const_iterator it = cs->begin(), it_end = cs->end(); it != it_end; ++it){

                const size_t color_id = *it / um_km_sz;
                const size_t km_dist = *it - (color_id * um_km_sz);

                if ((km_dist >= um.dist) && (km_dist < end)){

                    fake_um.dist = km_dist - um.dist;

                    new_cs.add(fake_um, color_id);
                }
            }
        }
    }

    return new_cs;
}

ColorSet* ColoredCDBG::getColorSet(const UnitigMap<HashID>& um) {

    ColorSet* color_set = nullptr;

    if (!um.isEmpty && (color_sets != nullptr)){

        const Kmer head = um.getHead();

        const uint8_t hash_id = um.getData()->get();

        if (hash_id == 0){

            KmerHashTable<ColorSet>::iterator it = km_overflow.find(head);

            if (it != km_overflow.end()) color_set = &(*it);
        }
        else color_set = &color_sets[head.hash(seeds[hash_id - 1]) % nb_color_sets];
    }

    return color_set;
}

const ColorSet* ColoredCDBG::getColorSet(const UnitigMap<HashID>& um) const {

    if (!um.isEmpty && (color_sets != nullptr)){

        const Kmer head = um.getHead();

        const uint8_t hash_id = um.getData()->get();

        if (hash_id == 0){

            KmerHashTable<ColorSet>::const_iterator it = km_overflow.find(head);

            if (it != km_overflow.end()) return &(*it);
        }
        else return &color_sets[head.hash(seeds[hash_id - 1]) % nb_color_sets];
    }

    return nullptr;
}

bool ColoredCDBG::write(const string output_filename, const size_t nb_threads, const bool verbose){

    if (CompactedDBG::write(output_filename, nb_threads, true, verbose)){

        if (verbose) cout << endl << "ColoredCDBG::write(): Writing colors to disk" << endl;

        const string out = output_filename + ".bfg_colors";

        FILE* fp = fopen(out.c_str(), "wb");

        if (fp == NULL) {

            cerr << "ColoredCDBG::write(): Could not open file " << out << " for writing color sets" << endl;
            return false;
        }
        else {

            fclose(fp);

            if (std::remove(out.c_str()) != 0) cerr << "ColoredCDBG::write(): Could not remove temporary file " << out << endl;
        }

        ofstream colorsfile_out;
        ostream colors_out(nullptr);

        colorsfile_out.open(out.c_str(), ios_base::out | ios_base::binary);
        colors_out.rdbuf(colorsfile_out.rdbuf());

        const size_t format_version = BFG_COLOREDCDBG_FORMAT_VERSION;
        const size_t k_ = getK();
        const size_t km_overflow_sz = km_overflow.size();

        //Write the file format version number
        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&format_version), sizeof(size_t));
        //Write k-mer length
        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&k_), sizeof(size_t));
        //Write number of different seeds for hash function
        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_seeds), sizeof(size_t));
        //Write number of color sets in the graph
        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_color_sets), sizeof(size_t));
        //Write number of (kmer, color set) overflowing
        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&km_overflow_sz), sizeof(size_t));

        for (size_t i = 0; (i < nb_seeds) && colors_out.good(); ++i){
            //Write the hash function seeds of the graph
            colors_out.write(reinterpret_cast<const char*>(&seeds[i]), sizeof(uint64_t));
        }

        for (size_t i = 0; (i < nb_color_sets) && colors_out.good(); ++i) color_sets[i].write(colors_out); //Write the color sets

        KmerHashTable<ColorSet>::const_iterator it = km_overflow.begin(), it_end = km_overflow.end();

        for (; (it != it_end) && colors_out.good(); ++it){

            it.getKey().write(colors_out);
            it->write(colors_out);
        }

        const bool ret = colors_out.good();

        colorsfile_out.close();

        return ret;
    }

    return false;
}

bool ColoredCDBG::readColorSets(const CCDBG_Build_opt& opt){

    if (opt.verbose) cout << endl << "ColoredCDBG::readColorSets(): Reading color sets from disk" << endl;

    if (opt.filename_colors_in.size() == 0){

        cerr << "ColoredCDBG::readColorSets(): No color sets file given in input." << endl;
        return false;
    }

    if (opt.filename_colors_in.size() > 1){

        cerr << "ColoredCDBG::readColorSets(): Bifrost cannot use multiple color sets files in input at the moment." << endl;
        return false;
    }

    FILE* fp = fopen(opt.filename_colors_in[0].c_str(), "rb");

    if (fp == NULL) {

        cerr << "ColoredCDBG::readColorSets(): Could not open file " << opt.filename_colors_in[0] << " for reading color sets" << endl;
        return false;
    }
    else fclose(fp);

    ifstream colorsfile_in;
    istream colors_in(nullptr);

    colorsfile_in.open(opt.filename_colors_in[0].c_str(), ios_base::in | ios_base::binary);
    colors_in.rdbuf(colorsfile_in.rdbuf());

    size_t format_version , k_, km_overflow_sz;

    //Write the file format version number
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&format_version), sizeof(size_t));
    //Write k-mer length
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&k_), sizeof(size_t));
    //Write number of different seeds for hash function
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_seeds), sizeof(size_t));
    //Write number of color sets in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_color_sets), sizeof(size_t));
    //Write number of (kmer, color set) overflowing
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&km_overflow_sz), sizeof(size_t));

    if (k_ != Kmer::k){

        cerr << "ColoredCDBG::readColorSets(): Length k is not the same as length of k-mers used to build the graph" << endl;
        cerr << "ColoredCDBG::readColorSets(): It is possible the given color set file is not the one corresponding to the given graph file" << endl;
        return false;
    }

    if (nb_seeds >= 256){

        cerr << "ColoredCDBG::readColorSets(): Does not support more than 255 hash seeds" << endl;
        return false;
    }

    if (nb_color_sets != size()){

        cerr << "ColoredCDBG::readColorSets(): Number of color sets is not the same as the number of unitigs in the graph" << endl;
        return false;
    }

    if (color_sets != nullptr) delete[] color_sets;

    km_overflow = KmerHashTable<ColorSet>(km_overflow_sz);

    color_sets = new ColorSet[nb_color_sets];

    for (size_t i = 0; (i < nb_seeds) && colors_in.good(); ++i){
        //Write the hash function seeds of the graph
        colors_in.read(reinterpret_cast<char*>(&seeds[i]), sizeof(uint64_t));
    }

    for (size_t i = 0; (i < nb_color_sets) && colors_in.good(); ++i) color_sets[i].read(colors_in);

    Kmer km;

    for (size_t i = 0; (i < km_overflow_sz) && colors_in.good(); ++i){

        km.read(colors_in);

        std::pair<KmerHashTable<ColorSet>::iterator, bool> p = km_overflow.insert(km, ColorSet());

        p.first->read(colors_in);
    }

    const bool ret = colors_in.good();

    colorsfile_in.close();

    return true;
}

uint64_t ColoredCDBG::getHash(const UnitigMap<HashID>& um) const {

    if (!um.isEmpty && (color_sets != nullptr)){

        const Kmer head = um.getHead();
        const uint8_t hash_id = um.getData()->get();

        if (hash_id != 0) return head.hash(seeds[hash_id - 1]) % nb_color_sets;
    }

    return 0;
}

void ColoredCDBG::checkColors(const CCDBG_Build_opt& opt) {

    cout << "ColoredCDBG::checkColors(): Start" << endl;

    size_t file_id = 0;

    string s;

    KmerHashTable<tiny_vector<size_t, 1>> km_h;

    FastqFile FQ(opt.filename_seq_in);

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

    file_id = 0;

    FastqFile FQ2(opt.filename_seq_in);

    while (FQ2.read_next(s, file_id) >= 0){

        for (KmerIterator it_km(s.c_str()), it_km_end; it_km != it_km_end; ++it_km) {

            UnitigMap<HashID> um = find(it_km->first);

            if (um.isEmpty){

                cerr << "ColoredCDBG::checkColors(): K-mer " << it_km->first.toString() << " is not found in the graph" << endl;
                exit(1);
            }

            ColorSet* cs = getColorSet(um);

            if (cs == nullptr){

                cerr << "ColoredCDBG::checkColors(): K-mer " << it_km->first.toString() << " has no color set associated" << endl;
                exit(1);
            }

            KmerHashTable<tiny_vector<size_t, 1>>::const_iterator it = km_h.find(it_km->first.rep());

            if (it == km_h.end()){

                cerr << "ColoredCDBG::checkColors(): K-mer " << it_km->first.toString() << " was not inserted in the hash table of k-mers" << endl;
                exit(1);
            }

            const tiny_vector<size_t, 1>& tv = *it;
            const size_t tv_nb_max_elem = tv.size() * 64;

            for (size_t i = 0; i < std::min(opt.filename_seq_in.size(), tv_nb_max_elem); ++i){

                const bool color_pres_graph = cs->contains(um, i);
                const bool color_pres_hasht = ((tv[i/64] >> (i%64)) & 0x1) == 0x1;

                if (color_pres_graph != color_pres_hasht){

                    cerr << "ColoredCDBG::checkColors(): Current color is " << file_id << ": " << opt.filename_seq_in[file_id] << endl;
                    cerr << "ColoredCDBG::checkColors(): K-mer " << it_km->first.toString() << " for color " << i << ": " << opt.filename_seq_in[i] << endl;
                    cerr << "ColoredCDBG::checkColors(): Full unitig: " << um.toString() << endl;
                    cerr << "ColoredCDBG::checkColors(): Present in graph: " << color_pres_graph << endl;
                    cerr << "ColoredCDBG::checkColors(): Present in hash table: " << color_pres_hasht << endl;

                    exit(1);
                }
            }
        }
    }

    FQ2.close();

    cout << "ColoredCDBG::checkColors(): Checked all colors of all k-mers: everything is fine" << endl;
    cout << "ColoredCDBG::checkColors(): Number of k-mers in the graph: " << km_h.size() << endl;
}
