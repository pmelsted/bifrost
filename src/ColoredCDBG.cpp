#include "ColoredCDBG.hpp"

/** Colored and Compacted de Bruijn graph constructor (set up an empty colored cdBG).
* @param kmer_length is the length k of k-mers used in the graph (each unitig is of length at least k).
* @param minimizer_length is the length g of minimizers (g < k) used in the graph.
*/
ColoredCDBG::ColoredCDBG(int kmer_length, int minimizer_length) :   CompactedDBG(kmer_length, minimizer_length),
                                                                    color_sets(nullptr), invalid(false), nb_seeds(0),
                                                                    nb_color_sets(0) {

    std::random_device rd; //Seed
    std::default_random_engine generator(rd()); //Random number generator
    std::uniform_int_distribution<long long unsigned> distribution(0,0xFFFFFFFFFFFFFFFF); //Distribution on which to apply the generator

    //Initialize the hash function seeds for
    for (int i = 0; i < 256; ++i) seeds[i] = distribution(generator);

    invalid = CompactedDBG::isInvalid();
}

/** Colored and compacted de Bruijn graph copy constructor (copy a colored cdBG).
* This function is expensive in terms of time and memory as the content of a colored and compacted
* de Bruijn graph is copied. After the call to this function, the same graph exists twice in memory.
* @param o is a constant reference to the colored and compacted de Bruijn graph to copy.
*/
ColoredCDBG::ColoredCDBG(const ColoredCDBG& o) :    CompactedDBG(o), invalid(o.invalid), nb_seeds(o.nb_seeds),
                                                    nb_color_sets(o.nb_color_sets), color_names(o.color_names),
                                                    color_sets(nullptr) {

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    if ((o.color_sets != nullptr) && (o.nb_color_sets != 0)){

        color_sets = new ColorSet[nb_color_sets];

        copy(o.color_sets, o.color_sets + nb_color_sets, color_sets);
    }
}

/** Colored and compacted de Bruijn graph move constructor (move a colored cdBG).
* The content of o is moved ("transfered") to a new colored and compacted de Bruijn graph.
* The colored and compacted de Bruijn graph referenced by o will be empty after the call to this constructor.
* @param o is a reference on a reference to the colored and compacted de Bruijn graph to move.
*/
ColoredCDBG::ColoredCDBG(ColoredCDBG&& o) : CompactedDBG(o), invalid(o.invalid), nb_seeds(o.nb_seeds),
                                            nb_color_sets(o.nb_color_sets), color_names(move(o.color_names)),
                                            color_sets(o.color_sets) {

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    o.color_sets = nullptr;
    o.clear();
}

/** Colored and compacted de Bruijn graph destructor.
*/
ColoredCDBG::~ColoredCDBG() {

    empty();
}

/** Clear the graph: empty the graph and reset its parameters.
*/
void ColoredCDBG::clear(){

    invalid = true;

    nb_color_sets = 0;
    nb_seeds = 0;

    CompactedDBG::clear();

    empty();
}

/** Empty the graph (does not reset its parameters).
*/
void ColoredCDBG::empty(){

    if (color_sets != nullptr){

        delete[] color_sets;
        color_sets = nullptr;
    }

    color_names.clear();

    CompactedDBG::empty();
}

/** Colored and compacted de Bruijn graph copy assignment operator (copy a colored cdBG).
* This function is expensive in terms of time and memory as the content of a colored and compacted
* de Bruijn graph is copied.  After the call to this function, the same graph exists twice in memory.
* @param o is a constant reference to the colored and compacted de Bruijn graph to copy.
* @return a reference to the colored and compacted de Bruijn which is the copy.
*/
ColoredCDBG& ColoredCDBG::operator=(const ColoredCDBG& o) {

    CompactedDBG::operator=(o);

    invalid = o.invalid;
    nb_seeds = o.nb_seeds;
    nb_color_sets = o.nb_color_sets;

    color_names = o.color_names;
    color_sets = nullptr;

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    if ((o.color_sets != nullptr) && (o.nb_color_sets != 0)){

        color_sets = new ColorSet[nb_color_sets];

        copy(o.color_sets, o.color_sets + nb_color_sets, color_sets);
    }

    return *this;
}

/** Colored and compacted de Bruijn graph move assignment operator (move a colored cdBG).
* The content of o is moved ("transfered") to a new colored and compacted de Bruijn graph.
* The colored and compacted de Bruijn graph referenced by o will be empty after the call to this operator.
* @param o is a reference on a reference to the colored and compacted de Bruijn graph to move.
* @return a reference to the colored and compacted de Bruijn which has (and owns) the content of o.
*/
ColoredCDBG& ColoredCDBG::operator=(ColoredCDBG&& o) {

    if (this != &o) {

        CompactedDBG::operator=(o);

        invalid = o.invalid;
        nb_seeds = o.nb_seeds;
        nb_color_sets = o.nb_color_sets;

        color_names = move(o.color_names);
        color_sets = o.color_sets;

        memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

        o.color_sets = nullptr;
        o.clear();
    }

    return *this;
}

/** Build the Colored and compacted de Bruijn graph (only the unitigs).
* A call to ColoredCDBG::mapColors is required afterwards to map colors to unitigs.
* @param opt is a structure from which the members are parameters of this function. See CCDBG_Build_opt.
* @return boolean indicating if the graph has been built successfully.
*/
bool ColoredCDBG::build(const CCDBG_Build_opt& opt){

    if (!invalid){

        CDBG_Build_opt opt_ = opt.getCDBG_Build_opt();

        invalid = !CompactedDBG::build(opt_);
    }
    else cerr << "ColoredCDBG::build(): Graph is invalid and cannot be built." << endl;

    return !invalid;
}

/** Map the colors to the unitigs. This is done by reading the input files and querying the graph.
* If a color filename is provided in opt.filename_colors_in, colors are loaded from that file instead.
* @param opt is a structure from which the members are parameters of this function. See CCDBG_Build_opt.
* @return boolean indicating if the colors have been mapped successfully.
*/
bool ColoredCDBG::mapColors(const CCDBG_Build_opt& opt){

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

/** Set a color for a unitig or a sub-unitig.
* @param um is a UnitigMap representing the mapping of a unitig for which the color must be added.
* The color will be added only for the sub-unitig mapped, i.e, unitig[um.dist..um.dist+um.len+k-1]
* @param color_id is the ID of the color to add.
* @return boolean indicating if the color was successfully set.
*/
bool ColoredCDBG::setColor(const UnitigMap<HashID>& um, const size_t color_id) {

    if (!invalid){

        if (color_id < color_names.size()){

            if (!um.isEmpty && (color_sets != nullptr)){

                ColorSet* color_set = getColorSet(um);

                if (color_set != nullptr){

                    color_set->add(um, color_id);

                    return true;
                }
            }
        }
        else cerr << "ColoredCDBG::setColor(): Color ID does not match any color inserted." << endl;
    }
    else cerr << "ColoredCDBG::setColor(): Graph is invalid, it is not possible to set a color for a unitig." << endl;

    return false;
}

/** Join two color sets (union). All colors of the color set matching the unitig mapped by um_src
* are added to the color set matching the unitig mapped  by um_dest. The reverse-complements of
* um_src and um_dest (UnitigMap<HashID>::strand) are considered. Any sub-unitig information such
* as UnitigMap<HashID>::dist or UnitigMap<HashID>::len is discarded.
* @param um_dest is a UnitigMap representing a mapping to a unitig. The unitig color set will
* contain the union of itself and the color set matching um_src.
* @param um_src is a UnitigMap representing a mapping to a unitig. The unitig color set will be
* joined with the color set matching um_dest.
* @return boolean indicating if the color sets have been joined successfully.
*/
bool ColoredCDBG::joinColors(const UnitigMap<HashID>& um_dest, const UnitigMap<HashID>& um_src) {

    if (!invalid){

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
    }
    else cerr << "ColoredCDBG::joinColors(): Graph is invalid, it is not possible to join color sets." << endl;

    return false;
}

/** Extract the color set matching a sub-unitig (see UnitigMap).
* @param um is a UnitigMap representing a mapping to a unitig from the graph.
* @return a new color set
*/
ColorSet ColoredCDBG::extractColors(const UnitigMap<HashID>& um) const {

    ColorSet new_cs;

    if (!invalid){

        if (!um.isEmpty && (color_sets != nullptr)){

            const ColorSet* cs = getColorSet(um);

            if (cs != nullptr){

                const size_t end = um.dist + um.len;
                const size_t um_km_sz = um.size - getK() + 1;

                UnitigMap<HashID> fake_um(0, 0, 1, um.len, false, false, um.strand, *(um.cdbg));

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
    }
    else cerr << "ColoredCDBG::extractColors(): Graph is invalid, no colors can be extracted." << endl;

    return new_cs;
}

/** Same as ColoredCDBG::extractColors but extract the color names instead.
* @param um is a UnitigMap representing a mapping to a unitig from the graph.
* @return a vector a string. Each string is the name of a color.
*/
vector<string> ColoredCDBG::extractColorNames(const UnitigMap<HashID>& um) const {

    vector<string> v_out;

    if (!invalid){

        if (!um.isEmpty && (color_sets != nullptr)){

            const ColorSet* cs = getColorSet(um);

            if (cs != nullptr){

                size_t prev_color_id = 0xffffffffffffffff;

                const size_t end = um.dist + um.len;
                const size_t um_km_sz = um.size - getK() + 1;

                for (ColorSet::const_iterator it = cs->begin(), it_end = cs->end(); it != it_end; ++it){

                    const size_t color_id = *it / um_km_sz;
                    const size_t km_dist = *it - (color_id * um_km_sz);

                    if ((km_dist >= um.dist) && (km_dist < end) && (color_id != prev_color_id)){

                        v_out.push_back(color_names[color_id]);
                    }

                    prev_color_id = color_id;
                }
            }
        }
    }
    else cerr << "ColoredCDBG::extractColors(): Graph is invalid, no colors can be extracted." << endl;

    return v_out;
}

/** Get the color set of a unitig.
* @param um is a UnitigMap representing a mapping to a unitig from the graph.
* @return a constant pointer to the color set matching um. If no such color set
* is found, the pointer is nullptr.
*/
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

/** Write a colored and compacted de Bruijn graph to disk.
* @param prefix_output_filename is a string which is the prefix of the filename for the two files that are
* going to be written to disk. If this prefix is "XXX", two files "XXX.gfa" and "XXX.bfg_colors" will be
* written to disk.
* @param nb_threads is the number of threads that can be used to write the graph to disk.
* @param verbose is a boolean indicating if information message are printed during writing (true) or not (false).
* @return a boolean indicating if the graph was successfully written.
*/
bool ColoredCDBG::write(const string prefix_output_filename, const size_t nb_threads, const bool verbose){

    if (CompactedDBG::write(prefix_output_filename, nb_threads, true, verbose)){

        if (verbose) cout << endl << "ColoredCDBG::write(): Writing colors to disk" << endl;

        const string out = prefix_output_filename + ".bfg_colors";

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
        const size_t nb_colors = color_names.size();

        //Write the file format version number
        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&format_version), sizeof(size_t));
        //Write k-mer length
        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&k_), sizeof(size_t));
        //Write number of different seeds for hash function
        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_seeds), sizeof(size_t));
       //Write number of colors in the graph
        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_colors), sizeof(size_t));
        //Write number of color sets in the graph
        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_color_sets), sizeof(size_t));
        //Write number of (kmer, color set) overflowing
        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&km_overflow_sz), sizeof(size_t));

        for (size_t i = 0; (i < nb_seeds) && colors_out.good(); ++i){
            //Write the hash function seeds of the graph
            colors_out.write(reinterpret_cast<const char*>(&seeds[i]), sizeof(uint64_t));
        }

        for (size_t i = 0; (i < nb_colors) && colors_out.good(); ++i){
            //Write the color names of the graph
            colors_out.write(color_names[i].c_str(), color_names[i].length() + 1);
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

void ColoredCDBG::initColorSets(const CCDBG_Build_opt& opt, const size_t max_nb_hash){

    const size_t nb_locks = opt.nb_threads * 256;

    mutex mutex_km_overflow;

    std::atomic_flag* cs_locks = new std::atomic_flag[nb_locks];

    for (size_t i = 0; i < nb_locks; ++i) cs_locks[i].clear();

    nb_seeds = max_nb_hash;
    nb_color_sets = size();

    color_sets = new ColorSet[nb_color_sets];

    color_names = opt.filename_seq_in;

    auto worker_function = [&](ColoredCDBG::iterator a, ColoredCDBG::iterator b){

        int i;

        uint64_t h_v, id_lock;

        for (auto& unitig = a; unitig != b; ++unitig){

            const Kmer head = unitig->getHead();

            for (i = 0; i < max_nb_hash; ++i){

                h_v = head.hash(seeds[i]) % nb_color_sets; // Hash to which we can possibly put our colorset for current kmer
                id_lock = h_v % nb_locks; // Lock ID for this hash

                while (cs_locks[id_lock].test_and_set(std::memory_order_acquire)); // Set the corresponding lock

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

void ColoredCDBG::buildColorSets(const size_t nb_threads){

    const size_t nb_locks = nb_threads * 256;

    const int k_ = getK();

    const size_t chunk_size = 1000;

    size_t prev_file_id = 0;

    bool next_file = true;

    FileParser fp(color_names);

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

                        um.len = 1 + um.lcp(read_color.first.c_str(), it_km->second + k_, um.strand ? um.dist + k_ : um.dist - 1, !um.strand);

                        if (!um.isShort && !um.isAbundant && !um.strand) um.dist -= um.len - 1;

                        it_km += um.len - 1;
                    }

                    HashID* hid = um.getData();

                    if (hid->get() == 0){

                        unique_lock<mutex> lock(mutex_km_overflow);

                        setColor(um, read_color.second);
                    }
                    else {

                        const uint64_t id_lock = getHash(um) % nb_locks;

                        while (cs_locks[id_lock].test_and_set(std::memory_order_acquire)); // Set the corresponding lock

                        setColor(um, read_color.second);

                        cs_locks[id_lock].clear(std::memory_order_release);
                    }
                }
            }
        }
    };

    /*auto reading_function = [&](vector<pair<string, size_t>>& v_read_color) {

        string s;

        size_t reads_now = 0;
        size_t file_id = prev_file_id;

        const size_t chunk = chunk_size * 100;

        while (reads_now < chunk) {

            if (fp.read(s, file_id)) {

                v_read_color.emplace_back(make_pair(s, file_id));

                reads_now += s.length();
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
    };*/

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
        vector<vector<pair<string, size_t>>> reads_colors(nb_threads);

        mutex mutex_file;

        while (next_file){

            stop = false;

            //cout << "prev_file_id = " << prev_file_id << endl;

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

            for (size_t i = 0; i < nb_color_sets; ++i) color_sets[i].optimize();

            for (KmerHashTable<ColorSet>::iterator it = km_overflow.begin(), it_end = km_overflow.end(); it != it_end; ++it) it->optimize();
        }
    }

    for (size_t i = 0; i < nb_color_sets; ++i) color_sets[i].optimize();

    for (KmerHashTable<ColorSet>::iterator it = km_overflow.begin(), it_end = km_overflow.end(); it != it_end; ++it) it->optimize();

    fp.close();

    delete[] cs_locks;
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

    Kmer km;

    size_t format_version , k_, km_overflow_sz, nb_colors;

    ifstream colorsfile_in;
    istream colors_in(nullptr);

    char* buffer = nullptr;

    colorsfile_in.open(opt.filename_colors_in[0].c_str(), ios_base::in | ios_base::binary);
    colors_in.rdbuf(colorsfile_in.rdbuf());

    //Read the file format version number
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&format_version), sizeof(size_t));
    //Read k-mer length
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&k_), sizeof(size_t));
    //Read number of different seeds for hash function
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_seeds), sizeof(size_t));
    //Read number of colors in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_colors), sizeof(size_t));
    //Read number of color sets in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_color_sets), sizeof(size_t));
    //Read number of (kmer, color set) overflowing
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

    color_names.clear();

    km_overflow = KmerHashTable<ColorSet>(km_overflow_sz);

    color_sets = new ColorSet[nb_color_sets];

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

    for (size_t i = 0; (i < nb_color_sets) && colors_in.good(); ++i) color_sets[i].read(colors_in);

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
