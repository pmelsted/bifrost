#ifndef BFG_COLOREDCDBG_HPP
#define BFG_COLOREDCDBG_HPP

#include <iostream>
#include <random>

#include "CompactedDBG.hpp"

class Color : public CDBG_Data_t<Color> {

    public:

        Color(const uint8_t c = 0) : color_id(c) {}

        void join(const Color& data, CompactedDBG<Color>& cdbg);
        void split(const size_t pos, const size_t len, Color& new_data, CompactedDBG<Color>& cdbg) const;

        inline uint8_t getColor() const { return color_id; }
        inline void setColor(const uint8_t col) { color_id = col; }

    private:

        uint8_t color_id;
};

class ColoredCDBG : public CompactedDBG<Color>{

    public:

        ColoredCDBG(int kmer_length = DEFAULT_K, int minimizer_length = DEFAULT_G) : CompactedDBG(kmer_length, minimizer_length), color_sets(nullptr) {

            std::random_device rd; //Seed
            std::default_random_engine generator(rd()); //Random number generator
            std::uniform_int_distribution<long long unsigned> distribution(0,0xFFFFFFFFFFFFFFFF); //Distribution on which to apply the generator

            for (int i = 0; i < 256; ++i) seeds[i] = distribution(generator);
        }

        ~ColoredCDBG() {

            if (color_sets != nullptr){

                delete[] color_sets;
                color_sets = nullptr;
            }
        }

        bool build(const CDBG_Build_opt& opt){

            CompactedDBG::build(opt);

            createMNPHF();

            mapColors(opt);

            return true;
        }

    private:

        void createMNPHF(const size_t max_nb_hash = 255){

            int i;

            const size_t sz_cdbg = size();
            const size_t sz_collisions = (sz_cdbg + 63) / 64;

            size_t nb_unitig_not_hashed = 0;

            uint64_t* collisions = new uint64_t[sz_collisions]();

            for (auto& unitig : *this){

                const Kmer head = unitig.getHead();

                for (i = 0; i < max_nb_hash; ++i){

                    const uint64_t h_v = head.hash(seeds[i]) % sz_cdbg;

                    if ((collisions[h_v >> 6] & (1ULL << (h_v & 0x3f))) == 0) break;
                }

                const Color col(static_cast<uint8_t>(i == max_nb_hash ? 0 : i + 1));
                unitig.setData(&col);

                if (i == max_nb_hash) km_overflow.insert(head, 0);
                else {

                    const uint64_t h_v = head.hash(seeds[i]) % sz_cdbg;
                    collisions[h_v >> 6] |= 1ULL << (h_v & 0x3f);
                }
            }

            delete[] collisions;

            color_sets = new uint64_t[sz_collisions * 64]();

            cout << "Number of unitigs not hashed is " << km_overflow.size() << " on " << sz_cdbg << " unitigs." << endl;
        }

        void mapColors(const CDBG_Build_opt& opt){

            const int k_ = getK();

            size_t file_id = 0;

            const size_t sz_cdbg = size();

            string s;

            vector<string> readv;

            FastqFile FQ(opt.fastx_filename_in);

            // Main worker thread
            auto worker_function = [&](vector<string>::iterator a, vector<string>::iterator b, vector<pair<Kmer, uint64_t>>* v_out) {

                // for each input
                for (auto x = a; x != b; ++x) {

                    for (KmerIterator it_km(x->c_str()), it_km_end; it_km != it_km_end; ++it_km) {

                        const UnitigMap<Color> um = find(it_km->first);

                        if (!um.isEmpty) {

                            const Kmer head = um.getHead();
                            const uint8_t color_id = um.getData()->getColor();

                            if (color_id == 0) v_out->push_back(make_pair(head, file_id));
                            else {

                                __sync_fetch_and_or(&color_sets[head.hash(seeds[color_id - 1]) % sz_cdbg], 1ULL << file_id);

                                if (um.strand || (um.dist != 0)){

                                    it_km += um.lcp(x->c_str(), it_km->second + k_, um.strand ? um.dist + k_ : um.dist - 1, um.strand);
                                }
                            }
                        }
                    }
                }
            };

            vector<pair<Kmer, uint64_t>> km_overflow_edit[opt.nb_threads];

            bool done = false;

            while (!done) {

                size_t reads_now = 0;

                while (reads_now < opt.read_chunksize) {

                    if (FQ.read_next(s, file_id) >= 0){

                        readv.push_back(s);
                        ++reads_now;
                    }
                    else {

                        done = true;
                        break;
                    }
                }

                // run parallel code
                vector<thread> workers;

                auto rit = readv.begin();
                size_t batch_size = readv.size() / opt.nb_threads;
                size_t leftover   = readv.size() % opt.nb_threads;

                for (size_t i = 0; i < opt.nb_threads; ++i) {

                    size_t jump = batch_size + (i < leftover ? 1 : 0);
                    auto rit_end(rit);

                    advance(rit_end, jump);
                    workers.push_back(thread(worker_function, rit, rit_end, &km_overflow_edit[i]));

                    rit = rit_end;
                }

                assert(rit == readv.end());

                for (auto& t : workers) t.join();

                for (auto& v : km_overflow_edit) {

                    for (const auto& p : v) {

                        KmerHashTable<uint64_t>::iterator it = km_overflow.find(p.first);

                        if (it != km_overflow.end()) *it |= 1ULL << p.second;
                    }

                    v.clear();
                }

                readv.clear();
            }

            FQ.close();
        }

        uint64_t seeds[256];

        uint64_t* color_sets;

        KmerHashTable<uint64_t> km_overflow;
};

#include "ColoredCDBG.tpp"

/*class ColoredCDBG {

    public:

        ColoredCDBG(int kmer_length = DEFAULT_K, int minimizer_length = DEFAULT_G) : cdbg(kmer_length, minimizer_length), color_sets(nullptr) {

            std::random_device rd; //Seed
            std::default_random_engine generator(rd()); //Random number generator
            std::uniform_int_distribution<long long unsigned> distribution(0,0xFFFFFFFFFFFFFFFF); //Distribution on which to apply the generator

            for (int i = 0; i < 256; ++i) seeds[i] = distribution(generator);
        }

        ~ColoredCDBG() {

            if (color_sets != nullptr){

                delete[] color_sets;
                color_sets = nullptr;
            }
        }

        bool build(const CDBG_Build_opt& opt){

            cdbg.build(opt);

            createMNPHF();

            mapColors(opt);

            return true;
        }

    private:

        void createMNPHF(const size_t max_nb_hash = 255){

            int i;

            const size_t sz_cdbg = cdbg.size();
            const size_t sz_collisions = (sz_cdbg + 63) / 64;

            size_t nb_unitig_not_hashed = 0;

            uint64_t* collisions = new uint64_t[sz_collisions]();

            for (auto& unitig : cdbg){

                const Kmer head = unitig.getHead();

                for (i = 0; i < max_nb_hash; ++i){

                    const uint64_t h_v = head.hash(seeds[i]) % sz_cdbg;

                    if ((collisions[h_v >> 6] & (1ULL << (h_v & 0x3f))) == 0) break;
                }

                const Color col(static_cast<uint8_t>(i == max_nb_hash ? 0 : i + 1));
                unitig.setData(&col);

                if (i == max_nb_hash) km_overflow.insert(head, 0);
                else {

                    const uint64_t h_v = head.hash(seeds[i]) % sz_cdbg;
                    collisions[h_v >> 6] |= 1ULL << (h_v & 0x3f);
                }
            }

            delete[] collisions;

            color_sets = new uint64_t[sz_collisions * 64]();

            cout << "Number of unitigs not hashed is " << km_overflow.size() << " on " << sz_cdbg << " unitigs." << endl;
        }

        void mapColors(const CDBG_Build_opt& opt){

            const int k_ = cdbg.getK();

            size_t file_id = 0;

            const size_t sz_cdbg = cdbg.size();

            string s;

            vector<string> readv;

            FastqFile FQ(opt.fastx_filename_in);

            // Main worker thread
            auto worker_function = [&](vector<string>::iterator a, vector<string>::iterator b, vector<pair<Kmer, uint64_t>>* v_out) {

                // for each input
                for (auto x = a; x != b; ++x) {

                    for (KmerIterator it_km(x->c_str()), it_km_end; it_km != it_km_end; ++it_km) {

                        const UnitigMap<Color> um = cdbg.find(it_km->first);

                        if (!um.isEmpty) {

                            const Kmer head = um.getHead();
                            const uint8_t color_id = um.getData()->getColor();

                            if (color_id == 0) v_out->push_back(make_pair(head, file_id));
                            else {

                                __sync_fetch_and_or(&color_sets[head.hash(seeds[color_id - 1]) % sz_cdbg], 1ULL << file_id);

                                if (um.strand || (um.dist != 0)){

                                    it_km += um.lcp(x->c_str(), it_km->second + k_, um.strand ? um.dist + k_ : um.dist - 1, um.strand);
                                }
                            }
                        }
                    }
                }
            };

            vector<pair<Kmer, uint64_t>> km_overflow_edit[opt.nb_threads];

            bool done = false;

            while (!done) {

                size_t reads_now = 0;

                while (reads_now < opt.read_chunksize) {

                    if (FQ.read_next(s, file_id) >= 0){

                        readv.push_back(s);
                        ++reads_now;
                    }
                    else {

                        done = true;
                        break;
                    }
                }

                // run parallel code
                vector<thread> workers;

                auto rit = readv.begin();
                size_t batch_size = readv.size() / opt.nb_threads;
                size_t leftover   = readv.size() % opt.nb_threads;

                for (size_t i = 0; i < opt.nb_threads; ++i) {

                    size_t jump = batch_size + (i < leftover ? 1 : 0);
                    auto rit_end(rit);

                    advance(rit_end, jump);
                    workers.push_back(thread(worker_function, rit, rit_end, &km_overflow_edit[i]));

                    rit = rit_end;
                }

                assert(rit == readv.end());

                for (auto& t : workers) t.join();

                for (auto& v : km_overflow_edit) {

                    for (const auto& p : v) {

                        KmerHashTable<uint64_t>::iterator it = km_overflow.find(p.first);

                        if (it != km_overflow.end()) *it |= 1ULL << p.second;
                    }

                    v.clear();
                }

                readv.clear();
            }

            FQ.close();
        }

        uint64_t seeds[256];

        uint64_t* color_sets;

        KmerHashTable<uint64_t> km_overflow;

        CompactedDBG<Color> cdbg;
};*/

#endif
