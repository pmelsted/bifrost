#ifndef BFG_COLOREDCDBG_HPP
#define BFG_COLOREDCDBG_HPP

#include <iostream>
#include <random>

#include "Color.hpp"
#include "CompactedDBG.hpp"

#define BFG_COLOREDCDBG_FORMAT_VERSION 1

struct CCDBG_Build_opt {

    bool reference_mode;
    bool verbose;

    size_t nb_threads;
    size_t read_chunksize;
    size_t unitig_size;
    size_t nb_unique_kmers;
    size_t nb_non_unique_kmers;
    size_t nb_bits_unique_kmers_bf;
    size_t nb_bits_non_unique_kmers_bf;

    string inFilenameBBF;
    string outFilenameBBF;

    vector<string> filename_seq_in;
    vector<string> filename_colors_in;

    // The following members are not used by CompactedDBG<T>::build
    // but you can set them to use them as parameters for other functions
    // such as CompactedDBG<T>::simplify or CompactedDBG<T>::write.

    size_t k, g;

    bool clipTips;
    bool deleteIsolated;
    bool useMercyKmers;

    bool outputGFA;
    bool outputColors;

    string prefixFilenameOut;

    CCDBG_Build_opt() : nb_threads(1), k(DEFAULT_K), g(DEFAULT_G), nb_unique_kmers(0), nb_non_unique_kmers(0),
                        nb_bits_unique_kmers_bf(14), nb_bits_non_unique_kmers_bf(14), read_chunksize(10000),
                        unitig_size(1000000), verbose(false), clipTips(false), deleteIsolated(false), useMercyKmers(false),
                        outputGFA(true), outputColors(false), reference_mode(true), inFilenameBBF(""), outFilenameBBF("") {}

    CDBG_Build_opt getCDBG_Build_opt() const {

        CDBG_Build_opt cdbg_opt;

        cdbg_opt.reference_mode = reference_mode;
        cdbg_opt.filename_in = filename_seq_in;

        cdbg_opt.verbose = verbose;

        cdbg_opt.nb_threads = nb_threads;
        cdbg_opt.read_chunksize = read_chunksize;
        cdbg_opt.unitig_size = unitig_size;
        cdbg_opt.nb_unique_kmers = nb_unique_kmers;
        cdbg_opt.nb_non_unique_kmers = nb_non_unique_kmers;
        cdbg_opt.nb_bits_unique_kmers_bf = nb_bits_unique_kmers_bf;
        cdbg_opt.nb_bits_non_unique_kmers_bf = nb_bits_non_unique_kmers_bf;

        cdbg_opt.inFilenameBBF = inFilenameBBF;
        cdbg_opt.outFilenameBBF = outFilenameBBF;

        cdbg_opt.clipTips = clipTips;
        cdbg_opt.deleteIsolated = deleteIsolated;
        cdbg_opt.useMercyKmers = useMercyKmers;

        cdbg_opt.outputGFA = outputGFA;

        cdbg_opt.prefixFilenameOut = prefixFilenameOut;

        return cdbg_opt;
    }
};

class ColoredCDBG : public CompactedDBG<HashID> {

    friend class HashID;

    public:

        ColoredCDBG(int kmer_length = DEFAULT_K, int minimizer_length = DEFAULT_G);
        ~ColoredCDBG();

        void clear();
        void empty();

        bool build(const CCDBG_Build_opt& opt);

        bool mapColors(const CCDBG_Build_opt& opt);

        bool setColor(const UnitigMap<HashID>& um, size_t color);
        bool joinColors(const UnitigMap<HashID>& um_dest, const UnitigMap<HashID>& um_src);
        ColorSet extractColors(const UnitigMap<HashID>& um) const;

        bool write(const string output_filename, const size_t nb_threads, const bool verbose);

    private:

        void initColorSets(const CCDBG_Build_opt& opt, const size_t max_nb_hash = 31);
        void buildColorSets(const CCDBG_Build_opt& opt);
        bool readColorSets(const CCDBG_Build_opt& opt);

        ColorSet* getColorSet(const UnitigMap<HashID>& um);
        const ColorSet* getColorSet(const UnitigMap<HashID>& um) const;

        void checkColors(const CCDBG_Build_opt& opt);

        uint64_t getHash(const UnitigMap<HashID>& um) const;

        uint64_t seeds[256];

        size_t nb_seeds;
        size_t nb_color_sets;

        ColorSet* color_sets;

        KmerHashTable<ColorSet> km_overflow;
};

#endif
