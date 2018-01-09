#ifndef BFG_COLOREDCDBG_HPP
#define BFG_COLOREDCDBG_HPP

#include <iostream>
#include <random>

#include "Color.hpp"
#include "CompactedDBG.hpp"

class ColoredCDBG : public CompactedDBG<HashID> {

    friend class HashID;

    public:

        ColoredCDBG(int kmer_length = DEFAULT_K, int minimizer_length = DEFAULT_G);
        ~ColoredCDBG();

        bool build(const CDBG_Build_opt& opt);

        bool setColor(const UnitigMap<HashID>& um, size_t color);
        bool joinColors(const UnitigMap<HashID>& um_dest, const UnitigMap<HashID>& um_src);
        ColorSet extractColors(const UnitigMap<HashID>& um, const size_t pos, const size_t len) const;

    private:

        void initColorSets(const size_t max_nb_hash = 127);
        void mapColors(const CDBG_Build_opt& opt);

        ColorSet* getColorSet(const UnitigMap<HashID>& um);
        const ColorSet* getColorSet(const UnitigMap<HashID>& um) const;

        uint64_t seeds[256];

        size_t nb_color_sets;

        ColorSet* color_sets;

        KmerHashTable<ColorSet> km_overflow;
};

#endif
