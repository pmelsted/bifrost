#ifndef BFG_COLOREDCDBG_HPP
#define BFG_COLOREDCDBG_HPP

#include <iostream>
#include <random>

#include "Color.hpp"
#include "CompactedDBG.hpp"

class ColoredCDBG : public CompactedDBG<Color> {

    public:

        ColoredCDBG(int kmer_length = DEFAULT_K, int minimizer_length = DEFAULT_G);
        ~ColoredCDBG();

        bool build(const CDBG_Build_opt& opt);

        ColorSet* getColorSet(const UnitigMap<Color>& um);

        void setColor(const UnitigMap<Color>& um, size_t color);
        void setColors(const UnitigMap<Color>& um_dest, const UnitigMap<Color>& um_src);

    private:

        void initColorSets(const size_t max_nb_hash = 255);
        void mapColors(const CDBG_Build_opt& opt);

        uint64_t seeds[256];

        size_t nb_color_sets;

        ColorSet* color_sets;

        KmerHashTable<ColorSet> km_overflow;
};

#endif
