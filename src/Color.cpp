#include "Color.hpp"
#include "ColoredCDBG.hpp"

Color::Color(const uint8_t c) : color_id(c) {}

void Color::join(const UnitigMap<Color>& um_dest, const UnitigMap<Color>& um_src){

    ColoredCDBG* colored_cdbg_dest = static_cast<ColoredCDBG*>(um_dest.cdbg);

    colored_cdbg_dest->setColors(um_dest, um_src);
}

void Color::split(const UnitigMap<Color>& um_split, const size_t pos, const size_t len, Color& data_dest) const {

    data_dest.color_id = color_id;
}

ColorSet::ColorSet() : asBits(0) {}

ColorSet::ColorSet(const size_t color) : asBits(0) {

    if (color < maxBitVectorIDs) asBits |= 1ULL << (color + 2);
    else asBits = (color << 2) | localSingleColor;
}

ColorSet::~ColorSet() {

    //releasePointer();
}

void ColorSet::insertColor(const size_t color) { //Thread safe

    const uintptr_t flag = asBits & flagMask;

    if ((flag == localBitVectorColor) || (flag == localSingleColor)){

            //__sync_fetch_and_or(&color_sets[head.hash(seeds[color_id - 1]) % nb_color_sets], 1ULL << file_id);
    }
    else { //Pointer mode or multi k-mer mode, not supported yet
    }
}
