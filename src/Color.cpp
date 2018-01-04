#include "Color.hpp"
#include "ColoredCDBG.hpp"

Color::Color(const uint8_t c) : color_id(c) {}

void Color::getLock() {

    uint8_t color_id_unlocked = color_id & 0x7f;
    uint8_t color_id_locked = color_id | 0x80;

    while (!__sync_bool_compare_and_swap(&color_id, color_id_unlocked, color_id_locked)) {

        color_id_unlocked = color_id & 0x7f;
        color_id_locked = color_id | 0x80;
    }
}

void Color::join(const UnitigMap<Color>& um_dest, const UnitigMap<Color>& um_src){

    ColoredCDBG* colored_cdbg_dest = static_cast<ColoredCDBG*>(um_dest.cdbg);

    colored_cdbg_dest->setColors(um_dest, um_src);
}

void Color::split(const UnitigMap<Color>& um_split, const size_t pos, const size_t len, Color& data_dest) const {

    data_dest.color_id = color_id;
}

ColorSet::ColorSet() : setBits(localBitVectorColor) {}

ColorSet::ColorSet(const size_t color) : setBits(localBitVectorColor) { add(color); }

ColorSet::~ColorSet() {

    releasePointer();
}

void ColorSet::releasePointer(){

    if ((setBits & flagMask) == pointerCompressedBitmap){

        delete getPointer();
        setBits = localBitVectorColor;
    }
}

void ColorSet::add(const size_t color) {

    const uintptr_t flag = setBits & flagMask;

    if (flag == pointerCompressedBitmap) getPointer()->add(color);
    else if (flag == localBitVectorColor){

        if (color < maxBitVectorIDs) setBits |= (1ULL << (color + 2)) | localBitVectorColor;
        else if (setBits == localBitVectorColor) setBits = (color << 2) | localSingleColor;
        else {

            uintptr_t setBits_tmp = setBits >> 2;

            setPointer = new Roaring;

            for (size_t i = 0; setBits_tmp != 0; ++i, setBits_tmp >>= 1) {

                if ((setBits_tmp & 0x1) != 0) setPointer->add(i);
            }

            setPointer->add(color);

            setBits &= pointerMask;
        }
    }
    else if (flag == localSingleColor){

        const uintptr_t setBits_tmp = setBits >> 2;

        if ((color < maxBitVectorIDs) && (setBits_tmp < maxBitVectorIDs)){

            setBits = (1ULL << (color + 2)) | localBitVectorColor;
            setBits |= 1ULL << (setBits_tmp + 2);
        }
        else {

            setPointer = new Roaring;

            setPointer->add(setBits_tmp);
            setPointer->add(color);

            setBits &= pointerMask;
        }
    }
}

bool ColorSet::contains(const size_t color) const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == pointerCompressedBitmap) return getConstPointer()->contains(color);
    else if (flag == localBitVectorColor){

        if (color < maxBitVectorIDs) return ((setBits >> (color + 2)) & 0x1) != 0;
    }
    else if (flag == localSingleColor) return color == (setBits >> 2);

    return false;
}

// Union
void ColorSet::join(const ColorSet& cs) {

    const uintptr_t flag_this = setBits & flagMask;
    const uintptr_t flag_cs = cs.setBits & flagMask;

    if ((flag_this == pointerCompressedBitmap) && (flag_cs == pointerCompressedBitmap)) *getPointer() |= *(cs.getConstPointer());
    else if ((flag_this == localBitVectorColor) && (flag_cs == localBitVectorColor)) setBits |= cs.setBits;
    else {

        for (ColorSet::const_iterator it = cs.begin(); it != cs.end(); ++it) add(*it);
    }
}

size_t ColorSet::size() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == pointerCompressedBitmap) return getConstPointer()->cardinality();
    else if (flag == localBitVectorColor) return __builtin_popcount(setBits & pointerMask);
    else if (flag == localSingleColor) return 1;

    return 0;
}
