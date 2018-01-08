#include "Color.hpp"
#include "ColoredCDBG.hpp"

HashID::HashID(const uint8_t hid) : hash_id(hid) {}

void HashID::lock() {

    uint8_t hash_id_unlocked = hash_id & 0x7f;
    uint8_t hash_id_locked = hash_id | 0x80;

    while (!__sync_bool_compare_and_swap(&hash_id, hash_id_unlocked, hash_id_locked)) {

        hash_id_unlocked = hash_id & 0x7f;
        hash_id_locked = hash_id | 0x80;
    }
}

void HashID::join(const UnitigMap<HashID>& um_dest, const UnitigMap<HashID>& um_src){

    ColoredCDBG* colored_cdbg_dest = static_cast<ColoredCDBG*>(um_dest.cdbg);

    colored_cdbg_dest->joinColors(um_dest, um_src);
}

void HashID::split(const UnitigMap<HashID>& um_split, const size_t pos, const size_t len, HashID& data_dest) const {

    //TODO
}








ColorSet::ColorSet() : setBits(localBitVectorColor) {}

ColorSet::~ColorSet() {

    releasePointer();
}

void ColorSet::releasePointer(){

    if ((setBits & flagMask) == ptrCompressedBitmap){

        delete getPtrBitmap();
        setBits = localBitVectorColor;
    }
}

void ColorSet::add(const UnitigMap<HashID>& um, const size_t color_id) {

    const uintptr_t flag = setBits & flagMask;
    size_t color_id_start = um.size * color_id + color_id;
    const size_t color_id_end = color_id_start + (std::min(um.size, um.dist + um.len) - um.dist);

    if (flag == ptrCompressedBitmap){

        Bitmap* bitmap = getPtrBitmap();

        for (; color_id_start < color_id_end; ++color_id_start) bitmap->add(color_id_start);
    }
    else if (flag == localBitVectorColor){

        if ((setBits == localBitVectorColor) && (um.len == 1)) setBits = (color_id_start << 2) | localSingleColor;
        else if ((color_id_end - 1) < maxBitVectorIDs){

            for (color_id_start += 2; color_id_start < color_id_end; ++color_id_start) setBits |= 1ULL << color_id_start;
        }
        else {

            uintptr_t setBits_tmp = setBits >> 2;

            setPointer = new Roaring;

            for (size_t i = 0; setBits_tmp != 0; ++i, setBits_tmp >>= 1) {

                if ((setBits_tmp & 0x1) != 0) setPointer->add(i);
            }

            for (; color_id_start < color_id_end; ++color_id_start) setPointer->add(color_id_start);

            setBits &= pointerMask;
        }
    }
    else if (flag == localSingleColor){

        const uintptr_t setBits_tmp = setBits >> 2;

        if ((setBits_tmp < maxBitVectorIDs) && ((color_id_end - 1) < maxBitVectorIDs)){

            setBits = (1ULL << (setBits_tmp + 2)) | localBitVectorColor;

            for (color_id_start += 2; color_id_start < color_id_end; ++color_id_start) setBits |= 1ULL << color_id_start;
        }
        else {

            setPointer = new Roaring;

            setPointer->add(setBits_tmp);

            for (; color_id_start < color_id_end; ++color_id_start) setPointer->add(color_id_start);

            setBits &= pointerMask;
        }
    }
}

bool ColorSet::contains(const UnitigMap<HashID>& um, const size_t color_id) const {

    const uintptr_t flag = setBits & flagMask;
    size_t color_id_start = um.size * color_id + color_id;
    const size_t color_id_end = color_id_start + (std::min(um.size, um.dist + um.len) - um.dist);

    if (flag == ptrCompressedBitmap){

        const Bitmap* bitmap = getConstPtrBitmap();

        for (; color_id_start < color_id_end; ++color_id_start){

            if (!bitmap->contains(color_id_start)) return false;
        }

        return true;
    }
    else if (flag == localBitVectorColor){

        if ((color_id_end - 1) < maxBitVectorIDs){

            uintptr_t setBits_tmp = setBits >> (color_id_start + 2);

            for (; color_id_start < color_id_end; ++color_id_start, setBits_tmp >>= 1){

                if ((setBits_tmp & 0x1) == 0) return false;
            }

            return true;
        }
    }
    else if ((flag == localSingleColor) && (um.len == 1)) return color_id_start == (setBits >> 2);

    return false;
}

size_t ColorSet::size() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrCompressedBitmap) return getConstPtrBitmap()->cardinality();
    else if (flag == localBitVectorColor) return __builtin_popcount(setBits & pointerMask);
    else if (flag == localSingleColor) return 1;

    return 0;
}
