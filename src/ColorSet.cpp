#include "ColorSet.hpp"

/** Color set constructor (set up an empty color set). The color set is
* initialized as "unoccupied": the color set is "free" to be used,
* it is not associated with any unitig.
*/
ColorSet::ColorSet() : setBits(unoccupied) {}

/** Color set copy constructor. After the call to this constructor,
* the same color set exists twice in memory.
* @param o is the color set to copy.
*/
ColorSet::ColorSet(const ColorSet& o){

    const uintptr_t flag = o.setBits & flagMask;

    if (flag == ptrCompressedBitmap){

        setPointer = new Bitmap(*(o.getConstPtrBitmap()));

        setBits &= pointerMask;
    }
    else setBits = o.setBits;
}

/** Color set move constructor. After the call to this constructor,
* the color set to move is empty, its content has been transfered (moved)
* to the new color set.
* @param o is the color set to move.
*/
ColorSet::ColorSet(ColorSet&& o){

    setBits = o.setBits;
    o.setBits = unoccupied;
}

/** Color set destructor.
*/
ColorSet::~ColorSet() {

    releasePointer();
}

/** Color set copy assignment operator. After the call to this operator,
* the same color set exists twice in memory.
* @param o is the color set to copy.
* @return a reference to the current color set being a copy of o.
*/
ColorSet& ColorSet::operator=(const ColorSet& o){

    const uintptr_t flag = o.setBits & flagMask;

    releasePointer();

    if (flag == ptrCompressedBitmap){

        setPointer = new Bitmap(*(o.getConstPtrBitmap()));

        setBits &= pointerMask;
    }
    else setBits = o.setBits;

    return *this;
}

/** Color set move assignment operator. After the call to this operator,
* the color set to move is empty, its content has been transfered (moved)
* to another color set.
* @param o is the color set to move.
* @return a reference to the current color set having (owning) the content of o.
*/
ColorSet& ColorSet::operator=(ColorSet&& o){

    if (this != &o) {

        releasePointer();

        setBits = o.setBits;
        o.setBits = unoccupied;
    }

    return *this;
}

/** Empty a color set
*/
void ColorSet::empty(){

    releasePointer();
    setBits = localBitVectorColor;
}

/** Empty a color set and set it as "unoccupied": the color set is "free" to be used,
* it is not associated with any unitig.
*/
void ColorSet::setUnoccupied(){

    releasePointer();
    setBits = unoccupied;
}

/** Add a color in the current color set for a unitig or a sub-unitig.
* @param um is a UnitigMap representing the mapping of a unitig for which the color must be added.
* The color will be added only for the sub-unitig mapped, i.e, unitig[um.dist..um.dist+um.len+k-1]
* @param color_id is the ID of the color to add.
*/
void ColorSet::add(const UnitigMap<HashID>& um, const size_t color_id) {

    const size_t um_km_sz = um.size - Kmer::k + 1;
    size_t color_id_start = um_km_sz * color_id + um.dist;
    const size_t color_id_end = color_id_start + std::min(um_km_sz - um.dist, um.len);

    if ((setBits & flagMask) == unoccupied) setBits = localBitVectorColor;

    const uintptr_t flag = setBits & flagMask;

    if (flag == localBitVectorColor){

        if ((setBits == localBitVectorColor) && (um.len == 1)) setBits = (color_id_start << 2) | localSingleColor;
        else if ((color_id_end - 1) < maxBitVectorIDs){

            for (; color_id_start < color_id_end; ++color_id_start) setBits |= 1ULL << (color_id_start + 2);
        }
        else {

            uintptr_t setBits_tmp = setBits >> 2;

            setPointer = new Bitmap;

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

            for (; color_id_start < color_id_end; ++color_id_start) setBits |= 1ULL << (color_id_start + 2);
        }
        else {

            setPointer = new Bitmap;

            setPointer->add(setBits_tmp);

            for (; color_id_start < color_id_end; ++color_id_start) setPointer->add(color_id_start);

            setBits &= pointerMask;
        }
    }
    else { // flag == ptrCompressedBitmap

        Bitmap* bitmap = getPtrBitmap();

        for (; color_id_start < color_id_end; ++color_id_start) bitmap->add(color_id_start);
    }
}

/** Check if a color is present on a unitig or a sub-unitig.
* @param um is a UnitigMap representing the mapping of a unitig for which the color presence must
* be checked. The color will be checked only for the sub-unitig mapped, i.e,
* unitig[um.dist..um.dist+um.len+k-1]
* @param color_id is the ID of the color to check.
*/
bool ColorSet::contains(const UnitigMap<HashID>& um, const size_t color_id) const {

    const uintptr_t flag = setBits & flagMask;
    const size_t um_km_sz = um.size - Kmer::k + 1;
    size_t color_id_start = um_km_sz * color_id + um.dist;
    const size_t color_id_end = color_id_start + std::min(um_km_sz - um.dist, um.len);

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

    return false; //flag == unoccupied
}

/** Reverse the order of the colors. This function must be used if a unitig is reverse-complemented.
* @param len_unitig is the length of the unitig to which this color set is associated.
* @return a new color set which is the reverse of the current one.
*/
ColorSet ColorSet::reverse(const size_t len_unitig) const {

    const size_t len_unitig_km = len_unitig - Kmer::k + 1;

    ColorSet new_cs;
    ColorSet::const_iterator it = begin(), it_end = end();

    for (; it != it_end; ++it){

        const size_t color_id = *it / len_unitig_km;
        const size_t km_dist = *it - (color_id * len_unitig_km);
        const size_t new_km_dist = len_unitig_km - km_dist - 1;

        new_cs.add(len_unitig_km * color_id + new_km_dist);
    }

    return new_cs;
}

/** Get the number of colors (total number of colors associated with a unitig).
* @return number of colors (total number of colors associated with a unitig) in
* this color set.
*/
size_t ColorSet::size() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrCompressedBitmap) return getConstPtrBitmap()->cardinality();
    else if (flag == localBitVectorColor) return __builtin_popcount(setBits & pointerMask);
    else if (flag == localSingleColor) return 1;

    return 0; //flag == unoccupied
}

/** Write a color set.
* @param stream_out is an out stream to which the color set will be written. It must be
* opened prior to the call of this function and it won't be closed by this function.
* @return a boolean indicating if the write was successful.
*/
bool ColorSet::write(ostream& stream_out) const {

    if (stream_out.good()){

        const uintptr_t flag = setBits & flagMask;

        if (flag == ptrCompressedBitmap){

            const uint32_t expected_sz = getConstPtrBitmap()->getSizeInBytes();

            const uintptr_t flag_expected_sz = (static_cast<uintptr_t>(expected_sz) << 32) | flag;

            char* serialized = new char[expected_sz];

            getConstPtrBitmap()->write(serialized);

            stream_out.write(reinterpret_cast<const char*>(&flag_expected_sz), sizeof(uintptr_t));
            stream_out.write(serialized, expected_sz);

            delete[] serialized;
        }
        else stream_out.write(reinterpret_cast<const char*>(&setBits), sizeof(uintptr_t));

        return true;
    }

    return false;
}

bool ColorSet::read(istream& stream_in) {

    if (stream_in.good()){

        stream_in.read(reinterpret_cast<char*>(&setBits), sizeof(uintptr_t));

        const uintptr_t flag = setBits & flagMask;

        if (flag == ptrCompressedBitmap){

            const uint32_t expected_sz = static_cast<uint32_t>(setBits >> 32);

            char* serialized = new char[expected_sz];
            setPointer = new Bitmap;

            stream_in.read(serialized, expected_sz);

            *setPointer = std::move(Bitmap::read(serialized));

            setBits &= pointerMask;

            delete[] serialized;
        }

        return true;
    }

    return false;
}

// PRIVATE
void ColorSet::add(const size_t color_id) {

    if ((setBits & flagMask) == unoccupied) setBits = localBitVectorColor;

    const uintptr_t flag = setBits & flagMask;

    if (flag == localBitVectorColor){

        if (setBits == localBitVectorColor) setBits = (color_id << 2) | localSingleColor;
        else if (color_id < maxBitVectorIDs) setBits |= 1ULL << (color_id + 2);
        else {

            uintptr_t setBits_tmp = setBits >> 2;

            setPointer = new Bitmap;

            for (size_t i = 0; setBits_tmp != 0; ++i, setBits_tmp >>= 1) {

                if ((setBits_tmp & 0x1) != 0) setPointer->add(i);
            }

            setPointer->add(color_id);

            setBits &= pointerMask;
        }
    }
    else if (flag == localSingleColor){

        const uintptr_t setBits_tmp = setBits >> 2;

        if ((setBits_tmp < maxBitVectorIDs) && (color_id < maxBitVectorIDs)){

            setBits = (1ULL << (setBits_tmp + 2)) | (1ULL << (color_id + 2)) | localBitVectorColor;
        }
        else {

            setPointer = new Bitmap;

            setPointer->add(setBits_tmp);
            setPointer->add(color_id);

            setBits &= pointerMask;
        }
    }
    else getPtrBitmap()->add(color_id); // flag == ptrCompressedBitmap
}
