#ifndef BFG_COLOR_TCC
#define BFG_COLOR_TCC

/** Color set constructor (set up an empty color set). The color set is
* initialized as "unoccupied": the color set is "free" to be used,
* it is not associated with any unitig.
*/
template<typename U>
UnitigColors<U>::UnitigColors() : setBits(unoccupied) {}

/** Color set copy constructor. After the call to this constructor,
* the same color set exists twice in memory.
* @param o is the color set to copy.
*/
template<typename U>
UnitigColors<U>::UnitigColors(const UnitigColors& o){

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
template<typename U>
UnitigColors<U>::UnitigColors(UnitigColors&& o){

    setBits = o.setBits;
    o.setBits = unoccupied;
}

/** Color set destructor.
*/
template<typename U>
UnitigColors<U>::~UnitigColors() {

    releasePointer();
}

/** Color set copy assignment operator. After the call to this operator,
* the same color set exists twice in memory.
* @param o is the color set to copy.
* @return a reference to the current color set being a copy of o.
*/
template<typename U>
UnitigColors<U>& UnitigColors<U>::operator=(const UnitigColors& o){

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
template<typename U>
UnitigColors<U>& UnitigColors<U>::operator=(UnitigColors&& o){

    if (this != &o) {

        releasePointer();

        setBits = o.setBits;
        o.setBits = unoccupied;
    }

    return *this;
}

/** Empty a color set
*/
template<typename U>
void UnitigColors<U>::empty(){

    releasePointer();
    setBits = localBitVectorColor;
}

/** Add a color in the current color set for a unitig or a sub-unitig.
* @param um is a UnitigMap representing the mapping of a unitig for which the color must be added.
* The color will be added only for the sub-unitig mapped, i.e, unitig[um.dist..um.dist+um.len+k-1]
* @param color_id is the ID of the color to add.
*/
template<typename U>
void UnitigColors<U>::add(const const_UnitigColorMap<U>& um, const size_t color_id) {

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
template<typename U>
bool UnitigColors<U>::contains(const const_UnitigColorMap<U>& um, const size_t color_id) const {

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
template<typename U>
UnitigColors<U> UnitigColors<U>::reverse(const const_UnitigColorMap<U>& um) const {

    const size_t len_unitig_km = um.size - Kmer::k + 1;

    UnitigColors<U> new_cs;
    UnitigColors<U>::const_iterator it = begin(), it_end = end();

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
template<typename U>
size_t UnitigColors<U>::size() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrCompressedBitmap) return getConstPtrBitmap()->cardinality();
    else if (flag == localBitVectorColor) return __builtin_popcountll(setBits & pointerMask);
    else if (flag == localSingleColor) return 1;

    return 0; //flag == unoccupied
}

/** Write a color set.
* @param stream_out is an out stream to which the color set will be written. It must be
* opened prior to the call of this function and it won't be closed by this function.
* @return a boolean indicating if the write was successful.
*/
template<typename U>
bool UnitigColors<U>::write(ostream& stream_out) const {

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

template<typename U>
bool UnitigColors<U>::read(istream& stream_in) {

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
template<typename U>
void UnitigColors<U>::add(const size_t color_id) {

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

template<typename U>
size_t UnitigColors<U>::getSizeInBytes() const {

    if ((setBits & flagMask) == ptrCompressedBitmap) return sizeof(UnitigColors<U>) + getPtrBitmap()->getSizeInBytes();
    return sizeof(UnitigColors);
}

/*template<typename U>
void UnitigColors<U>::test(const const_UnitigColorMap<U>& um) {

    UnitigColors<U> new_cs;
    UnitigColors<U>::const_iterator it = begin(), it_end = end(), prev_it = begin(), it_tmp;

    const size_t len_unitig_km = um.size - Kmer::k + 1;
    const size_t newlen_unitig_km = len_unitig_km + 1;

    size_t prev_color_id = 0;
    size_t count_km_per_color = 0;

    while (it != it_end){ // Iterate over colors

        const size_t color_id = *it / len_unitig_km; // Current color

        it_tmp = it;
        ++it_tmp;

        // If the current color is not the same as the previous one
        // or we have reached  the end of the UnitigColors
        if ((it_tmp == it_end) || (color_id != prev_color_id)){

            // For the current color, if only half of the k-mers in the unitigs have it
            if (count_km_per_color >= (len_unitig_km/2)){

                new_cs.add(newlen_unitig_km * prev_color_id);

                size_t prev_km_dist = 0xffffffffffffffff;

                while (prev_it != it){

                    const size_t km_dist = *prev_it - (prev_color_id * len_unitig_km);

                    it_tmp = prev_it;
                    ++it_tmp;

                    if ((it_tmp == it) || (km_dist != (prev_km_dist + 1))){

                        for (++prev_km_dist; prev_km_dist != km_dist; ++prev_km_dist){

                            new_cs.add(newlen_unitig_km * prev_color_id + (prev_km_dist + 1));
                        }
                    }
                    else ++prev_km_dist;

                    ++prev_it;
                }
            }
            else {

                for (; prev_it != it; ++prev_it){

                    const size_t new_km_dist = *prev_it - (prev_color_id * len_unitig_km);

                    new_cs.add(newlen_unitig_km * prev_color_id + (new_km_dist + 1));
                }
            }

            prev_color_id = color_id;
            count_km_per_color = 1;
        }
        else ++count_km_per_color;

        ++it;
    }

    *this = new_cs;
}*/

#endif
