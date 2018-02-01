#include "Color.hpp"
#include "ColoredCDBG.hpp"

HashID::HashID(const uint8_t hid) : hash_id(hid) {}

void HashID::join(const UnitigMap<HashID>& um_dest, const UnitigMap<HashID>& um_src){

    ColoredCDBG* colored_cdbg = static_cast<ColoredCDBG*>(um_dest.cdbg);

    ColorSet* cs_dest = colored_cdbg->getColorSet(um_dest);
    ColorSet* cs_src = colored_cdbg->getColorSet(um_src);

    const HashID hid(0);

    if (cs_dest != nullptr){ // If a colorset exists for um_dest

        const Kmer new_head = um_dest.strand ? um_dest.getHead() : um_dest.getTail().rep();

        colored_cdbg->joinColors(um_dest, um_src); // Join the color sets

        // TODO: Insert in tombstone if available
        // TODO: if new_head = head, do not insert in overflow

        // Insert new colorset with corresponding head into overflow of k-mers
        colored_cdbg->km_overflow.insert(new_head, ColorSet(*cs_dest));

        cs_dest->setUnoccupied(); // Set um_dest colorset as a tombstone

        um_dest.setData(&hid);
    }

    if (cs_src != nullptr){

        cs_src->setUnoccupied();
        // TODO: if new_head = head, do not insert in overflow
        um_src.setData(&hid);
    }
}

void HashID::sub(const UnitigMap<HashID>& um, HashID& data_dest, const bool last_extraction) const {

    ColoredCDBG* colored_cdbg = static_cast<ColoredCDBG*>(um.cdbg);

    ColorSet cs = colored_cdbg->extractColors(um);

    if (cs.size() != 0){

        const Kmer km = um.getKmer(um.dist);

        // TODO: Insert in tombstone if available
        colored_cdbg->km_overflow.insert(km, cs);

        if (last_extraction){

            const HashID hid(0);

            colored_cdbg->getColorSet(um)->setUnoccupied();

            um.setData(&hid);
        }
    }
}

string HashID::serialize() const { return std::to_string(hash_id); }






ColorSet::ColorSet() : setBits(unoccupied) {}

ColorSet::~ColorSet() {

    releasePointer();
}

ColorSet::ColorSet(const ColorSet& o){

    const uintptr_t flag = o.setBits & flagMask;

    if (flag == ptrCompressedBitmap){

        setPointer = new Bitmap(*(o.getConstPtrBitmap()));

        setBits &= pointerMask;
    }
    else setBits = o.setBits;
}

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

void ColorSet::empty(){

    releasePointer();
    setBits = localBitVectorColor;
}

void ColorSet::setUnoccupied(){

    releasePointer();
    setBits = unoccupied;
}

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

size_t ColorSet::size() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrCompressedBitmap) return getConstPtrBitmap()->cardinality();
    else if (flag == localBitVectorColor) return __builtin_popcount(setBits & pointerMask);
    else if (flag == localSingleColor) return 1;

    return 0; //flag == unoccupied
}

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
