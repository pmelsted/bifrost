#include "ColorSet.hpp"

UnitigColors::UnitigColors() : setBits(unoccupied) {}

UnitigColors::UnitigColors(const UnitigColors& o){

    const uintptr_t flag = o.setBits & flagMask;

    if (flag == ptrCompressedBitmap){

        setPointer = new Bitmap(*(o.getConstPtrBitmap()));

        setBits &= pointerMask;
    }
    else setBits = o.setBits;
}

UnitigColors::UnitigColors(UnitigColors&& o){

    setBits = o.setBits;
    o.setBits = unoccupied;
}

UnitigColors::~UnitigColors() {

    releasePointer();
}

UnitigColors& UnitigColors::operator=(const UnitigColors& o){

    const uintptr_t flag = o.setBits & flagMask;

    releasePointer();

    if (flag == ptrCompressedBitmap){

        setPointer = new Bitmap(*(o.getConstPtrBitmap()));

        setBits &= pointerMask;
    }
    else setBits = o.setBits;

    return *this;
}

UnitigColors& UnitigColors::operator=(UnitigColors&& o){

    if (this != &o) {

        releasePointer();

        setBits = o.setBits;
        o.setBits = unoccupied;
    }

    return *this;
}

void UnitigColors::empty(){

    releasePointer();
    setBits = localBitVectorColor;
}

size_t UnitigColors::getSizeInBytes() const {

    if ((setBits & flagMask) == ptrCompressedBitmap) return getPtrBitmap()->getSizeInBytes() + sizeof(Bitmap) + sizeof(UnitigColors);
    return sizeof(UnitigColors);
}

void UnitigColors::add(const UnitigMapBase& um, const size_t color_id) {

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

void UnitigColors::remove(const UnitigMapBase& um, const size_t color_id) {

    const uintptr_t flag = setBits & flagMask;

    if (flag == unoccupied) return;

    const size_t um_km_sz = um.size - Kmer::k + 1;
    size_t color_id_start = um_km_sz * color_id + um.dist;
    const size_t color_id_end = color_id_start + std::min(um_km_sz - um.dist, um.len);

    if (flag == localBitVectorColor){

        uintptr_t mask = 0;

        for (; color_id_start < min(maxBitVectorIDs, color_id_end); ++color_id_start) mask |= 1ULL << (color_id_start + 2);

        setBits &= ~mask;

        if (setBits == localBitVectorColor) setBits = localSingleColor;
    }
    else if (flag == localSingleColor){

        const uintptr_t setBits_tmp = setBits >> 2;

        if ((setBits_tmp >= color_id_start) && (setBits_tmp <= color_id_end)) setBits = localSingleColor;
    }
    else { // flag == ptrCompressedBitmap

        Bitmap* bitmap = getPtrBitmap();

        for (; color_id_start < color_id_end; ++color_id_start) bitmap->remove(color_id_start);

        if (bitmap->cardinality() == 0) empty();
    }
}

void UnitigColors::add(const size_t color_id) { // PRIVATE

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

bool UnitigColors::contains(const UnitigMapBase& um, const size_t color_id) const {

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

size_t UnitigColors::size() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrCompressedBitmap) return getConstPtrBitmap()->cardinality();
    else if (flag == localBitVectorColor) return __builtin_popcountll(setBits & pointerMask);
    else if (flag == localSingleColor) return 1;

    return 0; //flag == unoccupied
}

bool UnitigColors::write(ostream& stream_out) const {

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

bool UnitigColors::read(istream& stream_in) {

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

UnitigColors UnitigColors::reverse(const UnitigMapBase& um) const {

    const size_t len_unitig_km = um.size - Kmer::k + 1;

    UnitigColors new_cs;
    UnitigColors::const_iterator it = begin(), it_end = end();

    for (; it != it_end; ++it){

        const size_t color_id = it->getColorID(len_unitig_km);
        const size_t km_dist = it->getKmerPosition(len_unitig_km);
        const size_t new_km_dist = len_unitig_km - km_dist - 1;

        new_cs.add(len_unitig_km * color_id + new_km_dist);
    }

    return new_cs;
}

void UnitigColors::merge(const UnitigColors& cs){

    const uintptr_t flag = setBits & flagMask;
    const uintptr_t flag_cs = cs.setBits & flagMask;

    if (flag_cs == unoccupied) return;

    if (flag == flag_cs){

        if (flag == localBitVectorColor) setBits |= cs.setBits;
        else if (flag == ptrCompressedBitmap) *getPtrBitmap() |= *(cs.getConstPtrBitmap());
        else add(cs.setBits >> 2);
    }
    else if (flag_cs == ptrCompressedBitmap) {

        const Bitmap* bmp = cs.getConstPtrBitmap();

        for (Bitmap::const_iterator it = bmp->begin(), it_end = bmp->end(); it != it_end; ++it) add(*it);
    }
    else if (flag_cs == localBitVectorColor){

        uintptr_t setBits_tmp = cs.setBits >> 2;

        for (size_t it = 0; setBits_tmp != 0; ++it, setBits_tmp >>= 1){

            if ((setBits_tmp & 0x1) == 1) add(it);
        }
    }
    else add(cs.setBits >> 2); // flag_cs = localSingleColor
}

const size_t UnitigColors::maxBitVectorIDs = 62; // 64 bits - 2 bits for the color set type

const uintptr_t UnitigColors::ptrCompressedBitmap = 0x0;
const uintptr_t UnitigColors::localBitVectorColor = 0x1;
const uintptr_t UnitigColors::localSingleColor = 0x2;
const uintptr_t UnitigColors::unoccupied = 0x3;

const uintptr_t UnitigColors::flagMask = 0x3;
const uintptr_t UnitigColors::pointerMask = 0xfffffffffffffffc;
