#include "ColorSet.hpp"

UnitigColors::UnitigColors() : setBits(localBitVector) {}

UnitigColors::UnitigColors(const UnitigColors& o){

    const uintptr_t flag = o.setBits & flagMask;

    if (flag == ptrBitmap){

        setPointer = new Bitmap(*(o.getConstPtrBitmap()));
        setBits = (setBits & pointerMask) | ptrBitmap;
    }
    else if (flag == localTinyBitmap){

        t_bmp = o.t_bmp;
        setBits = (setBits & pointerMask) | localTinyBitmap;
    }
    else setBits = o.setBits;
}

UnitigColors::UnitigColors(UnitigColors&& o){

    setBits = o.setBits;
    o.setBits = localBitVector;
}

UnitigColors::~UnitigColors() {

    releaseMemory();
}

UnitigColors& UnitigColors::operator=(const UnitigColors& o){

    const uintptr_t flag = o.setBits & flagMask;

    releaseMemory();

    if (flag == ptrBitmap){

        setPointer = new Bitmap(*(o.getConstPtrBitmap()));
        setBits = (setBits & pointerMask) | ptrBitmap;
    }
    else if (flag == localTinyBitmap){

        t_bmp = o.t_bmp;
        setBits = (setBits & pointerMask) | localTinyBitmap;
    }
    else setBits = o.setBits;

    return *this;
}

UnitigColors& UnitigColors::operator=(UnitigColors&& o){

    if (this != &o) {

        releaseMemory();

        setBits = o.setBits;
        o.setBits = localBitVector;
    }

    return *this;
}

void UnitigColors::empty(){

    releaseMemory();
    setBits = localBitVector;
}

size_t UnitigColors::getSizeInBytes() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getPtrBitmap()->getSizeInBytes() + sizeof(Bitmap) + sizeof(UnitigColors);
    else if (flag == localTinyBitmap) return t_bmp.getSizeInBytes();
    return sizeof(UnitigColors);
}

void UnitigColors::add(const UnitigMapBase& um, const size_t color_id) {

    const size_t um_km_sz = um.size - Kmer::k + 1;
    size_t color_id_start = um_km_sz * color_id + um.dist;
    const size_t color_id_end = color_id_start + std::min(um_km_sz - um.dist, um.len);

    uintptr_t flag = setBits & flagMask;

    if (flag == localSingleInt){

        const uintptr_t setBits_tmp = setBits >> 2;

        if ((setBits_tmp < maxBitVectorIDs) && ((color_id_end - 1) < maxBitVectorIDs)){

            setBits = (1ULL << (setBits_tmp + 2)) | localBitVector;
        }
        else {

            t_bmp = TinyBitmap();

            if (t_bmp.add(setBits_tmp)) setBits = (setBits & pointerMask) | localTinyBitmap;
            else {

                t_bmp.empty();

                setPointer = new Bitmap;
                setPointer->add(setBits_tmp);

                setBits = (setBits & pointerMask) | ptrBitmap;
            }
        }

        flag = setBits & flagMask;
    }

    if (flag == localBitVector){

        if ((setBits == localBitVector) && (um.len == 1)) setBits = (color_id_start << 2) | localSingleInt;
        else if ((color_id_end - 1) < maxBitVectorIDs){

            for (; color_id_start != color_id_end; ++color_id_start) setBits |= 1ULL << (color_id_start + 2);
        }
        else {

            uintptr_t setBits_tmp_tb = setBits >> 2;
            uintptr_t setBits_tmp_cr = setBits >> 2;

            t_bmp = TinyBitmap();

            bool add_ok = true;

            for (size_t i = 0; (setBits_tmp_tb != 0) && add_ok; ++i, setBits_tmp_tb >>= 1) {

                if (setBits_tmp_tb & 0x1) add_ok = t_bmp.add(i);
            }

            if (add_ok) setBits = (setBits & pointerMask) | localTinyBitmap;
            else {

                t_bmp.empty();

                setPointer = new Bitmap;

                for (size_t i = 0; setBits_tmp_cr != 0; ++i, setBits_tmp_cr >>= 1) {

                    if (setBits_tmp_cr & 0x1) setPointer->add(i);
                }

                setBits = (setBits & pointerMask) | ptrBitmap;
            }
        }

        flag = setBits & flagMask;
    }

    if (flag == localTinyBitmap) {

        bool add_ok = true;

        while ((color_id_start < color_id_end) && add_ok) color_id_start += (add_ok = t_bmp.add(color_id_start));

        if (add_ok) setBits = (setBits & pointerMask) | localTinyBitmap;
        else {

            const size_t sz_t_bmp = t_bmp.size();

            uint32_t* values = new uint32_t[sz_t_bmp];

            size_t i = 0;

            for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it, ++i) values[i] = *it;

            t_bmp.empty();

            setPointer = new Bitmap;
            setPointer->addMany(sz_t_bmp, values);

            setBits = (setBits & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            delete[] values;
        }
    }

    if (flag == ptrBitmap) {

        Bitmap* bitmap = getPtrBitmap();

        for (; color_id_start < color_id_end; ++color_id_start) bitmap->add(color_id_start);

        bitmap->runOptimize();
    }
}

void UnitigColors::add(const size_t color_id) { // PRIVATE

    uintptr_t flag = setBits & flagMask;

    if (flag == localSingleInt){

        const uintptr_t setBits_tmp = setBits >> 2;

        if ((setBits_tmp < maxBitVectorIDs) && (color_id < maxBitVectorIDs)){

            setBits = (1ULL << (setBits_tmp + 2)) | (1ULL << (color_id + 2)) | localBitVector;
        }
        else {

            t_bmp = TinyBitmap();

            if (t_bmp.add(setBits_tmp)) setBits = (setBits & pointerMask) | localTinyBitmap;
            else {

                t_bmp.empty();

                setPointer = new Bitmap;
                setPointer->add(setBits_tmp);

                setBits = (setBits & pointerMask) | ptrBitmap;
            }
        }

        flag = setBits & flagMask;
    }

    if (flag == localBitVector){

        if (setBits == localBitVector) setBits = (color_id << 2) | localSingleInt;
        else if (color_id < maxBitVectorIDs) setBits |= 1ULL << (color_id + 2);
        else {

            uintptr_t setBits_tmp_tb = setBits >> 2;
            uintptr_t setBits_tmp_cr = setBits >> 2;

            t_bmp = TinyBitmap();

            bool add_ok = true;

            for (size_t i = 0; (setBits_tmp_tb != 0) && add_ok; ++i, setBits_tmp_tb >>= 1) {

                if (setBits_tmp_tb & 0x1) add_ok = t_bmp.add(i);
            }

            if (add_ok) setBits = (setBits & pointerMask) | localTinyBitmap;
            else {

                t_bmp.empty();

                setPointer = new Bitmap;

                for (size_t i = 0; setBits_tmp_cr != 0; ++i, setBits_tmp_cr >>= 1) {

                    if (setBits_tmp_cr & 0x1) setPointer->add(i);
                }

                setBits = (setBits & pointerMask) | ptrBitmap;
            }
        }

        flag = setBits & flagMask;
    }

    if (flag == localTinyBitmap) {

        bool add_ok = t_bmp.add(color_id);

        if (add_ok) setBits = (setBits & pointerMask) | localTinyBitmap;
        else {

            const size_t sz_t_bmp = t_bmp.size();

            uint32_t* values = new uint32_t[sz_t_bmp];

            size_t i = 0;

            for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it, ++i) values[i] = *it;

            t_bmp.empty();

            setPointer = new Bitmap;
            setPointer->addMany(sz_t_bmp, values);

            setBits = (setBits & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            delete[] values;
        }
    }

    if (flag == ptrBitmap) getPtrBitmap()->add(color_id); // flag == ptrBitmap
}

void UnitigColors::remove(const UnitigMapBase& um, const size_t color_id) {

    uintptr_t flag = setBits & flagMask;

    const size_t um_km_sz = um.size - Kmer::k + 1;
    size_t color_id_start = um_km_sz * color_id + um.dist;
    const size_t color_id_end = color_id_start + std::min(um_km_sz - um.dist, um.len);

    if (flag == localBitVector){

        uintptr_t mask = 0;

        for (; color_id_start < min(maxBitVectorIDs, color_id_end); ++color_id_start) mask |= 1ULL << (color_id_start + 2);

        setBits &= ~mask;
    }
    else if (flag == localSingleInt){

        const uintptr_t color_id = setBits >> 2;

        if ((color_id >= color_id_start) && (color_id < color_id_end)) setBits = localBitVector;
    }
    else if (flag == localTinyBitmap){

        bool rm_ok = true;

        while ((color_id_start < color_id_end) && rm_ok) color_id_start += (rm_ok = t_bmp.remove(color_id_start));

        if (rm_ok){

            const size_t card = t_bmp.size();

            if (card == 0) empty();
            else if (card == 1){

                const uint32_t color_id = *(t_bmp.begin());

                empty();
                add(color_id);
            }
            else if ((card <= maxBitVectorIDs) && (t_bmp.maximum() < maxBitVectorIDs)){

                UnitigColors new_uc;

                for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) new_uc.add(*it);

                *this = move(new_uc);
            }
            else setBits = (setBits & pointerMask) | localTinyBitmap;
        }
        else {

            const size_t sz_t_bmp = t_bmp.size();

            uint32_t* values = new uint32_t[sz_t_bmp];

            size_t i = 0;

            for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it, ++i) values[i] = *it;

            t_bmp.empty();

            setPointer = new Bitmap;
            setPointer->addMany(sz_t_bmp, values);

            setBits = (setBits & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            delete[] values;
        }
    }

    if (flag == ptrBitmap) {

        Bitmap* bitmap = getPtrBitmap();

        for (; color_id_start < color_id_end; ++color_id_start) bitmap->remove(color_id_start);

        const size_t card = bitmap->cardinality();

        if (card == 0) empty();
        else if (card == 1){

            uint32_t color_id;

            bitmap->select(0, &color_id);

            empty();
            add(color_id);
        }
        else if ((card <= maxBitVectorIDs) && (bitmap->maximum() < maxBitVectorIDs)){

            UnitigColors new_uc;

            const_iterator it = begin(), it_end = end();

            for (; it != it_end; ++it) new_uc.add(it->ck_id);

            *this = move(new_uc);
        }
        else if ((setBits & flagMask) == ptrBitmap) bitmap->runOptimize();
    }
}

bool UnitigColors::contains(const UnitigMapBase& um, const size_t color_id) const {

    const uintptr_t flag = setBits & flagMask;
    const size_t um_km_sz = um.size - Kmer::k + 1;
    size_t color_id_start = um_km_sz * color_id + um.dist;
    const size_t color_id_end = color_id_start + std::min(um_km_sz - um.dist, um.len);

    if (flag == ptrBitmap){

        const Bitmap* bitmap = getConstPtrBitmap();

        for (; color_id_start < color_id_end; ++color_id_start){

            if (!bitmap->contains(color_id_start)) return false;
        }
    }

    if (flag == localTinyBitmap){

        for (; color_id_start < color_id_end; ++color_id_start){

            if (!t_bmp.contains(color_id_start)) return false;
        }
    }

    if (flag == localBitVector){

        if ((color_id_end - 1) < maxBitVectorIDs){

            uintptr_t setBits_tmp = setBits >> (color_id_start + 2);

            for (; color_id_start < color_id_end; ++color_id_start, setBits_tmp >>= 1){

                if ((setBits_tmp & 0x1) == 0) return false;
            }
        }
        else return false;
    }

    if (flag == localSingleInt) return (um.len == 1) && (color_id_start == (setBits >> 2));

    return true;
}

bool UnitigColors::contains(const size_t color_km_id) const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->contains(color_km_id);
    if (flag == localTinyBitmap) return t_bmp.contains(color_km_id);
    if (flag == localSingleInt) return (color_km_id == (setBits >> 2));

    if (color_km_id < maxBitVectorIDs){

        const uintptr_t setBits_tmp = 0x1ULL << (color_km_id + 2);

        return ((setBits & setBits_tmp) != 0);
    }

    return false;
}

UnitigColors::ColorKmer_ID UnitigColors::maximum() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return ColorKmer_ID(getConstPtrBitmap()->maximum());
    if (flag == localTinyBitmap) return ColorKmer_ID(t_bmp.maximum());
    if (flag == localSingleInt) return ColorKmer_ID(setBits >> 2);

    return ColorKmer_ID(maxBitVectorIDs - __builtin_clzll(setBits | ~pointerMask) - 1);
}

size_t UnitigColors::size() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->cardinality();
    if (flag == localTinyBitmap) return t_bmp.size();
    if (flag == localBitVector) return __builtin_popcountll(setBits & pointerMask);

    return 1;
}

size_t UnitigColors::size(const UnitigMapBase& um, const size_t color_id) const {

    const size_t length_unitig_km = um.size - Kmer::k + 1;
    const size_t start_pos = color_id * length_unitig_km;

    size_t end_pos = start_pos + length_unitig_km;

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) {

        const Bitmap* bmp = getConstPtrBitmap();

        if (start_pos == 0) return bmp->rank(end_pos - 1);

        return bmp->rank(end_pos - 1) - bmp->rank(start_pos - 1);
    }

    if (flag == localTinyBitmap) return t_bmp.size(start_pos, end_pos);

    if (flag == localBitVector){

        if (start_pos < maxBitVectorIDs){

            uintptr_t mask = 0;

            end_pos = min(end_pos, maxBitVectorIDs);

            for (size_t i = start_pos; i != end_pos; ++i) mask |= 0x1ULL << (i + 2);

            return __builtin_popcountll(setBits & mask);
        }
    }
    else {

        const uintptr_t ck_id = setBits >> 2;

        if ((ck_id >= start_pos) && (ck_id < end_pos)) return 1;
    }

    return 0;
}

bool UnitigColors::write(ostream& stream_out) const {

    if (stream_out.good()){

        const uintptr_t flag = setBits & flagMask;

        if (flag == ptrBitmap){

            const uint32_t expected_sz = getConstPtrBitmap()->getSizeInBytes();

            const uintptr_t flag_expected_sz = (static_cast<uintptr_t>(expected_sz) << 32) | flag;

            char* serialized = new char[expected_sz];

            getConstPtrBitmap()->write(serialized);

            stream_out.write(reinterpret_cast<const char*>(&flag_expected_sz), sizeof(uintptr_t));
            stream_out.write(serialized, expected_sz);

            delete[] serialized;
        }
        else if (flag == localTinyBitmap){

            stream_out.write(reinterpret_cast<const char*>(&flag), sizeof(uintptr_t));
            t_bmp.write(stream_out);
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

        if (flag == ptrBitmap){

            const uint32_t expected_sz = static_cast<uint32_t>(setBits >> 32);

            char* serialized = new char[expected_sz];
            setPointer = new Bitmap;

            stream_in.read(serialized, expected_sz);

            *setPointer = std::move(Bitmap::read(serialized));

            setBits = (setBits & pointerMask) | ptrBitmap;

            delete[] serialized;
        }
        else if (flag == localTinyBitmap){

            t_bmp = TinyBitmap();

            t_bmp.read(stream_in);

            setBits = (setBits & pointerMask) | localTinyBitmap;
        }

        return true;
    }

    return false;
}

UnitigColors::const_iterator UnitigColors::begin() const {

    const_iterator it(this, true);
    return ++it;
}

UnitigColors::const_iterator UnitigColors::end() const {

    return const_iterator(this, false);
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

UnitigColors::ColorKmer_ID::ColorKmer_ID() : ck_id(0xffffffffffffffff) {}

UnitigColors::ColorKmer_ID::ColorKmer_ID(const size_t id) : ck_id(id) {}

UnitigColors::ColorKmer_ID& UnitigColors::ColorKmer_ID::operator=(const size_t ck_id_){

    ck_id = ck_id_;

    return *this;
}

size_t UnitigColors::ColorKmer_ID::getColorID(const size_t length_unitig_km) const {

    if (ck_id == 0xffffffffffffffff){

        cerr << "ColorKmer_ID:getColorID(): Invalid ColorKmer_ID" << endl;
        return ck_id;
    }

    return (ck_id / length_unitig_km);
}

size_t UnitigColors::ColorKmer_ID::getKmerPosition(const size_t length_unitig_km) const {

    if (ck_id == 0xffffffffffffffff){

        cerr << "ColorKmer_ID:getKmerPosition(): Invalid ColorKmer_ID" << endl;
        return ck_id;
    }

    return (ck_id % length_unitig_km);
}




UnitigColors::UnitigColors_const_iterator::UnitigColors_const_iterator() :  cs(nullptr), flag(localBitVector), it_setBits(0),
                                                                            cs_sz(0), it_roar(empty_roar.end()) {}

UnitigColors::UnitigColors_const_iterator::UnitigColors_const_iterator(const UnitigColors* cs_, const bool beg) : cs(cs_),
                                                                        it_setBits(0xffffffffffffffff), ck_id(0), it_roar(empty_roar.end()) {

    flag = cs->setBits & flagMask;

    if (flag == ptrBitmap){

        it_roar = beg ? cs->getConstPtrBitmap()->begin() : cs->getConstPtrBitmap()->end();
        cs_sz = cs->size();
    }
    else if (flag == localTinyBitmap){

        it_t_bmp = beg ? cs->t_bmp.begin() : cs->t_bmp.end();
        cs_sz = cs->size();
    }
    else cs_sz = (flag == localSingleInt) ? 1 : maxBitVectorIDs;

    if (!beg) it_setBits = cs_sz;
}

UnitigColors::UnitigColors_const_iterator& UnitigColors::UnitigColors_const_iterator::operator=(const UnitigColors_const_iterator& o) {

    cs = o.cs;
    flag = o.flag;
    it_setBits = o.it_setBits;
    ck_id = o.ck_id;
    cs_sz = o.cs_sz;
    it_roar = o.it_roar;
    it_t_bmp = o.it_t_bmp;

    return *this;
}

const UnitigColors::ColorKmer_ID& UnitigColors::UnitigColors_const_iterator::operator*() const { return ck_id; }

const UnitigColors::ColorKmer_ID* UnitigColors::UnitigColors_const_iterator::operator->() const { return &ck_id; }

UnitigColors::UnitigColors_const_iterator UnitigColors::UnitigColors_const_iterator::operator++(int) {

    UnitigColors_const_iterator tmp(*this);
    operator++();
    return tmp;
}

UnitigColors::UnitigColors_const_iterator& UnitigColors::UnitigColors_const_iterator::operator++() {

    if (it_setBits == cs_sz) return *this;

    ++it_setBits;

    if (flag == ptrBitmap) {

        if (it_setBits != 0) ++it_roar;
        if (it_roar != cs->getConstPtrBitmap()->end()) ck_id = *it_roar;
    }
    else if (flag == localTinyBitmap) {

        if (it_setBits != 0) ++it_t_bmp;
        if (it_t_bmp != cs->t_bmp.end()) ck_id = *it_t_bmp;
    }
    else if (flag == localBitVector){

        while (it_setBits < maxBitVectorIDs){

            if (((cs->setBits >> (it_setBits + 2)) & 0x1) != 0){

                ck_id = it_setBits;
                break;
            }

            ++it_setBits;
        }
    }
    else if (flag == localSingleInt) ck_id = cs->setBits >> 2;

    return *this;
}

UnitigColors::UnitigColors_const_iterator& UnitigColors::UnitigColors_const_iterator::nextColor(const size_t length_unitig_km) {

    if (it_setBits == cs_sz) return *this;

    const size_t nextPos = (ck_id.getColorID(length_unitig_km) + 1) * length_unitig_km;

    if (flag == ptrBitmap) {

        it_roar.equalorlarger(nextPos);

        if (it_roar != cs->getConstPtrBitmap()->end()) ck_id = *it_roar;
        else it_setBits = cs_sz;
    }
    else if (flag == localTinyBitmap) {

        while ((it_t_bmp != cs->t_bmp.end()) && (*it_t_bmp < nextPos)) ++it_t_bmp;

        if (it_t_bmp != cs->t_bmp.end()) ck_id = *it_t_bmp;
        else it_setBits = cs_sz;
    }
    else if (flag == localBitVector){

        it_setBits = min(maxBitVectorIDs, nextPos);

        while (it_setBits < maxBitVectorIDs){

            if (((cs->setBits >> (it_setBits + 2)) & 0x1) != 0){

                ck_id = it_setBits;
                break;
            }

            ++it_setBits;
        }
    }
    else ++it_setBits;

    return *this;
}

bool UnitigColors::UnitigColors_const_iterator::operator==(const UnitigColors_const_iterator& o) const {

    return  (cs == o.cs) && (flag == o.flag) && (cs_sz == o.cs_sz) &&
            ((flag == ptrBitmap) ? (it_roar == o.it_roar) : ((flag == localTinyBitmap) ? (it_t_bmp == o.it_t_bmp) : (it_setBits == o.it_setBits)));
}

bool UnitigColors::UnitigColors_const_iterator::operator!=(const UnitigColors_const_iterator& o) const {

    return !operator==(o);
}

const size_t UnitigColors::maxBitVectorIDs = 62; // 64 bits - 2 bits for the color set type

const uintptr_t UnitigColors::localTinyBitmap = 0x0;
const uintptr_t UnitigColors::localBitVector = 0x1;
const uintptr_t UnitigColors::localSingleInt = 0x2;
const uintptr_t UnitigColors::ptrBitmap = 0x3;

const uintptr_t UnitigColors::flagMask = 0x3;
const uintptr_t UnitigColors::pointerMask = 0xFFFFFFFFFFFFFFFC;
