#include "ColorSet.hpp"

UnitigColors::UnitigColors() : setBits(localBitVector) {}

UnitigColors::UnitigColors(const UnitigColors& o) {

    const uintptr_t flag = o.setBits & flagMask;

    if (flag == ptrUnitigColors){

        const UnitigColors* uc_o = o.getConstPtrUnitigColors();

        UnitigColors* setPtrUC = new UnitigColors[2];

        setPtrUC[0] = uc_o[0];
        setPtrUC[1] = uc_o[1];

        setBits = (reinterpret_cast<uintptr_t>(setPtrUC) & pointerMask) | ptrUnitigColors;
    }
    else if (flag == ptrBitmap){

        Bitmap* setPtrBmp = new Bitmap;

        setPtrBmp->r = o.getConstPtrBitmap()->r;

        setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
    }
    else if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = o.getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        TinyBitmap t_bmp_cpy(t_bmp);

        t_bmp.detach();

        setBits = (reinterpret_cast<uintptr_t>(t_bmp_cpy.detach()) & pointerMask) | localTinyBitmap;
    }
    else setBits = o.setBits;
}

UnitigColors::UnitigColors(UnitigColors&& o) : setBits(o.setBits) {

    o.setBits = localBitVector;
}

UnitigColors::~UnitigColors() {

    releaseMemory();
}

UnitigColors& UnitigColors::operator=(const UnitigColors& o){

    if (this != &o) {

        const uintptr_t flag = o.setBits & flagMask;
        const uintptr_t this_flag = setBits & flagMask;

        if (flag == ptrUnitigColors){

            const UnitigColors* uc_o = o.getConstPtrUnitigColors();

            UnitigColors* setPtrUC = nullptr;

            if (this_flag == ptrUnitigColors) setPtrUC = getPtrUnitigColors();
            else {

                releaseMemory();
                setPtrUC = new UnitigColors[2];
            }

            setPtrUC[0] = uc_o[0];
            setPtrUC[1] = uc_o[1];

            setBits = (reinterpret_cast<uintptr_t>(setPtrUC) & pointerMask) | ptrUnitigColors;
        }
        else if (flag == ptrBitmap){

            Bitmap* setPtrBmp = nullptr;

            if (this_flag == ptrBitmap) setPtrBmp = getPtrBitmap();
            else {

                releaseMemory();
                setPtrBmp = new Bitmap;
            }

            setPtrBmp->r = o.getConstPtrBitmap()->r;

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
        }
        else if (flag == localTinyBitmap){

            releaseMemory();

            uint16_t* setPtrTinyBmp = o.getPtrTinyBitmap();
            TinyBitmap t_bmp(&setPtrTinyBmp);
            TinyBitmap t_bmp_cpy(t_bmp);

            t_bmp.detach();

            setBits = (reinterpret_cast<uintptr_t>(t_bmp_cpy.detach()) & pointerMask) | localTinyBitmap;
        }
        else setBits = o.setBits;
    }

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

    if (flag == ptrUnitigColors){

        return sizeof(UnitigColors) + getConstPtrUnitigColors()[0].getSizeInBytes() + getConstPtrUnitigColors()[1].getSizeInBytes();
    }

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.getSizeInBytes() + sizeof(Bitmap) + sizeof(UnitigColors);

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const size_t ret = t_bmp.getSizeInBytes();

        t_bmp.detach();

        return ret;
    }

    return sizeof(UnitigColors);
}

void UnitigColors::add(const UnitigMapBase& um, const size_t color_id) {

    uintptr_t flag = setBits & flagMask;

    if (flag == ptrUnitigColors){

        UnitigColors* uc = getPtrUnitigColors();

        if (!uc[0].contains(color_id)) uc[1].add(um, color_id);
    }
    else {

        const size_t um_km_sz = um.size - Kmer::k + 1;
        size_t color_id_start = um_km_sz * color_id + um.dist;
        const size_t color_id_end = color_id_start + std::min(um_km_sz - um.dist, um.len);

        if (flag == localSingleInt){

            const uintptr_t setBits_tmp = setBits >> shiftMaskBits;

            if ((setBits_tmp < maxBitVectorIDs) && ((color_id_end - 1) < maxBitVectorIDs)){

                setBits = (1ULL << (setBits_tmp + shiftMaskBits)) | localBitVector;
            }
            else {

                TinyBitmap t_bmp;

                if (t_bmp.add(setBits_tmp)) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
                else {

                    Bitmap* setPtrBmp = new Bitmap;

                    t_bmp.empty();
                    setPtrBmp->r.add(setBits_tmp);

                    setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
                }
            }

            flag = setBits & flagMask;
        }

        if (flag == localBitVector){

            if ((setBits == localBitVector) && (um.len == 1)) setBits = (color_id_start << shiftMaskBits) | localSingleInt;
            else if ((color_id_end - 1) < maxBitVectorIDs){

                for (; color_id_start != color_id_end; ++color_id_start) setBits |= 1ULL << (color_id_start + shiftMaskBits);
            }
            else {

                uintptr_t setBits_tmp_tb = setBits >> shiftMaskBits;
                uintptr_t setBits_tmp_cr = setBits >> shiftMaskBits;

                TinyBitmap t_bmp;

                bool add_ok = true;

                for (size_t i = 0; (setBits_tmp_tb != 0) && add_ok; ++i, setBits_tmp_tb >>= 1) {

                    if (setBits_tmp_tb & 0x1) add_ok = t_bmp.add(i);
                }

                if (add_ok) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
                else {

                    Bitmap* setPtrBmp = new Bitmap;

                    t_bmp.empty();

                    for (size_t i = 0; setBits_tmp_cr != 0; ++i, setBits_tmp_cr >>= 1) {

                        if (setBits_tmp_cr & 0x1) setPtrBmp->r.add(i);
                    }

                    setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
                }
            }

            flag = setBits & flagMask;
        }

        if (flag == localTinyBitmap) {

            bool add_ok = true;

            uint16_t* setPtrTinyBmp = getPtrTinyBitmap();

            TinyBitmap t_bmp(&setPtrTinyBmp);

            while ((color_id_start < color_id_end) && add_ok) color_id_start += (add_ok = t_bmp.add(color_id_start));

            if (add_ok) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
            else {

                const size_t sz_t_bmp = t_bmp.size();

                uint32_t* values = new uint32_t[sz_t_bmp];

                size_t i = 0;

                Bitmap* setPtrBmp = new Bitmap;

                for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it, ++i) values[i] = *it;

                t_bmp.empty();
                setPtrBmp->r.addMany(sz_t_bmp, values);

                setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
                flag = ptrBitmap;

                delete[] values;
            }
        }

        if (flag == ptrBitmap) {

            Bitmap* bitmap = getPtrBitmap();

            for (; color_id_start < color_id_end; ++color_id_start) bitmap->r.add(color_id_start);

            bitmap->r.runOptimize();
        }
    }
}

void UnitigColors::add(const size_t color_id) { // PRIVATE

    uintptr_t flag = setBits & flagMask;

    if (flag == ptrUnitigColors){

        UnitigColors* uc = getPtrUnitigColors();

        if (!uc[0].contains(color_id)) uc[1].add(color_id);
    }
    else {

        if (flag == localSingleInt){

            const uintptr_t setBits_tmp = setBits >> shiftMaskBits;

            if ((setBits_tmp < maxBitVectorIDs) && (color_id < maxBitVectorIDs)){

                setBits = (1ULL << (setBits_tmp + shiftMaskBits)) | (1ULL << (color_id + shiftMaskBits)) | localBitVector;
            }
            else {

                TinyBitmap t_bmp;

                if (t_bmp.add(setBits_tmp)) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
                else {

                    Bitmap* setPtrBmp = new Bitmap;

                    t_bmp.empty();
                    setPtrBmp->r.add(setBits_tmp);

                    setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
                }
            }

            flag = setBits & flagMask;
        }

        if (flag == localBitVector){

            if (setBits == localBitVector) setBits = (color_id << shiftMaskBits) | localSingleInt;
            else if (color_id < maxBitVectorIDs) setBits |= 1ULL << (color_id + shiftMaskBits);
            else {

                uintptr_t setBits_tmp_tb = setBits >> shiftMaskBits;
                uintptr_t setBits_tmp_cr = setBits >> shiftMaskBits;

                TinyBitmap t_bmp;

                bool add_ok = true;

                for (size_t i = 0; (setBits_tmp_tb != 0) && add_ok; ++i, setBits_tmp_tb >>= 1) {

                    if (setBits_tmp_tb & 0x1) add_ok = t_bmp.add(i);
                }

                if (add_ok) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
                else {

                    Bitmap* setPtrBmp = new Bitmap;

                    t_bmp.empty();

                    for (size_t i = 0; setBits_tmp_cr != 0; ++i, setBits_tmp_cr >>= 1) {

                        if (setBits_tmp_cr & 0x1) setPtrBmp->r.add(i);
                    }

                    setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
                }
            }

            flag = setBits & flagMask;
        }

        if (flag == localTinyBitmap) {

            uint16_t* setPtrTinyBmp = getPtrTinyBitmap();

            TinyBitmap t_bmp(&setPtrTinyBmp);

            if (t_bmp.add(color_id)) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
            else {

                const size_t sz_t_bmp = t_bmp.size();

                uint32_t* values = new uint32_t[sz_t_bmp];

                size_t i = 0;

                Bitmap* setPtrBmp = new Bitmap;

                for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it, ++i) values[i] = *it;

                t_bmp.empty();
                setPtrBmp->r.addMany(sz_t_bmp, values);

                setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
                flag = ptrBitmap;

                delete[] values;
            }
        }

        if (flag == ptrBitmap) getPtrBitmap()->r.add(color_id); // flag == ptrBitmap
    }
}

void UnitigColors::remove(const UnitigMapBase& um, const size_t color_id) {

    uintptr_t flag = setBits & flagMask;

    const size_t um_km_sz = um.size - Kmer::k + 1;

    if (flag == ptrUnitigColors){

        UnitigColors* uc = getPtrUnitigColors();

        if (uc[0].contains(color_id)){

            const UnitigMapBase fake_um(0, 1, Kmer::k, true);

            uc[0].remove(fake_um, color_id);

            if (um.dist != 0){

                const UnitigMapBase fake_um_start(0, um.dist, um.size, um.strand);

                uc[1].add(fake_um_start, color_id);
            }

            if ((um_km_sz - um.len - um.dist) != 0){

                const UnitigMapBase fake_um_end(um.dist + um.len, um_km_sz - um.len - um.dist, um.size, um.strand);

                uc[1].add(fake_um_end, color_id);
            }
        }
        else uc[1].remove(um, color_id);

        return;
    }

    size_t color_id_start = um_km_sz * color_id + um.dist;
    const size_t color_id_end = color_id_start + std::min(um_km_sz - um.dist, um.len);

    if (flag == localBitVector){

        uintptr_t mask = 0;

        for (; color_id_start < min(maxBitVectorIDs, color_id_end); ++color_id_start) mask |= 1ULL << (color_id_start + shiftMaskBits);

        setBits &= ~mask;
    }
    else if (flag == localSingleInt){

        const uintptr_t color_id = setBits >> shiftMaskBits;

        if ((color_id >= color_id_start) && (color_id < color_id_end)) setBits = localBitVector;
    }
    else if (flag == localTinyBitmap){

        bool rm_ok = true;

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();

        TinyBitmap t_bmp(&setPtrTinyBmp);

        while ((color_id_start < color_id_end) && rm_ok) color_id_start += (rm_ok = t_bmp.remove(color_id_start));

        if (rm_ok){

            const size_t card = t_bmp.size();

            if (card == 0){

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
                empty();
            }
            else if (card == 1){

                const uint32_t color_id = *(t_bmp.begin());

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;

                empty();
                add(color_id);
            }
            else if ((card <= maxBitVectorIDs) && (t_bmp.maximum() < maxBitVectorIDs)){

                UnitigColors new_uc;

                for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) new_uc.add(*it);

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;

                *this = move(new_uc);
            }
            else setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        }
        else {

            const size_t sz_t_bmp = t_bmp.size();

            uint32_t* values = new uint32_t[sz_t_bmp];

            size_t i = 0;

            Bitmap* setPtrBmp = new Bitmap;

            for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it, ++i) values[i] = *it;

            t_bmp.empty();
            setPtrBmp->r.addMany(sz_t_bmp, values);

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            delete[] values;
        }
    }

    if (flag == ptrBitmap) {

        Bitmap* bitmap = getPtrBitmap();

        for (; color_id_start < color_id_end; ++color_id_start) bitmap->r.remove(color_id_start);

        const size_t card = bitmap->r.cardinality();

        if (card == 0) empty();
        else if (card == 1){

            uint32_t color_id;

            bitmap->r.select(0, &color_id);

            empty();
            add(color_id);
        }
        else if ((card <= maxBitVectorIDs) && (bitmap->r.maximum() < maxBitVectorIDs)){

            UnitigColors new_uc;

            const_iterator it = begin(um), it_end = end();

            for (; it != it_end; ++it) new_uc.add(it.getColorID() * um_km_sz + it.getKmerPosition());

            *this = move(new_uc);
        }
        else if ((setBits & flagMask) == ptrBitmap) bitmap->r.runOptimize();
    }
}

bool UnitigColors::contains(const UnitigMapBase& um, const size_t color_id) const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrUnitigColors){

        const UnitigColors* uc = getConstPtrUnitigColors();

        return (uc[0].contains(color_id) || uc[1].contains(um, color_id));
    }

    const size_t um_km_sz = um.size - Kmer::k + 1;
    size_t color_id_start = um_km_sz * color_id + um.dist;
    const size_t color_id_end = color_id_start + std::min(um_km_sz - um.dist, um.len);

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.containsRange(color_id_start, color_id_end - 1);

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const bool ret = t_bmp.containsRange(color_id_start, color_id_end - 1);

        t_bmp.detach();

        return ret;
    }

    if (flag == localBitVector){

        if ((color_id_end - 1) < maxBitVectorIDs){

            uintptr_t setBits_tmp = setBits >> (color_id_start + shiftMaskBits);

            for (; color_id_start < color_id_end; ++color_id_start, setBits_tmp >>= 1){

                if ((setBits_tmp & 0x1) == 0) return false;
            }
        }
        else return false;
    }

    if (flag == localSingleInt) return (um.len == 1) && (color_id_start == (setBits >> shiftMaskBits));

    return true;
}

bool UnitigColors::contains(const size_t color_km_id) const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrUnitigColors){

        const UnitigColors* uc = getConstPtrUnitigColors();

        return (uc[0].contains(color_km_id) || uc[1].contains(color_km_id));
    }

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.contains(color_km_id);

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const bool ret = t_bmp.contains(color_km_id);

        t_bmp.detach();

        return ret;
    }

    if (flag == localSingleInt) return (color_km_id == (setBits >> shiftMaskBits));

    if (color_km_id < maxBitVectorIDs){

        const uintptr_t setBits_tmp = 0x1ULL << (color_km_id + shiftMaskBits);

        return ((setBits & setBits_tmp) != 0);
    }

    return false;
}

size_t UnitigColors::colorMax(const UnitigMapBase& um) const {

    const uintptr_t flag = setBits & flagMask;
    const size_t length_unitig_km = um.size - Kmer::k + 1;

    if (flag == ptrUnitigColors){

        const UnitigColors* uc = getConstPtrUnitigColors();
        const UnitigMapBase fake_um(0, 1, Kmer::k, um.strand);

        return max(uc[0].colorMax(fake_um), uc[1].colorMax(um));
    }

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.maximum() / length_unitig_km;

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const size_t ret = t_bmp.maximum() / length_unitig_km;

        t_bmp.detach();

        return ret;
    }

    if (flag == localSingleInt) return (setBits >> shiftMaskBits) / length_unitig_km;

    const int nb_lead_0 = __builtin_clzll(setBits | flagMask);

    return (maxBitVectorIDs - nb_lead_0 - (nb_lead_0 != maxBitVectorIDs)) / length_unitig_km;
}

size_t UnitigColors::size(const UnitigMapBase& um) const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrUnitigColors){

        const UnitigColors* uc = getConstPtrUnitigColors();

        return uc[0].size() * (um.size - Kmer::k + 1) + uc[1].size();
    }

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.cardinality();

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const size_t ret = t_bmp.size();

        t_bmp.detach();

        return ret;
    }

    if (flag == localBitVector) return __builtin_popcountll(setBits & pointerMask);

    return 1;
}

size_t UnitigColors::size(const UnitigMapBase& um, const size_t color_id) const {

    const size_t length_unitig_km = um.size - Kmer::k + 1;

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrUnitigColors){

        const UnitigColors* uc = getConstPtrUnitigColors();

        return static_cast<size_t>(uc[0].contains(color_id)) * length_unitig_km + uc[1].size(um, color_id);
    }

    const size_t start_pos = color_id * length_unitig_km;

    size_t end_pos = start_pos + length_unitig_km;

    if (flag == ptrBitmap) {

        const Bitmap* bmp = getConstPtrBitmap();

        if (start_pos == 0) return bmp->r.rank(end_pos - 1);

        return bmp->r.rank(end_pos - 1) - bmp->r.rank(start_pos - 1);
    }

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const size_t ret = t_bmp.size(start_pos, end_pos);

        t_bmp.detach();

        return ret;
    }

    if (flag == localBitVector){

        if (start_pos < maxBitVectorIDs){

            uintptr_t mask = 0;

            end_pos = min(end_pos, maxBitVectorIDs);

            for (size_t i = start_pos; i != end_pos; ++i) mask |= 0x1ULL << (i + shiftMaskBits);

            return __builtin_popcountll(setBits & mask);
        }
    }
    else {

        const uintptr_t ck_id = setBits >> shiftMaskBits;

        if ((ck_id >= start_pos) && (ck_id < end_pos)) return 1;
    }

    return 0;
}

size_t UnitigColors::size() const { //Private

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrUnitigColors){

        const UnitigColors* uc = getConstPtrUnitigColors();

        return uc[0].size() + uc[1].size();
    }

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.cardinality();

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const size_t ret = t_bmp.size();

        t_bmp.detach();

        return ret;
    }

    if (flag == localBitVector) return __builtin_popcountll(setBits & pointerMask);

    return 1;
}

bool UnitigColors::optimizeFullColors(const UnitigMapBase& um, const size_t color_start){

    const uintptr_t flag = setBits & flagMask;

    if ((um.size > Kmer::k) && (flag != localBitVector) && (flag != localSingleInt)){

        const size_t color_maximum = colorMax(um);

        UnitigColors* uc = nullptr;

        UnitigColors new_uc;

        if (flag == ptrUnitigColors){

            const UnitigColors* this_uc = getConstPtrUnitigColors();

            for (size_t color = color_start; color <= color_maximum; ++color){

                if (!this_uc[0].contains(color) && this_uc[1].contains(um, color)){

                    if (uc == nullptr){

                        new_uc = *this;
                        uc = new_uc.getPtrUnitigColors();
                    }

                    uc[0].add(color);
                    uc[1].remove(um, color);
                }
            }
        }
        else {

            for (size_t color = color_start; color <= color_maximum; ++color){

                if (contains(um, color)){

                    if (uc == nullptr){

                        uc = new UnitigColors[2];
                        uc[1] = *this;
                    }

                    uc[0].add(color);
                    uc[1].remove(um, color);
                }
            }
        }

        if (uc != nullptr){

            new_uc.setBits = (reinterpret_cast<uintptr_t>(uc) & pointerMask) | ptrUnitigColors;

            if (new_uc.getSizeInBytes() < getSizeInBytes()){

                *this = move(new_uc);
                return true;
            }
        }
    }

    return false;
}

bool UnitigColors::write(ostream& stream_out) const {

    if (stream_out.good()){

        const uintptr_t flag = setBits & flagMask;

        if (flag == ptrUnitigColors){

            stream_out.write(reinterpret_cast<const char*>(&flag), sizeof(uintptr_t));

            const UnitigColors* uc = getConstPtrUnitigColors();

            bool ret = uc[0].write(stream_out);

            return (ret ? uc[1].write(stream_out) : ret);
        }
        else if (flag == ptrBitmap){

            const uint32_t expected_sz = getConstPtrBitmap()->r.getSizeInBytes();

            const uintptr_t flag_expected_sz = (static_cast<uintptr_t>(expected_sz) << 32) | flag;

            char* serialized = new char[expected_sz];

            getConstPtrBitmap()->r.write(serialized);

            stream_out.write(reinterpret_cast<const char*>(&flag_expected_sz), sizeof(uintptr_t));
            stream_out.write(serialized, expected_sz);

            delete[] serialized;
        }
        else if (flag == localTinyBitmap){

            uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
            TinyBitmap t_bmp(&setPtrTinyBmp);

            stream_out.write(reinterpret_cast<const char*>(&flag), sizeof(uintptr_t));

            t_bmp.write(stream_out);
            t_bmp.detach();
        }
        else stream_out.write(reinterpret_cast<const char*>(&setBits), sizeof(uintptr_t));

        return true;
    }

    return false;
}

bool UnitigColors::read(istream& stream_in) {

    if (stream_in.good()){

        empty();

        stream_in.read(reinterpret_cast<char*>(&setBits), sizeof(uintptr_t));

        const uintptr_t flag = setBits & flagMask;

        if (flag == ptrUnitigColors){

            UnitigColors* setPtrUC = new UnitigColors[2];

            setPtrUC[0].read(stream_in);
            setPtrUC[1].read(stream_in);

            setBits = (reinterpret_cast<uintptr_t>(setPtrUC) & pointerMask) | ptrUnitigColors;
        }
        else if (flag == ptrBitmap){

            const uint32_t expected_sz = static_cast<uint32_t>(setBits >> 32);

            char* serialized = new char[expected_sz];
            stream_in.read(serialized, expected_sz);

            Bitmap* setPtrBmp = new Bitmap;

            setPtrBmp->r = std::move(Roaring::read(serialized));

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;

            delete[] serialized;
        }
        else if (flag == localTinyBitmap){

            TinyBitmap t_bmp;

            t_bmp.read(stream_in);

            setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        }

        return true;
    }

    return false;
}

UnitigColors::const_iterator UnitigColors::begin(const UnitigMapBase& um) const {

    const_iterator it(this, um.size - Kmer::k + 1, true);
    return ++it;
}

UnitigColors::const_iterator UnitigColors::begin(const size_t len_km_sz) const {

    const_iterator it(this, len_km_sz, true);
    return ++it;
}

UnitigColors::const_iterator UnitigColors::end() const {

    return const_iterator(this, 0, false);
}

UnitigColors UnitigColors::reverse(const UnitigMapBase& um) const {

    UnitigColors new_cs;

    if ((setBits & flagMask) == ptrUnitigColors){

        const UnitigColors* uc = getConstPtrUnitigColors();

        UnitigColors* setPtrUC = new UnitigColors[2];

        setPtrUC[0] = uc[0];
        setPtrUC[1] = uc[1].reverse(um);

        new_cs.setBits = (reinterpret_cast<uintptr_t>(setPtrUC) & pointerMask) | ptrUnitigColors;
    }
    else {

        const size_t len_unitig_km = um.size - Kmer::k + 1;

        UnitigColors::const_iterator it = begin(um), it_end = end();

        for (; it != it_end; ++it){

            const size_t new_km_dist = len_unitig_km - it.getKmerPosition() - 1;

            new_cs.add(len_unitig_km * it.getColorID() + new_km_dist);
        }
    }

    return new_cs;
}

UnitigColors::UnitigColors_const_iterator::UnitigColors_const_iterator() :  cs(nullptr), it_uc(nullptr), flag(localBitVector), it_setBits(0),
                                                                            ck_id(0xffffffffffffffff, 0xffffffffffffffff),
                                                                            cs_sz(0), um_sz(0), it_roar(empty_roar.end()) {}

UnitigColors::UnitigColors_const_iterator::UnitigColors_const_iterator(const UnitigColors* cs_, const size_t len_unitig_km, const bool beg) :
                                                                        cs(cs_), it_uc(nullptr), ck_id(0xffffffffffffffff, 0xffffffffffffffff),
                                                                        it_setBits(0xffffffffffffffff), um_sz(len_unitig_km), it_roar(empty_roar.end()) {

    flag = cs->setBits & flagMask;

    if (flag == ptrUnitigColors){

        const UnitigColors* uc = cs->getConstPtrUnitigColors();

        it_uc = new UnitigColors_const_iterator[2];

        it_uc[0] = beg ? uc[0].begin(1) : uc[0].end();
        it_uc[1] = beg ? uc[1].begin(len_unitig_km) : uc[1].end();

        cs_sz = cs->size();
    }
    else if (flag == ptrBitmap){

        it_roar = beg ? cs->getConstPtrBitmap()->r.begin() : cs->getConstPtrBitmap()->r.end();
        cs_sz = cs->size();
    }
    else if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = cs->getPtrTinyBitmap();

        t_bmp = &setPtrTinyBmp;
        it_t_bmp = beg ? t_bmp.begin() : t_bmp.end();
        cs_sz = cs->size();
    }
    else cs_sz = (flag == localSingleInt) ? 1 : maxBitVectorIDs;

    if (!beg) it_setBits = cs_sz;
}

UnitigColors::UnitigColors_const_iterator::UnitigColors_const_iterator(const UnitigColors_const_iterator& o) : cs(o.cs), flag(o.flag),
                                                                        it_setBits(o.it_setBits), cs_sz(o.cs_sz), um_sz(o.um_sz), ck_id(o.ck_id),
                                                                        it_uc(nullptr), it_roar(o.it_roar), it_t_bmp(o.it_t_bmp) {

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = cs->getPtrTinyBitmap();
        t_bmp = &setPtrTinyBmp;
    }
    else if (flag == ptrUnitigColors){

        it_uc = new UnitigColors_const_iterator[2];

        it_uc[0] = o.it_uc[0];
        it_uc[1] = o.it_uc[1];
    }
}

UnitigColors::UnitigColors_const_iterator::~UnitigColors_const_iterator() {

    t_bmp.detach();

    if (it_uc != nullptr) delete[] it_uc;
}

UnitigColors::UnitigColors_const_iterator& UnitigColors::UnitigColors_const_iterator::operator=(const UnitigColors_const_iterator& o) {

    cs = o.cs;

    flag = o.flag;

    it_setBits = o.it_setBits;
    cs_sz = o.cs_sz;
    um_sz = o.um_sz;

    ck_id = o.ck_id;

    it_roar = o.it_roar;
    it_t_bmp = o.it_t_bmp;

    it_uc = nullptr;

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = cs->getPtrTinyBitmap();
        t_bmp = &setPtrTinyBmp;
    }
    else if (flag == ptrUnitigColors){

        it_uc = new UnitigColors_const_iterator[2];

        it_uc[0] = o.it_uc[0];
        it_uc[1] = o.it_uc[1];
    }

    return *this;
}

pair<size_t, size_t> UnitigColors::UnitigColors_const_iterator::operator*() const { return ck_id; }

size_t UnitigColors::UnitigColors_const_iterator::getKmerPosition() const { return ck_id.first; }

size_t UnitigColors::UnitigColors_const_iterator::getColorID() const { return ck_id.second; }

UnitigColors::UnitigColors_const_iterator UnitigColors::UnitigColors_const_iterator::operator++(int) {

    UnitigColors_const_iterator tmp(*this);
    operator++();
    return tmp;
}

UnitigColors::UnitigColors_const_iterator& UnitigColors::UnitigColors_const_iterator::operator++() {

    if (it_setBits == cs_sz) return *this;

    ++it_setBits;

    if (flag == ptrUnitigColors) {

        if (it_setBits != 0){ // Must increment one iterator

            if (!it_uc[0].isInvalid() && (ck_id.second == it_uc[0].ck_id.second)){ // We were iterating over a full color

                ++ck_id.first;

                if (ck_id.first == um_sz) ++it_uc[0]; // If the next k-mer position is past the last
                else {

                    --it_setBits;
                    return *this;
                }
            }
            else if (!it_uc[1].isInvalid()) ++it_uc[1];
        }

        if (it_uc[0].isInvalid()) ck_id = it_uc[1].ck_id; // If first iterator finished, consider second iterator
        else if (it_uc[1].isInvalid() || (it_uc[0].ck_id.second < it_uc[1].ck_id.second)){

            ck_id.second = it_uc[0].ck_id.second;
            ck_id.first = 0;
        }
        else ck_id = it_uc[1].ck_id;
    }
    else if (flag == ptrBitmap) {

        if (it_setBits != 0) ++it_roar;
        if (it_roar != cs->getConstPtrBitmap()->r.end()){

            ck_id.first = *it_roar % um_sz;
            ck_id.second = *it_roar / um_sz;
        }
    }
    else if (flag == localTinyBitmap) {

        if (it_setBits != 0) ++it_t_bmp;
        if (it_t_bmp != t_bmp.end()){

            ck_id.first = *it_t_bmp % um_sz;
            ck_id.second = *it_t_bmp / um_sz;
        }
    }
    else if (flag == localBitVector){

        while (it_setBits < maxBitVectorIDs){

            if (((cs->setBits >> (it_setBits + shiftMaskBits)) & 0x1) != 0){

                ck_id.first = it_setBits % um_sz;
                ck_id.second = it_setBits / um_sz;

                break;
            }

            ++it_setBits;
        }
    }
    else {

        ck_id.first = cs->setBits >> shiftMaskBits;

        ck_id.second = ck_id.first / um_sz;
        ck_id.first %= um_sz;
    }

    return *this;
}

UnitigColors::UnitigColors_const_iterator& UnitigColors::UnitigColors_const_iterator::nextColor() {

    if (it_setBits == cs_sz) return *this;

    const size_t nextColor = ck_id.second + 1;
    const size_t nextPos = nextColor * um_sz;

    if (flag == ptrUnitigColors) {

        while (!it_uc[0].isInvalid() && (it_uc[0].ck_id.second < nextColor)) ++it_uc[0];
        while (!it_uc[1].isInvalid() && (it_uc[1].ck_id.second < nextColor)) it_uc[1].nextColor();

        if (it_uc[0].isInvalid()) ck_id = it_uc[1].ck_id;
        else if (it_uc[1].isInvalid() || (it_uc[0].ck_id.second < it_uc[1].ck_id.second)){

            ck_id.second = it_uc[0].ck_id.second;
            ck_id.first = ck_id.second * um_sz;
        }
        else ck_id = it_uc[1].ck_id;
    }
    else if (flag == ptrBitmap) {

        it_roar.equalorlarger(nextPos);

        if (it_roar != cs->getConstPtrBitmap()->r.end()){

            ck_id.first = *it_roar % um_sz;
            ck_id.second = *it_roar / um_sz;
        }
        else it_setBits = cs_sz;
    }
    else if (flag == localTinyBitmap) {

        while ((it_t_bmp != t_bmp.end()) && (*it_t_bmp < nextPos)) ++it_t_bmp;

        if (it_t_bmp != t_bmp.end()){

            ck_id.first = *it_t_bmp % um_sz;
            ck_id.second = *it_t_bmp / um_sz;
        }
        else it_setBits = cs_sz;
    }
    else if (flag == localBitVector){

        it_setBits = min(maxBitVectorIDs, nextPos);

        while (it_setBits < maxBitVectorIDs){

            if ((cs->setBits >> (it_setBits + shiftMaskBits)) & 0x1){

                ck_id.first = it_setBits % um_sz;
                ck_id.second = it_setBits / um_sz;

                break;
            }

            ++it_setBits;
        }
    }
    else ++it_setBits;

    return *this;
}

bool UnitigColors::UnitigColors_const_iterator::operator==(const UnitigColors_const_iterator& o) const {

    if ((cs == o.cs) && (flag == o.flag) && (cs_sz == o.cs_sz)){

        if (flag == ptrUnitigColors) return (it_uc[0] == o.it_uc[0]) && (it_uc[1] == o.it_uc[1]);
        if (flag == ptrBitmap) return (it_roar == o.it_roar);
        if (flag == localTinyBitmap) return (it_t_bmp == o.it_t_bmp);

        return (it_setBits == o.it_setBits);
    }

    return false;
}

bool UnitigColors::UnitigColors_const_iterator::operator!=(const UnitigColors_const_iterator& o) const {

    return !operator==(o);
}

const size_t UnitigColors::maxBitVectorIDs = 61; // 64 bits - 3 bits for the color set type

const uintptr_t UnitigColors::localTinyBitmap = 0x0;
const uintptr_t UnitigColors::localBitVector = 0x1;
const uintptr_t UnitigColors::localSingleInt = 0x2;
const uintptr_t UnitigColors::ptrBitmap = 0x3;
const uintptr_t UnitigColors::ptrUnitigColors = 0x4;

const size_t UnitigColors::shiftMaskBits = 3;

const uintptr_t UnitigColors::flagMask = 0x7;
const uintptr_t UnitigColors::pointerMask = 0xFFFFFFFFFFFFFFF8;
