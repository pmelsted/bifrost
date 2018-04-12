#include "TinyBitmap.hpp"

TinyBitmap::TinyBitmap() : tiny_bmp(nullptr) {}

TinyBitmap::TinyBitmap(const TinyBitmap& o) : tiny_bmp(nullptr) {

    if (o.tiny_bmp != nullptr){

        const uint16_t sz = sizes[o.getIndexSize()];

        tiny_bmp = new uint16_t[sz];

        std::copy(o.tiny_bmp, o.tiny_bmp + sz, tiny_bmp);
    }
}

TinyBitmap& TinyBitmap::operator=(const TinyBitmap& o) {

    if (this != &o){

        empty();

        if (o.tiny_bmp != nullptr){

            const uint16_t sz = sizes[o.getIndexSize()];

            tiny_bmp = new uint16_t[sz];

            std::copy(o.tiny_bmp, o.tiny_bmp + sz, tiny_bmp);
        }
    }

    return *this;
}

TinyBitmap::~TinyBitmap() { empty(); }

void TinyBitmap::empty() {

    if (tiny_bmp != nullptr){

        delete[] tiny_bmp;
        tiny_bmp = nullptr;
    }
}

bool TinyBitmap::add(const uint32_t val){

    //cout << "TinyBitmap::add(" << val << ")" << endl;

    const uint16_t val_div = val >> 16;
    const uint16_t val_mod = val & 0xFFFF;

    if (tiny_bmp == nullptr){

        tiny_bmp = new uint16_t[sizes[0]](); // Size = 0, offset = 0; bmp_mode
        tiny_bmp[2] = val_div;
    }

    if (getOffset() != val_div) return false;

    uint16_t sz_idx = getIndexSize();
    uint16_t sz = sizes[sz_idx];
    uint16_t mode = getMode();
    uint16_t cardinality = getCardinality();

    /*cout << "TinyBitmap::add(): mode = " << mode << endl;
    cout << "TinyBitmap::add(): cardinality = " << cardinality << endl;*/

    // Compute if inserting new value triggers an increase of the container size
    if (((mode == bmp_mode) && (val_mod >= ((sz - 3) << 4))) || ((mode != bmp_mode) && (cardinality == (sz - 3 - (mode == rle_list_mode))))){

        // Means that in its current mode, container size must be increased to add a value
        // We need to compute if which mode has the smaller container size

        bool inc = true;

        if ((mode != bmp_mode) && contains(val)) return true;

        if (runOptimize() != 0){

            sz_idx = getIndexSize();
            sz = sizes[sz_idx];
            mode = getMode();
            cardinality = getCardinality();

            if (cardinality != (sz - 4)) inc = false;
        }

        if (inc){

            const uint16_t max_val_mode = std::max(val_mod, static_cast<const uint16_t>(maximum() & 0xFFFF)) + 1;
            const uint16_t nb_uint_bmp = sizes[getIndexLargerSize((max_val_mode >> 4) + ((max_val_mode & 0xF) != 0) + 3)];

            bool res;

                /*cout << "TinyBitmap::add(): max_val_mode = " << max_val_mode << endl;
                cout << "TinyBitmap::add(): nb_uint_bmp = " << nb_uint_bmp << endl;*/

            if (mode == rle_list_mode){

                const uint16_t nb_val = size();
                const uint16_t nb_val_rle_list = sizes[sz_idx + 1];
                const uint16_t nb_val_list = sizes[getIndexLargerSize(nb_val + 4)];
                const uint16_t nb_val_min = (nb_val > (0xFFFF - 48)) ? 0xFFFF : std::min(nb_val_rle_list, std::min(nb_val_list, nb_uint_bmp));

                /*cout << "TinyBitmap::add(): Size must be increased" << endl;
                cout << "TinyBitmap::add(): nb_val_rle_list = " << nb_val_rle_list << endl;
                cout << "TinyBitmap::add(): nb_val_list = " << nb_val_list << endl;
                cout << "TinyBitmap::add(): nb_uint_bmp = " << nb_uint_bmp << endl;
                cout << "TinyBitmap::add(): nb_val_min = " << nb_val_min << endl;*/

                if (nb_val_min > sizes[nb_sizes - 1]) return false;

                res = (nb_val_rle_list == nb_val_min);
                res = res ? change_sz(nb_val_min) : switch_mode(nb_val_min, (nb_val_list <= nb_uint_bmp) ? list_mode : bmp_mode);

                cardinality = getCardinality();
            }
            else {

                const uint16_t nb_val_list = sizes[getIndexLargerSize(cardinality + 4)];

                /*cout << "TinyBitmap::add(): Size must be increased" << endl;
                cout << "TinyBitmap::add(): max_val_mode = " << max_val_mode << endl;
                cout << "TinyBitmap::add(): nb_uint_bmp = " << nb_uint_bmp << endl;
                cout << "TinyBitmap::add(): nb_val_list = " << nb_val_list << endl;*/

                if (mode == bmp_mode) res = (nb_uint_bmp <= nb_val_list) ? change_sz(nb_uint_bmp) : switch_mode(nb_val_list, list_mode);
                else res = (nb_val_list <= nb_uint_bmp) ? change_sz(nb_val_list) : switch_mode(nb_uint_bmp, bmp_mode);
            }

            if (!res) return false;

            mode = getMode();
        }
    }

    //cout << "TinyBitmap::add(): mode = " << mode << endl;

    if (mode == bmp_mode){ // Bitmap mode

        uint16_t& div = tiny_bmp[(val_mod >> 4) + 3];

        const uint16_t mod = 1U << (val_mod & 0xF);

        tiny_bmp[1] += ((div & mod) == 0); // += (1 << 8) if not already set (increase cardinality in header), 0 otherwise
        div |= mod; // Insert new value
    }
    else if (mode == list_mode) { // Binary search

        if (cardinality == 0){

            ++tiny_bmp[1];
            tiny_bmp[3] = val_mod;
        }
        else {

            uint16_t imid;

            uint16_t imin = 3;
            uint16_t imax = cardinality + 2;

            while (imin < imax){

                imid = (imin + imax) >> 1;

                if (tiny_bmp[imid] < val_mod) imin = imid + 1;
                else imax = imid;
            }

            if (tiny_bmp[imin] != val_mod){

                if (tiny_bmp[imin] < val_mod){

                    std::memmove(&tiny_bmp[imin + 2], &tiny_bmp[imin + 1], (cardinality + 2 - imin) * sizeof(uint16_t)); // Shift values
                    tiny_bmp[imin + 1] = val_mod; // Insert value
                }
                else {

                    std::memmove(&tiny_bmp[imin + 1], &tiny_bmp[imin], (cardinality + 3 - imin) * sizeof(uint16_t)); // Shift values
                    tiny_bmp[imin] = val_mod; // Insert value
                }

                ++tiny_bmp[1]; // Increase cardinality
            }
        }
    }
    else { // Binary search

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 1;

        while (imin < imax){

            imid = (imin + imax) >> 1;
            imid -= ((imid & 0x1) == 0);

            if (tiny_bmp[imid + 1] < val_mod) imin = imid + 2;
            else imax = imid;
        }

        if ((tiny_bmp[imin] < val_mod) || (val_mod > tiny_bmp[imin + 1])){

            if (val_mod > tiny_bmp[imin + 1]){

                if (val_mod == (tiny_bmp[imin + 1] + 1)) ++tiny_bmp[imin + 1];
                else {

                    std::memmove(&tiny_bmp[imin + 4], &tiny_bmp[imin + 2], (cardinality + 1 - imin) * sizeof(uint16_t)); // Shift values

                    tiny_bmp[imin + 2] = val_mod; // Insert start run
                    tiny_bmp[imin + 3] = val_mod; // Insert end run
                    tiny_bmp[1] += 2; // Increase cardinality
                }
            }
            else if (val_mod == (tiny_bmp[imin] - 1)) --tiny_bmp[imin];
            else {

                std::memmove(&tiny_bmp[imin + 2], &tiny_bmp[imin], (cardinality + 3 - imin) * sizeof(uint16_t)); // Shift values

                tiny_bmp[imin] = val_mod; // Insert start run
                tiny_bmp[imin + 1] = val_mod; // Insert end run
                tiny_bmp[1] += 2; // Increase cardinality
            }
        }
    }

    return true;
}

bool TinyBitmap::contains(const uint32_t val) const {

    // If not allocated or cardinality is 0, val is not present
    if ((tiny_bmp == nullptr) || (getCardinality() == 0) || ((val >> 16) != getOffset())) return false;

    const uint16_t mode = getMode();
    const uint16_t cardinality = getCardinality();

    const uint16_t val_mod = val & 0xFFFF;

    // Bitmap mode
    if (mode == bmp_mode){

        if (val_mod >= ((sizes[getIndexSize()] - 3) << 4)) return false;

        return ((tiny_bmp[(val_mod >> 4) + 3] & (1U << (val_mod & 0xF))) != 0);
    }
    else if (mode == list_mode) {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 2;

        while (imin < imax){

            imid = (imin + imax) >> 1;

            if (tiny_bmp[imid] < val_mod) imin = imid + 1;
            else imax = imid;
        }

        return (tiny_bmp[imin] == val_mod);
    }
    else {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 1;

        while (imin < imax){

            imid = (imin + imax) >> 1;
            imid -= ((imid & 0x1) == 0);

            if (tiny_bmp[imid + 1] < val_mod) imin = imid + 2;
            else imax = imid;
        }

        return ((val_mod >= tiny_bmp[imin]) && (val_mod <= tiny_bmp[imin + 1]));
    }

    return false;
}

void TinyBitmap::remove(const uint32_t val){

    // If not allocated or cardinality is 0, val is not present
    if ((tiny_bmp == nullptr) || (getCardinality() == 0) || ((val >> 16) != getOffset())) return;

    const uint16_t mode = getMode();
    const uint16_t cardinality = getCardinality();
    const uint16_t val_mod = val & 0xFFFF;

    bool try_decrease_sz = false;

    // Bitmap mode
    if (mode == bmp_mode){

        if (val_mod >= ((sizes[getIndexSize()] - 3) << 4)) return;

        tiny_bmp[(val_mod >> 4) + 3] &= ~(1U << (val_mod & 0xF));
        --tiny_bmp[1];
        try_decrease_sz = true;
    }
    else if (mode == list_mode) {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 2;

        while (imin < imax){

            imid = (imin + imax) >> 1;

            if (tiny_bmp[imid] < val_mod) imin = imid + 1;
            else imax = imid;
        }

        if (tiny_bmp[imin] == val_mod){

            std::memmove(&tiny_bmp[imin], &tiny_bmp[imin + 1], (cardinality + 2 - imin) * sizeof(uint16_t)); // Shift values
            --tiny_bmp[1];
            try_decrease_sz = true;
        }
    }
    else {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 1;

        while (imin < imax){

            imid = (imin + imax) >> 1;
            imid -= ((imid & 0x1) == 0);

            if (tiny_bmp[imid + 1] < val_mod) imin = imid + 2;
            else imax = imid;
        }

        if ((val_mod >= tiny_bmp[imin]) && (val_mod <= tiny_bmp[imin + 1])){

            if ((val_mod == tiny_bmp[imin]) && (val_mod == tiny_bmp[imin + 1])){ // The run is the value to delete

                std::memmove(&tiny_bmp[imin], &tiny_bmp[imin + 2], (cardinality + 3 - imin) * sizeof(uint16_t)); // Shift values
                tiny_bmp[1] -= 2;
            }
            else if (val_mod == tiny_bmp[imin]) ++tiny_bmp[imin];
            else if (val_mod == tiny_bmp[imin + 1]) --tiny_bmp[imin + 1];
            else {

                uint16_t sz = sizes[getIndexSize()];

                if (cardinality != (sz - 4)){ // There is enough space to insert a new run

                    std::memmove(&tiny_bmp[imin + 3], &tiny_bmp[imin + 1], (cardinality + 1 - imin) * sizeof(uint16_t)); // Shift values

                    tiny_bmp[imin + 1] = val_mod - 1; // Insert start run
                    tiny_bmp[imin + 2] = val_mod + 1; // Insert start run
                    tiny_bmp[1] += 2; // Increase cardinality
                }
                else {

                    switch_mode(size() + 3, list_mode);
                    remove(val);
                }
            }
        }
    }

    if (try_decrease_sz){ // Can only be bmp_mode || list_mode

        const uint16_t sz = sizes[getIndexSize()];
        const uint16_t max_val_mod = static_cast<const uint16_t>(maximum() & 0xFFFF) + 1;
        const uint16_t nb_uint_bmp = sizes[getIndexLargerSize((max_val_mod >> 4) + ((max_val_mod & 0xF) != 0) + 3)];
        const uint16_t nb_val_list = sizes[getIndexLargerSize(cardinality + 2)];

        if (std::min(nb_uint_bmp, nb_val_list) < sz){

            if (mode == bmp_mode) (nb_uint_bmp <= nb_val_list) ? change_sz(nb_uint_bmp) : switch_mode(nb_val_list, list_mode);
            else (nb_val_list <= nb_uint_bmp) ? change_sz(nb_val_list) : switch_mode(nb_uint_bmp, bmp_mode);
        }
    }
}

uint32_t TinyBitmap::maximum() const {

    if ((tiny_bmp == nullptr) || (getCardinality() == 0)) return 0; // If not allocated or cardinality is 0, return 0

    const uint32_t offset = getOffset() << 16;

    if (getMode() == bmp_mode){

        for (uint16_t i = sizes[getIndexSize()] - 1; i != 2; --i){

            for (uint16_t j = 15, e = tiny_bmp[i]; e != 0; --j, e <<= 1){

                if ((e & 0x8000) != 0) return offset | (((i - 3) << 4) + j);
            }
        }
    }

    return offset | tiny_bmp[getCardinality() + 2]; // mode = list_mode
}

size_t TinyBitmap::getSizeInBytes() const {

    if (tiny_bmp == nullptr) return sizeof(uint16_t*);
    return sizeof(uint16_t*) + sizes[getIndexSize()] * sizeof(uint16_t);
}

size_t TinyBitmap::size() const {

    if (tiny_bmp == nullptr) return 0;

    const uint16_t cardinality = getCardinality();

    if (getMode() == rle_list_mode){

        size_t card = 0;

        for (size_t i = 3; i < cardinality + 3; i += 2) card += tiny_bmp[i + 1] - tiny_bmp[i];

        return card + (cardinality >> 1);
    }

    return cardinality;
}

size_t TinyBitmap::runOptimize() {

    if (tiny_bmp != nullptr){

        const uint16_t mode = getMode();
        const uint16_t cardinality = getCardinality();

        if ((mode != rle_list_mode) && (cardinality != 0)){

            const uint16_t sz = sizes[getIndexSize()];

            if (mode == bmp_mode){

                uint16_t nb_run = 0;
                uint16_t prev_val_pres = 0xFFFF; // This value cannot exist in bitmap mode
                uint16_t cardinality_cpy = cardinality;

                for (uint16_t i = 3; (i != sz) && (cardinality_cpy != 0); ++i){

                    for (uint16_t j = (i - 3) << 4, e = tiny_bmp[i]; e != 0; e >>= 1, ++j){

                        if (e & 0x1){

                            nb_run += (j != (prev_val_pres + 1));
                            --cardinality_cpy;
                            prev_val_pres = j;
                        }
                    }
                }

                const uint16_t new_cardinality = nb_run << 1;
                const uint16_t new_sz_idx = getIndexLargerSize(new_cardinality + 3);
                const uint16_t new_sz = sizes[new_sz_idx];

                //cout << "TinyBitmap::runOptimize(): new_cardinality = " << new_cardinality << endl;

                if ((new_sz < sz) && (new_sz <= sizes[nb_sizes - 1])){

                    uint16_t k = 2;

                    cardinality_cpy = cardinality;
                    prev_val_pres = 0xFFFF;

                    uint16_t* tiny_bmp_new = new uint16_t[new_sz]();

                    for (uint16_t i = 3; (i != sz) && (cardinality_cpy != 0); ++i){

                        for (uint16_t j = (i - 3) << 4, e = tiny_bmp[i]; e != 0; e >>= 1, ++j){

                            if (e & 0x1){

                                if (j != (prev_val_pres + 1)){

                                    tiny_bmp_new[k++] = prev_val_pres;
                                    tiny_bmp_new[k++] = j;
                                }

                                --cardinality_cpy;
                                prev_val_pres = j;
                            }
                        }
                    }

                    tiny_bmp_new[k] = prev_val_pres;
                    tiny_bmp_new[0] = (new_sz_idx << 1) | rle_list_mode | bits_16;
                    tiny_bmp_new[1] = new_cardinality;
                    tiny_bmp_new[2] = getOffset();

                    delete[] tiny_bmp;
                    tiny_bmp = tiny_bmp_new;

                    return sz - new_sz;
                }
            }
            else {

                uint16_t nb_run = 1;

                for (size_t i = 4; i < cardinality + 3; ++i) nb_run += (tiny_bmp[i] != (tiny_bmp[i-1] + 1));

                const uint16_t new_cardinality = nb_run << 1;
                const uint16_t new_sz_idx = getIndexLargerSize(new_cardinality + 3);
                const uint16_t new_sz = sizes[new_sz_idx];

                if ((new_sz < sz) && (new_sz <= sizes[nb_sizes - 1])){

                    uint16_t k = 4;

                    uint16_t* tiny_bmp_new = new uint16_t[new_sz]();

                    tiny_bmp_new[0] = (new_sz_idx << 1) | rle_list_mode | bits_16;
                    tiny_bmp_new[1] = new_cardinality;
                    tiny_bmp_new[2] = getOffset();
                    tiny_bmp_new[3] = tiny_bmp[3];

                    for (size_t i = 4; i < cardinality + 3; ++i){

                        if (tiny_bmp[i] != (tiny_bmp[i-1] + 1)){

                            tiny_bmp_new[k++] = tiny_bmp[i-1];
                            tiny_bmp_new[k++] = tiny_bmp[i];
                        }
                    }

                    tiny_bmp_new[k] = tiny_bmp[cardinality + 2];

                    delete[] tiny_bmp;
                    tiny_bmp = tiny_bmp_new;

                    return sz - new_sz;
                }
            }
        }
    }

    return 0;
}

void TinyBitmap::print() const {

    if (tiny_bmp == nullptr) return;

    const uint16_t mode = getMode();
    const uint16_t sz = sizes[getIndexSize()];
    const uint16_t cardinality = getCardinality();

    if (mode == list_mode){

        for (size_t i = 3; i < cardinality + 3; ++i) cout << tiny_bmp[i] << endl;
    }
    else if (mode == rle_list_mode){

        for (size_t i = 3; i < cardinality + 3; i += 2){

            for (uint16_t j = tiny_bmp[i]; j <= tiny_bmp[i + 1]; ++j) cout << j << endl;
        }
    }
}

bool TinyBitmap::change_sz(const uint16_t sz_min) {

    if (sz_min > sizes[nb_sizes - 1]) return false;

    const bool is_allocated = (tiny_bmp != nullptr);

    const uint16_t sz = is_allocated ? sizes[getIndexSize()] : 0;

    const uint16_t new_sz_idx = getIndexLargerSize(sz_min);
    const uint16_t new_sz = sizes[new_sz_idx];

    if (is_allocated){

        uint16_t* tiny_bmp_new = new uint16_t[new_sz]();

        std::copy(tiny_bmp, tiny_bmp + (new_sz >= sz ? sz : sz_min), tiny_bmp_new);
        delete[] tiny_bmp;

        tiny_bmp = tiny_bmp_new;
        tiny_bmp[0] = (tiny_bmp[0] & ~sz_mask) | (new_sz_idx << 1);
    }
    else {

        tiny_bmp = new uint16_t[new_sz]();
        tiny_bmp[0] = bmp_mode | (new_sz_idx << 1);
    }

    return true;
}

bool TinyBitmap::switch_mode(const uint16_t sz_min, const uint16_t new_mode) {

    //cout << "TinyBitmap::switch_mode(" << sz_min << ", " << new_mode << ")" << endl;

    if (tiny_bmp == nullptr) return true;

    const uint16_t mode = getMode();
    const uint16_t sz = sizes[getIndexSize()];
    const uint16_t cardinality = getCardinality();
    const uint16_t offset = getOffset();

    uint16_t* tiny_bmp_new = nullptr;

    if ((mode == bmp_mode) && (new_mode == list_mode)) { // Switch to list mode

        const uint16_t new_sz_min = std::max(sz_min, static_cast<const uint16_t>(cardinality + 3));

        if (new_sz_min > sizes[nb_sizes - 1]) return false;

        std::swap(tiny_bmp, tiny_bmp_new);
        change_sz(new_sz_min);

        for (uint16_t i_bmp = 3, i_list = 3, card = cardinality; (i_bmp < sz) && (card != 0); ++i_bmp){

            for (uint16_t j = (i_bmp - 3) << 4; tiny_bmp_new[i_bmp] != 0; tiny_bmp_new[i_bmp] >>= 1, ++j){

                if (tiny_bmp_new[i_bmp] & 0x1){

                    tiny_bmp[i_list++] = j;
                    --card;
                }
            }
        }

        tiny_bmp[0] = (tiny_bmp[0] & sz_mask) | list_mode | bits_16;
        tiny_bmp[1] = cardinality;
        tiny_bmp[2] = offset;
    }
    else if ((mode == list_mode) && (new_mode == bmp_mode)) { // mode is list

        const uint16_t max_val_mod = static_cast<const uint16_t>((maximum() & 0xFFFF) + 1);
        const uint16_t new_sz_min = std::max(sz_min, static_cast<const uint16_t>((max_val_mod >> 4) + ((max_val_mod & 0xF) != 0) + 3));

        if (new_sz_min > sizes[nb_sizes - 1]) return false;

        std::swap(tiny_bmp, tiny_bmp_new);
        change_sz(new_sz_min);

        for (size_t i = 3; i < cardinality + 3; ++i) tiny_bmp[(tiny_bmp_new[i] >> 4) + 3] |= (1U << (tiny_bmp_new[i] & 0xF));

        tiny_bmp[0] = (tiny_bmp[0] & sz_mask) | bmp_mode | bits_16;
        tiny_bmp[1] = cardinality;
        tiny_bmp[2] = offset;
    }
    else if (mode == rle_list_mode){

        if (new_mode == bmp_mode){

            uint16_t new_card = 0;

            const uint16_t max_val_mod = static_cast<const uint16_t>((maximum() & 0xFFFF) + 1);
            const uint16_t new_sz_min = std::max(sz_min, static_cast<const uint16_t>((max_val_mod >> 4) + ((max_val_mod & 0xF) != 0) + 3));

            if (new_sz_min > sizes[nb_sizes - 1]) return false;

            std::swap(tiny_bmp, tiny_bmp_new);
            change_sz(new_sz_min);

            for (size_t i = 3; i < cardinality + 3; i += 2){

                for (uint16_t j = tiny_bmp_new[i]; j <= tiny_bmp_new[i + 1]; ++j, ++new_card) tiny_bmp[(j >> 4) + 3] |= (1U << (j & 0xF));
            }

            tiny_bmp[0] = (tiny_bmp[0] & sz_mask) | bmp_mode | bits_16;
            tiny_bmp[1] = new_card;
            tiny_bmp[2] = offset;
        }
        else if (new_mode == list_mode){

            const uint16_t new_card = size();
            const uint16_t new_sz_min = std::max(sz_min, static_cast<const uint16_t>(new_card + 3));

            if (new_sz_min > sizes[nb_sizes - 1]) return false;

            std::swap(tiny_bmp, tiny_bmp_new);
            change_sz(new_sz_min);

            for (size_t i = 3, k = 3; i < cardinality + 3; i += 2){

                for (uint16_t j = tiny_bmp_new[i]; j <= tiny_bmp_new[i + 1]; ++j, ++k) tiny_bmp[k] = j;
            }

            tiny_bmp[0] = (tiny_bmp[0] & sz_mask) | list_mode | bits_16;
            tiny_bmp[1] = new_card;
            tiny_bmp[2] = offset;
        }
    }

    if (tiny_bmp_new != nullptr) delete[] tiny_bmp_new;

    return true;
}

const uint16_t TinyBitmap::sz_mask = 0x003E;    // 0000000000111110
const uint16_t TinyBitmap::mode_mask = 0x00C0;  // 0000000011000000
const uint16_t TinyBitmap::bits_mask = 0x0001;  // 0000000000000001

const uint16_t TinyBitmap::bmp_mode = 0x0000;
const uint16_t TinyBitmap::list_mode = 0x0040;
const uint16_t TinyBitmap::rle_list_mode = 0x0080;

const uint16_t TinyBitmap::bits_16 = 0x0000;
const uint16_t TinyBitmap::bits_32 = 0x0001;

const uint16_t TinyBitmap::sz_min = 8;
const uint16_t TinyBitmap::sz_max = 4096;

const uint16_t TinyBitmap::sizes[] = {  8, 16, 32, 64, 96, 144, 216, 324,
                                        486, 730, 1096, 1370, 1712, 2140, 2676, 3346,
                                        4096, 0xFFFF };

const uint16_t TinyBitmap::nb_sizes = 17;
