#include "TinyBitmap.hpp"

TinyBitmap::TinyBitmap() : tiny_bmp(nullptr) {}

TinyBitmap::TinyBitmap(const TinyBitmap& o) : tiny_bmp(nullptr) {

    if (o.tiny_bmp != nullptr){

        const uint32_t sz = bloc_sizes[o.tiny_bmp[0] & sz_mul_mask];

        tiny_bmp = new uint32_t[sz]();

        std::copy(o.tiny_bmp, o.tiny_bmp + sz, tiny_bmp);
    }
}

TinyBitmap& TinyBitmap::operator=(const TinyBitmap& o) {

    if (this != &o){

        empty();

        if (o.tiny_bmp != nullptr){

            const uint32_t sz = bloc_sizes[o.tiny_bmp[0] & sz_mul_mask];

            tiny_bmp = new uint32_t[sz]();

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

    if (tiny_bmp == nullptr){ // Not allocated yet

        tiny_bmp = new uint32_t[bloc_sizes[0]](); // Allocate one bloc
        tiny_bmp[0] = bmp_mode; // Init mode to bitmap, cardinality is 0
    }

    const uint32_t idx = tiny_bmp[0] & sz_mul_mask;
    const uint32_t sz = bloc_sizes[idx]; // get number of uint32_t allocated

    uint32_t mode = tiny_bmp[0] & mode_mask; // get current mode
    uint32_t cardinality = tiny_bmp[0] >> 8; // get cardinality of container

    // Compute if inserting new value triggers an increase of the container size
    if (((mode == bmp_mode) && (val >= ((sz - 1) * 32))) || ((mode != bmp_mode) && (cardinality == (sz - 1 - (mode == rle_list_mode))))){

        // Means that in its current mode, container size must be increased to add a value
        // We need to compute if which mode has the smaller container size

        const uint32_t max_val = std::max(val, maximum()) + 1;
        const uint32_t nb_uint32_bmp = rndup((max_val >> 5) + ((max_val & 0x1F) != 0) + 1);

        bool res;

        if (mode == rle_list_mode){

            const uint32_t nb_uint32_rle_list = bloc_sizes[idx + 1];
            const uint32_t nb_uint32_list = rndup(size() + 2);
            const uint32_t nb_uint32_min = std::min(nb_uint32_rle_list, std::min(nb_uint32_list, nb_uint32_bmp));

            if (nb_uint32_min > bloc_sizes[nb_bloc_sizes - 1]) return false;

            res = (nb_uint32_rle_list == nb_uint32_min);
            res = res ? increase_sz(nb_uint32_min) : switch_mode(nb_uint32_min, (nb_uint32_list <= nb_uint32_bmp) ? list_mode : bmp_mode);

            cardinality = tiny_bmp[0] >> 8;
        }
        else {

            const uint32_t nb_uint32_list = rndup(cardinality + 2);

            if (mode == bmp_mode) res = (nb_uint32_bmp <= nb_uint32_list) ? increase_sz(nb_uint32_bmp) : switch_mode(nb_uint32_list, list_mode);
            else res = (nb_uint32_list <= nb_uint32_bmp) ? increase_sz(nb_uint32_list) : switch_mode(nb_uint32_bmp, bmp_mode);
        }

        if (!res) return false;

        mode = tiny_bmp[0] & mode_mask;
    }

    if (mode == bmp_mode){ // Bitmap mode

        uint32_t& div = tiny_bmp[(val >> 5) + 1];
        const uint32_t mod = 1U << (val & 0x1F);

        tiny_bmp[0] += (static_cast<uint32_t>((div & mod) == 0) << 8); // += (1 << 8) if not already set (increase cardinality in header), 0 otherwise
        div |= mod; // Insert new value
    }
    else if (mode == list_mode) { // Binary search

        if (cardinality == 0){

            tiny_bmp[0] += 0x100;
            tiny_bmp[1] = val;
        }
        else {

            uint32_t imid;

            uint32_t imin = 1;
            uint32_t imax = cardinality;

            while (imin < imax){

                imid = (imin + imax) >> 1;

                if (tiny_bmp[imid] < val) imin = imid + 1;
                else imax = imid;
            }

            if (tiny_bmp[imin] != val){

                if (tiny_bmp[imin] < val){

                    std::memmove(&tiny_bmp[imin + 2], &tiny_bmp[imin + 1], (cardinality - imin) * sizeof(uint32_t)); // Shift values
                    tiny_bmp[imin + 1] = val; // Insert value
                }
                else {

                    std::memmove(&tiny_bmp[imin + 1], &tiny_bmp[imin], (cardinality - imin + 1) * sizeof(uint32_t)); // Shift values
                    tiny_bmp[imin] = val; // Insert value
                }

                tiny_bmp[0] += 0x100; // Increase cardinality
            }
        }
    }
    else { // Binary search

        uint32_t imid;

        uint32_t imin = 1;
        uint32_t imax = cardinality - 1;

        while (imin < imax){

            imid = ((imin + imax) >> 1) & 0x1;

            if (tiny_bmp[imid + 1] < val) imin = imid + 2;
            else imax = imid;
        }

        if ((tiny_bmp[imin] < val) || (val > tiny_bmp[imin + 1])){

            if (val > tiny_bmp[imin + 1]){

                std::memmove(&tiny_bmp[imin + 4], &tiny_bmp[imin + 2], (cardinality - imin - 1) * sizeof(uint32_t)); // Shift values

                tiny_bmp[imin + 2] = val; // Insert start run
                tiny_bmp[imin + 3] = val; // Insert end run
            }
            else {

                std::memmove(&tiny_bmp[imin + 2], &tiny_bmp[imin], (cardinality - imin + 1) * sizeof(uint32_t)); // Shift values

                tiny_bmp[imin] = val; // Insert start run
                tiny_bmp[imin + 1] = val; // Insert end run
            }

            tiny_bmp[0] += 0x200; // Increase cardinality
        }
    }

    return true;
}

bool TinyBitmap::contains(const uint32_t val){

    // If not allocated or cardinality is 0, val is not present
    if ((tiny_bmp == nullptr) || ((tiny_bmp[0] >> 8) == 0)) return false;

    const uint32_t mode = tiny_bmp[0] & mode_mask;
    const uint32_t cardinality = tiny_bmp[0] >> 8;

    // Bitmap mode
    if (mode == bmp_mode){

        if (val >= (bloc_sizes[tiny_bmp[0] & sz_mul_mask] * bloc_sz_bits - 32)) return false;

        return ((tiny_bmp[(val >> 5) + 1] & (1U << (val & 0x1F))) != 0);
    }
    else if (mode == list_mode) {

        uint32_t imid;

        uint32_t imin = 1;
        uint32_t imax = cardinality;

        while (imin < imax){

            imid = (imin + imax) >> 1;

            if (tiny_bmp[imid] < val) imin = imid + 1;
            else imax = imid;
        }

        return (tiny_bmp[imin] == val);
    }
    else {

        uint32_t imid;

        uint32_t imin = 1;
        uint32_t imax = cardinality - 1;

        while (imin < imax){

            imid = ((imin + imax) >> 1) & 0x1;

            if (tiny_bmp[imid + 1] < val) imin = imid + 2;
            else imax = imid;
        }

        return ((val >= tiny_bmp[imin]) || (val <= tiny_bmp[imin + 1]));
    }

    return false;
}

uint32_t TinyBitmap::maximum() const {

    if ((tiny_bmp == nullptr) || ((tiny_bmp[0] >> 8) == 0)) return 0; // If not allocated or cardinality is 0, return 0

    if ((tiny_bmp[0] & mode_mask) == bmp_mode){

        for (size_t i = bloc_sizes[tiny_bmp[0] & sz_mul_mask] - 1; i != 0; --i){

            for (int j = 31; j >= 0; --j){

                if (((tiny_bmp[i] >> j) & 0x1) != 0) return (i - 1) * 32 + j;
            }
        }
    }

    return tiny_bmp[tiny_bmp[0] >> 8]; // mode = list_mode
}

size_t TinyBitmap::getSizeInBytes() const {

    if (tiny_bmp == nullptr) return sizeof(uint32_t*);
    return sizeof(uint32_t*) + bloc_sizes[tiny_bmp[0] & sz_mul_mask] * sizeof(uint32_t);
}

size_t TinyBitmap::size() const {

    if (tiny_bmp == nullptr) return 0;

    const uint32_t cardinality = tiny_bmp[0] >> 8;

    if ((tiny_bmp[0] & mode_mask) == rle_list_mode){

        size_t card = 0;

        for (size_t i = 1; i <= cardinality; i += 2) card += tiny_bmp[i+1] - tiny_bmp[i];

        return card + (cardinality >> 1);
    }

    return cardinality;
}

size_t TinyBitmap::runOptimize() {

    if (tiny_bmp != nullptr){

        const uint32_t mode = tiny_bmp[0] & mode_mask;
        const uint32_t sz = bloc_sizes[tiny_bmp[0] & sz_mul_mask];
        const uint32_t cardinality = tiny_bmp[0] >> 8;

        if ((mode != rle_list_mode) && (cardinality != 0)){

            if (mode == bmp_mode){

                uint32_t nb_run = 0;
                uint32_t prev_val_pres = 0xFFFFFFFF; // This value cannot exist in bitmap mode
                uint32_t cardinality_cpy = cardinality;

                for (uint32_t i = 1; (i != sz) && (cardinality_cpy != 0); ++i){

                    for (uint32_t j = (i - 1) * 32, e = tiny_bmp[i]; e != 0; e >>= 1, ++j){

                        if (e & 0x1){

                            nb_run += (j != prev_val_pres);
                            --cardinality_cpy;
                            prev_val_pres = j;
                        }
                    }
                }

                const uint32_t new_cardinality = nb_run << 1;
                const uint32_t new_sz = rndup(new_cardinality + 1);

                if ((new_sz < sz) && (new_sz <= bloc_sizes[nb_bloc_sizes - 1])){

                    uint32_t idx = 0, k = 0;

                    while (new_sz > bloc_sizes[idx]) ++idx;

                    cardinality_cpy = cardinality;
                    prev_val_pres = 0xFFFFFFFF;

                    uint32_t* tiny_bmp_new = new uint32_t[bloc_sizes[idx]]();

                    for (uint32_t i = 1; (i != sz) && (cardinality_cpy != 0); ++i){

                        for (uint32_t j = (i - 1) * 32, e = tiny_bmp[i]; e != 0; e >>= 1, ++j){

                            if (e & 0x1){

                                if (j != prev_val_pres){

                                    tiny_bmp_new[k++] = prev_val_pres;
                                    tiny_bmp_new[k++] = j;
                                }

                                --cardinality_cpy;
                                prev_val_pres = j;
                            }
                        }
                    }

                    tiny_bmp_new[k] = prev_val_pres;
                    tiny_bmp_new[0] = (new_cardinality << 8) | rle_list_mode | idx;

                    delete[] tiny_bmp;
                    tiny_bmp = tiny_bmp_new;

                    return sz - bloc_sizes[idx];
                }
            }
            else {

                uint32_t nb_run = 1;

                for (uint32_t i = 2; i <= cardinality; ++i) nb_run += (tiny_bmp[i] != tiny_bmp[i-1]);

                const uint32_t new_cardinality = nb_run << 1;
                const uint32_t new_sz = rndup(new_cardinality + 1);

                if ((new_sz < sz) && (new_sz <= bloc_sizes[nb_bloc_sizes - 1])){

                    uint32_t idx = 0, j = 2;

                    while (new_sz > bloc_sizes[idx]) ++idx;

                    uint32_t* tiny_bmp_new = new uint32_t[bloc_sizes[idx]]();

                    tiny_bmp_new[0] = (new_cardinality << 8) | rle_list_mode | idx;
                    tiny_bmp_new[1] = tiny_bmp[1];

                    for (uint32_t i = 2; i <= cardinality; ++i){

                        if (tiny_bmp[i] != tiny_bmp[i-1]){

                            tiny_bmp_new[j++] = tiny_bmp[i-1];
                            tiny_bmp_new[j++] = tiny_bmp[i];
                        }
                    }

                    tiny_bmp_new[j] = tiny_bmp[cardinality];

                    delete[] tiny_bmp;
                    tiny_bmp = tiny_bmp_new;

                    return sz - bloc_sizes[idx];
                }
            }
        }
    }

    return 0;
}

void TinyBitmap::print() const {

    if (tiny_bmp == nullptr) return;

    const uint32_t mode = tiny_bmp[0] & mode_mask;
    const uint32_t sz = bloc_sizes[tiny_bmp[0] & sz_mul_mask];
    const uint32_t cardinality = tiny_bmp[0] >> 8;

    if (mode == list_mode){

        for (size_t i = 1; i <= cardinality; ++i) cout << tiny_bmp[i] << endl;
    }
    else if (mode == rle_list_mode){

        for (size_t i = 1; i <= cardinality; i += 2){

            for (uint32_t j = tiny_bmp[i]; j <= tiny_bmp[i+1]; ++j) cout << j << endl;
        }
    }
}

bool TinyBitmap::increase_sz(const uint32_t sz_min) {

    //std::cout << "switch_mode(" << sz_min << ")" << std::endl;

    if (sz_min > bloc_sizes[nb_bloc_sizes - 1]) return false;

    const bool is_allocated = (tiny_bmp != nullptr);

    uint32_t idx = is_allocated ? (tiny_bmp[0] & sz_mul_mask) : 0;
    const uint32_t sz = bloc_sizes[idx];

    if (is_allocated && (sz >= sz_min)) return true;

    while (bloc_sizes[idx] < sz_min) ++idx;

    if (is_allocated){

        uint32_t* tiny_bmp_new = new uint32_t[bloc_sizes[idx]]();

        std::copy(tiny_bmp, tiny_bmp + sz, tiny_bmp_new);
        delete[] tiny_bmp;

        tiny_bmp = tiny_bmp_new;
        tiny_bmp[0] = (tiny_bmp[0] & ~sz_mul_mask) | idx;
    }
    else {

        tiny_bmp = new uint32_t[bloc_sizes[idx]]();
        tiny_bmp[0] = bmp_mode | idx;
    }

    return true;
}

bool TinyBitmap::switch_mode(const uint32_t sz_min, const uint32_t new_mode) {

    if (tiny_bmp == nullptr) return true;

    const uint32_t mode = tiny_bmp[0] & mode_mask;
    const uint32_t sz = bloc_sizes[tiny_bmp[0] & sz_mul_mask];

    uint32_t cardinality = tiny_bmp[0] >> 8;

    uint32_t* tiny_bmp_new = nullptr;

    if ((mode == bmp_mode) && (new_mode == list_mode)) { // Switch to list mode

        const uint32_t new_sz_min = std::max(sz_min, cardinality + 1);

        if (new_sz_min > bloc_sizes[nb_bloc_sizes - 1]) return false;

        tiny_bmp_new = new uint32_t[bloc_sizes[0]]();

        std::swap(tiny_bmp, tiny_bmp_new);
        increase_sz(new_sz_min);

        for (size_t i_bmp = 1, i_list = 1; (i_bmp < sz) && (cardinality != 0); ++i_bmp){

            for (size_t j = (i_bmp - 1) * 32; tiny_bmp_new[i_bmp] != 0; tiny_bmp_new[i_bmp] >>= 1, ++j){

                if ((tiny_bmp_new[i_bmp] & 0x1) != 0){

                    tiny_bmp[i_list++] = j;
                    --cardinality;
                }
            }
        }

        tiny_bmp[0] = (tiny_bmp[0] & sz_mul_mask) | list_mode | (tiny_bmp_new[0] & 0xFFFFFF00);
    }
    else if ((mode == list_mode) && (new_mode == bmp_mode)) { // mode is list

        const uint32_t max_val = maximum() + 1;
        const uint32_t new_sz_min = std::max(sz_min, (max_val >> 5) + ((max_val & 0x1F) != 0) + 1);

        if (new_sz_min > bloc_sizes[nb_bloc_sizes - 1]) return false;

        tiny_bmp_new = new uint32_t[bloc_sizes[0]]();

        std::swap(tiny_bmp, tiny_bmp_new);
        increase_sz(new_sz_min);

        for (size_t i = 1; i <= cardinality; ++i) tiny_bmp[(tiny_bmp_new[i] >> 5) + 1] |= (1U << (tiny_bmp_new[i] & 0x1F));

        tiny_bmp[0] = (tiny_bmp[0] & sz_mul_mask) | bmp_mode | (tiny_bmp_new[0] & 0xFFFFFF00);
    }
    else if (mode == rle_list_mode){

        if (new_mode == bmp_mode){

            uint32_t new_card = 0;

            const uint32_t max_val = maximum() + 1;
            const uint32_t new_sz_min = std::max(sz_min, (max_val >> 5) + ((max_val & 0x1F) != 0) + 1);

            if (new_sz_min > bloc_sizes[nb_bloc_sizes - 1]) return false;

            tiny_bmp_new = new uint32_t[bloc_sizes[0]]();

            std::swap(tiny_bmp, tiny_bmp_new);
            increase_sz(new_sz_min);

            for (size_t i = 1; i <= cardinality; i += 2){

                for (uint32_t j = tiny_bmp_new[i]; j <= tiny_bmp_new[i+1]; ++j, ++new_card) tiny_bmp[(j >> 5) + 1] |= (1U << (j & 0x1F));
            }

            tiny_bmp[0] = (tiny_bmp[0] & sz_mul_mask) | bmp_mode | (new_card << 8);
        }
        else if (new_mode == list_mode){

            const uint32_t new_card = size();
            const uint32_t new_sz_min = std::max(sz_min, new_card + 1);

            if (new_sz_min > bloc_sizes[nb_bloc_sizes - 1]) return false;

            tiny_bmp_new = new uint32_t[bloc_sizes[0]]();

            std::swap(tiny_bmp, tiny_bmp_new);
            increase_sz(new_sz_min);

            for (size_t i = 1, k = 1; i <= cardinality; i += 2){

                for (uint32_t j = tiny_bmp_new[i]; j <= tiny_bmp_new[i+1]; ++j, ++k) tiny_bmp[k] = j;
            }

            tiny_bmp[0] = (tiny_bmp[0] & sz_mul_mask) | list_mode | (new_card << 8);
        }
    }

    delete[] tiny_bmp_new;

    return true;
}

const uint32_t TinyBitmap::sz_mul_mask = 0xF; // 00001111
const uint32_t TinyBitmap::mode_mask = 0xF0; // 11110000

const uint32_t TinyBitmap::bmp_mode = 0x00;
const uint32_t TinyBitmap::list_mode = 0x10;
const uint32_t TinyBitmap::rle_list_mode = 0x20;

const uint32_t TinyBitmap::bloc_sz = 4; // uint32_t[4]
const uint32_t TinyBitmap::bloc_sz_bits = bloc_sz * 32; // uint32_t[4] is 128 bits
const uint32_t TinyBitmap::nb_blocks_max = 512; // 64 * uint32_t[4] is 8 kB

const uint32_t TinyBitmap::bloc_sizes[] = { bloc_sz, 2 * bloc_sz, 4 * bloc_sz, 8 * bloc_sz,
                                            16 * bloc_sz, 32 * bloc_sz, 64 * bloc_sz, 128 * bloc_sz,
                                            256 * bloc_sz, 512 * bloc_sz, 0xFFFFFFFF};

const uint32_t TinyBitmap::nb_bloc_sizes = 10;
