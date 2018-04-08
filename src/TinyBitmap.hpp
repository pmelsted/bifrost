#ifndef BFG_TINYBITMAP_HPP
#define BFG_TINYBITMAP_HPP

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdint.h>

class TinyBitmap {

    public:

        TinyBitmap();
        TinyBitmap(const TinyBitmap& o);

        ~TinyBitmap();

        TinyBitmap& operator=(const TinyBitmap& o);

        bool add(const uint32_t val);
        bool contains(const uint32_t val);

        uint32_t maximum() const;

        size_t getSizeInBytes() const;
        size_t size() const;

    private:

        bool increase_sz(const size_t sz_min);
        bool switch_mode();

        static const uint32_t sz_mul_mask;
        static const uint32_t mode_mask;

        static const uint32_t bloc_sz;
        static const uint32_t bloc_sz_bits;
        static const uint32_t nb_blocks_max;

        static const uint16_t bloc_sizes[];
        static const uint32_t nb_bloc_sizes;

        static const uint32_t bmp_mode;
        static const uint32_t list_mode;

        uint32_t* tiny_bmp;
};

#endif
