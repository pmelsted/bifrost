#ifndef BFG_TINYBITMAP_HPP
#define BFG_TINYBITMAP_HPP

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdint.h>

#include "Common.hpp"

class TinyBitmap {

    public:

        TinyBitmap();
        TinyBitmap(const TinyBitmap& o);

        ~TinyBitmap();

        TinyBitmap& operator=(const TinyBitmap& o);

        void empty();

        bool add(const uint32_t val);
        void remove(const uint32_t val);
        bool contains(const uint32_t val) const;

        uint32_t maximum() const;

        size_t getSizeInBytes() const;
        size_t size() const;

        size_t runOptimize();

        void print() const;

    private:

        bool change_sz(const uint16_t sz_min);
        bool switch_mode(const uint16_t sz_min, const uint16_t new_mode);

        inline uint16_t getIndexSize() const { return (tiny_bmp[0] & sz_mask) >> 1; }
        inline uint16_t getMode() const { return (tiny_bmp[0] & mode_mask); }
        inline uint16_t getBits() const { return (tiny_bmp[0] & bits_mask); }

        inline uint16_t getCardinality() const { return tiny_bmp[1]; }
        inline uint16_t getOffset() const { return tiny_bmp[2]; }

        static inline uint16_t getIndexLargerSize(const uint16_t sz_min) {

            uint16_t idx = 0;

            while (sizes[idx] < sz_min) ++idx;

            return idx;
        }

        static const uint16_t sz_mask;
        static const uint16_t mode_mask;
        static const uint16_t bits_mask;

        static const uint16_t bmp_mode;
        static const uint16_t list_mode;
        static const uint16_t rle_list_mode;

        static const uint16_t bits_16;
        static const uint16_t bits_32;

        static const uint16_t sz_min;
        static const uint16_t sz_max;

        static const uint16_t sizes[];
        static const uint16_t nb_sizes;

        uint16_t* tiny_bmp;
};

#endif
