#ifndef BFG_TINYBITMAP_HPP
#define BFG_TINYBITMAP_HPP

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdint.h>
#include <vector>

/* TinyBitmap is a compressed bitmap that mimics the behavior of a CRoaring container.
*  Its main purpose is to store a tiny set of unsigned integers, up to 65488 uint.
*  Main differences with a CRoaring bitmap are:
*  + For one value inserted, CRoaring allocates >100 bytes, TinyBitmap allocates 24.
*  + CRoaring containers have a variable size except in bitmap mode (fixed 8kB). TinyBitmap
*    container have variable size in all modes (bitmap, list, RLE).
*  + Accessing a value in a CRoaring container is >3 cache-miss. TinyBitmap is 1 cache-miss.
*  - CRoaring can store >65488 values.
*  - CRoaring is SIMD optimized.
*  - CRoaring has a lot more functions (set intersection, union, etc.).
*/

using namespace std;

class TinyBitmap {

    public:

        TinyBitmap();
        TinyBitmap(const TinyBitmap& o);
        TinyBitmap(TinyBitmap&& o);

        ~TinyBitmap();

        TinyBitmap& operator=(const TinyBitmap& o);
        TinyBitmap& operator=(TinyBitmap&& o);

        void empty();

        bool add(const uint32_t val);
        void remove(const uint32_t val);
        bool contains(const uint32_t val) const;

        uint32_t maximum() const;

        size_t getSizeInBytes() const;
        size_t size() const;

        size_t runOptimize();

        static bool test(const bool verbose = true);

    private:

        void print() const;

        bool change_sz(const uint16_t sz_min);
        bool switch_mode(const uint16_t sz_min, const uint16_t new_mode);

        inline uint16_t getSize() const { return (tiny_bmp[0] & sz_mask) >> 3; }
        inline uint16_t getMode() const { return (tiny_bmp[0] & mode_mask); }
        inline uint16_t getBits() const { return (tiny_bmp[0] & bits_mask); }

        inline uint16_t getCardinality() const { return tiny_bmp[1]; }
        inline uint16_t getOffset() const { return tiny_bmp[2]; }

        static inline uint16_t getNextSize(const uint16_t sz) {

            uint16_t idx = 0;

            while (sizes[idx] < sz) ++idx;

            return sizes[idx];
        }

        static const uint16_t sz_mask;
        static const uint16_t mode_mask;
        static const uint16_t bits_mask;

        static const uint16_t bmp_mode;
        static const uint16_t list_mode;
        static const uint16_t rle_list_mode;

        static const uint16_t bits_16;
        static const uint16_t bits_32;

        static const uint16_t sizes[];
        static const uint16_t nb_sizes;

        uint16_t* tiny_bmp;
};

#endif
