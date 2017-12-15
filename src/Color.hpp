#ifndef BFG_COLOR_HPP
#define BFG_COLOR_HPP

#include "CompactedDBG.hpp"

class Color : public CDBG_Data_t<Color> {

    public:

        Color(const uint8_t c = 0);

        void join(const UnitigMap<Color>& um_dest, const UnitigMap<Color>& um_src);
        void split(const UnitigMap<Color>& um_split, const size_t pos_split, const size_t len_split, Color& new_data) const;

        inline uint8_t getColor() const { return color_id; }
        inline void setColor(const uint8_t col) { color_id = col; }

    private:

        uint8_t color_id;
};

class ColorSet {

    public:

        ColorSet();
        ColorSet(const size_t color);
        ~ColorSet();

        void insertColor(const size_t color);

    private:

        static const size_t maxBitVectorIDs = 62; // 64 bits - 2 bits for the color set type

        // asBits and asPointer represent:
        // Flag 0 - A bit vector of 62 bits storing presence/absence of up to 62 colors (one bit = one color)
        // Flag 1 - A single integer which is a color
        // Flag 2 - A bit vector of 54 bits storing whether the k-mers of the corresponding unitig match different color sets
        // + 8 bits for the color ID of the next k-mer which has a different color set
        // Flag 3 - A pointer to a compressed bitmap containing colors

        static const uintptr_t localBitVectorColor = 0x0;
        static const uintptr_t localSingleColor = 0x1;
        static const uintptr_t localBitVectorKmer = 0x2;
        static const uintptr_t pointerCompressedBitmap = 0x3;

        static const uintptr_t flagMask = 0x3;
        static const uintptr_t pointerMask = 0xfffffffffffffffc;

        union {

            uintptr_t asBits;
            uint8_t* asPointer;
        };
};

#endif
