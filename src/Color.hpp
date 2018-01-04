#ifndef BFG_COLOR_HPP
#define BFG_COLOR_HPP

#include <roaring/roaring.hh>

#include "CompactedDBG.hpp"

class Color : public CDBG_Data_t<Color> {

    public:

        Color(const uint8_t c = 0);

        void join(const UnitigMap<Color>& um_dest, const UnitigMap<Color>& um_src);
        void split(const UnitigMap<Color>& um_split, const size_t pos_split, const size_t len_split, Color& new_data) const;

        void getLock();
        inline void releaseLock() { __sync_and_and_fetch(&color_id, 0x7f); }

        inline uint8_t getColorID() const { return color_id & 0x7f; }
        inline void setColorID(const uint8_t col) { color_id = col & 0x7f; }

    private:

        uint8_t color_id;
};

class ColorSet {

    public:

        class ColorSet_const_iterator : public std::iterator<std::forward_iterator_tag, size_t> {

            private:

                const ColorSet* cs;
                size_t color;

                size_t flag;
                size_t it_setBits;

                const Roaring empty_roar;
                Roaring::const_iterator it_roar;

            public:

                ColorSet_const_iterator() : cs(nullptr), flag(localBitVectorColor), it_setBits(maxBitVectorIDs), color(0), it_roar(empty_roar.end()) {}

                ColorSet_const_iterator(const ColorSet* cs_) : cs(cs_), it_setBits(0), color(0), it_roar(empty_roar.end()) {

                    flag = cs->setBits & flagMask;
                    if (flag == pointerCompressedBitmap) it_roar = cs->getConstPointer()->begin();
                }

                ColorSet_const_iterator& operator=(const ColorSet_const_iterator& o) {

                    cs = o.cs;
                    color = o.color;
                    flag = o.flag;
                    it_setBits = o.it_setBits;
                    it_roar = o.it_roar;

                    return *this;
                }

                size_t operator*() const { return color; }

                ColorSet_const_iterator operator++(int) {

                    ColorSet_const_iterator tmp(*this);
                    operator++();
                    return tmp;
                }

                ColorSet_const_iterator& operator++() {

                    if (flag == pointerCompressedBitmap) {

                        if (it_roar != cs->getConstPointer()->end()) color = *(++it_roar);
                    }
                    else if (flag == localBitVectorColor){

                        while (it_setBits < maxBitVectorIDs){

                            if (((cs->setBits >> (it_setBits + 2)) & 0x1) != 0){

                                color = it_setBits;
                                break;
                            }

                            ++it_setBits;
                        }
                    }
                    else if ((flag == localSingleColor) && (it_setBits < 1)){

                        color = cs->setBits >> 2;
                        ++it_setBits;
                    }

                    return *this;
                }

                bool operator==(const ColorSet_const_iterator& o) {

                    return  (cs == o.cs) && (color == o.color) && (flag == o.flag) &&
                            ((flag == pointerCompressedBitmap) ? (it_roar == o.it_roar) : (it_setBits == o.it_setBits));
                }

                bool operator!=(const ColorSet_const_iterator& o) {

                    return  (cs != o.cs) || (color != o.color) || (flag != o.flag) ||
                            ((flag == pointerCompressedBitmap) ? (it_roar != o.it_roar) : (it_setBits != o.it_setBits));
                }
        };

        typedef ColorSet_const_iterator const_iterator;

        ColorSet();
        ColorSet(const size_t color);
        ~ColorSet();

        void add(const size_t color);
        void join(const ColorSet& cs);
        bool contains(const size_t color) const;

        size_t size() const;

        const_iterator begin() const {

            const_iterator it(this);
            return ++it;
        }

        const_iterator end() const { return const_iterator(); }

    private:

        void releasePointer();

        inline Roaring* getPointer() const { return reinterpret_cast<Roaring*>(setBits & pointerMask); }
        inline const Roaring* getConstPointer() const { return reinterpret_cast<const Roaring*>(setBits & pointerMask); }

        static const size_t maxBitVectorIDs = 62; // 64 bits - 2 bits for the color set type

        // asBits and asPointer represent:
        // Flag 0 - A pointer to a compressed bitmap containing colors
        // Flag 1 - A bit vector of 62 bits storing presence/absence of up to 62 colors (one bit = one color)
        // Flag 2 - A single integer which is a color
        // Flag 3 - A bit vector of 54 bits storing whether the k-mers of the corresponding unitig match different color sets
        // + 8 bits for the color ID of the next k-mer which has a different color set

        static const uintptr_t pointerCompressedBitmap = 0x0;
        static const uintptr_t localBitVectorColor = 0x1;
        static const uintptr_t localSingleColor = 0x2;
        static const uintptr_t localBitVectorKmer = 0x3;

        static const uintptr_t flagMask = 0x3;
        static const uintptr_t pointerMask = 0xfffffffffffffffc;

        union {

            uintptr_t setBits;
            Roaring* setPointer;
        };
};

#endif
