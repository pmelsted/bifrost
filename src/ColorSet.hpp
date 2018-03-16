#ifndef BFG_COLOR_HPP
#define BFG_COLOR_HPP

#include <roaring/roaring.hh>

#include "CompactedDBG.hpp"

/** @file src/ColorSet.hpp
* Interface for the color sets used in ColoredCDBG.
* Code snippets using this interface are provided in snippets.hpp.
*/

template<typename Unitig_data_t> class ColoredCDBG;
template<typename Unitig_data_t> class DataAccessor;
template<typename Unitig_data_t> class DataStorage;

template<typename U> using UnitigColorMap = UnitigMap<DataAccessor<U>, DataStorage<U>>;
template<typename U> using const_UnitigColorMap = const_UnitigMap<DataAccessor<U>, DataStorage<U>>;

/** @class UnitigColors
* @brief Represent a color set for a unitig. The number of colors in such a set
* , i.e number of k-mers in unitig * number of color per k-mer, can't exceed 2^32
*/
template<typename Unitig_data_t = void>
class UnitigColors {

        typedef Roaring Bitmap;
        typedef Unitig_data_t U;

        template<typename U> friend class ColoredCDBG;
        template<typename U> friend class DataAccessor;
        template<typename U> friend class DataStorage;

    public:

        template<typename U>
        class UnitigColors_const_iterator : public std::iterator<std::forward_iterator_tag, size_t> {

            private:

                const UnitigColors<U>* cs;
                size_t color_id;

                size_t flag;
                size_t it_setBits;

                const Roaring empty_roar;
                Roaring::const_iterator it_roar;

            public:

                UnitigColors_const_iterator() : cs(nullptr), flag(localBitVectorColor), it_setBits(maxBitVectorIDs), color_id(0), it_roar(empty_roar.end()) {}

                UnitigColors_const_iterator(const UnitigColors<U>* cs_) : cs(cs_), it_setBits(0), color_id(0), it_roar(empty_roar.end()) {

                    flag = cs->setBits & flagMask;
                    if (flag == ptrCompressedBitmap) it_roar = cs->getConstPtrBitmap()->begin();
                    else if (flag == unoccupied){

                        flag = localBitVectorColor;
                        it_setBits = maxBitVectorIDs;
                    }
                }

                UnitigColors_const_iterator<U>& operator=(const UnitigColors_const_iterator& o) {

                    cs = o.cs;
                    color_id = o.color_id;
                    flag = o.flag;
                    it_setBits = o.it_setBits;
                    it_roar = o.it_roar;

                    return *this;
                }

                size_t operator*() const { return color_id; }

                UnitigColors_const_iterator<U> operator++(int) {

                    UnitigColors_const_iterator<U> tmp(*this);
                    operator++();
                    return tmp;
                }

                UnitigColors_const_iterator<U>& operator++() {

                    if (flag == ptrCompressedBitmap) {

                        if (it_roar != cs->getConstPtrBitmap()->end()) color_id = *(++it_roar);
                    }
                    else if (flag == localBitVectorColor){

                        while (it_setBits < maxBitVectorIDs){

                            if (((cs->setBits >> (it_setBits + 2)) & 0x1) != 0){

                                color_id = it_setBits;
                                break;
                            }

                            ++it_setBits;
                        }
                    }
                    else if ((flag == localSingleColor) && (it_setBits < 1)){

                        color_id = cs->setBits >> 2;
                        ++it_setBits;
                    }

                    return *this;
                }

                bool operator==(const UnitigColors_const_iterator& o) {

                    return  (cs == o.cs) && (color_id == o.color_id) && (flag == o.flag) &&
                            ((flag == ptrCompressedBitmap) ? (it_roar == o.it_roar) : (it_setBits == o.it_setBits));
                }

                bool operator!=(const UnitigColors_const_iterator& o) {

                    return  (cs != o.cs) || (color_id != o.color_id) || (flag != o.flag) ||
                            ((flag == ptrCompressedBitmap) ? (it_roar != o.it_roar) : (it_setBits != o.it_setBits));
                }
        };

        typedef UnitigColors_const_iterator<U> const_iterator;

        UnitigColors();
        UnitigColors(const UnitigColors& o); // Copy constructor
        UnitigColors(UnitigColors&& o); // Move  constructor

        ~UnitigColors();

        UnitigColors& operator=(const UnitigColors& o); // Copy assignment
        UnitigColors& operator=(UnitigColors&& o); // Move assignment

        void empty();

        void add(const const_UnitigColorMap<U>& um, const size_t color_id);
        bool contains(const const_UnitigColorMap<U>& um, const size_t color_id) const;

        size_t size() const;

        bool write(ostream& stream_out) const;
        bool read(istream& stream_in);

        /** Optimize the memory of a color set. Useful if multiple overlapping k-mers
        * of a unitig share the same color.
        */
        inline void optimize(){

            if ((setBits & flagMask) == ptrCompressedBitmap) getPtrBitmap()->runOptimize();
        }

        /** Create a constant iterator to the first color of the color set.
        * @return a constant iterator to the first color of the color set.
        */
        const_iterator begin() const {

            const_iterator it(this);
            return ++it;
        }

        /** Create a constant iterator to the the "past-the-last" color of the color set.
        * @return a constant iterator to the first the "past-the-last" of the color set.
        */
        const_iterator end() const { return const_iterator(); }

    private:

        inline void releasePointer(){

            if ((setBits & flagMask) == ptrCompressedBitmap) delete getPtrBitmap();
        }

        inline void setOccupied(){ if (isUnoccupied()) setBits = localBitVectorColor; }
        inline void setUnoccupied(){ releasePointer(); setBits = unoccupied; }

        inline bool isUnoccupied() const { return ((setBits & flagMask) == unoccupied); }
        inline bool isOccupied() const { return ((setBits & flagMask) != unoccupied); }

        void add(const size_t color_id);

        UnitigColors<U> reverse(const const_UnitigColorMap<U>& um) const;

        inline Bitmap* getPtrBitmap() const { return reinterpret_cast<Bitmap*>(setBits & pointerMask); }
        inline const Bitmap* getConstPtrBitmap() const { return reinterpret_cast<const Bitmap*>(setBits & pointerMask); }

        static const size_t maxBitVectorIDs = 62; // 64 bits - 2 bits for the color set type

        // asBits and asPointer represent:
        // Flag 0 - A pointer to a compressed bitmap containing colors
        // Flag 1 - A bit vector of 62 bits storing presence/absence of up to 62 colors (one bit = one color)
        // Flag 2 - A single integer which is a color
        // Flag 3 - Unoccupied color set (not associated with any unitig)

        static const uintptr_t ptrCompressedBitmap = 0x0;
        static const uintptr_t localBitVectorColor = 0x1;
        static const uintptr_t localSingleColor = 0x2;
        static const uintptr_t unoccupied = 0x3;

        static const uintptr_t flagMask = 0x3;
        static const uintptr_t pointerMask = 0xfffffffffffffffc;

        union {

            uintptr_t setBits;
            Bitmap* setPointer;
        };
};

#include "ColorSet.tcc"

#endif
