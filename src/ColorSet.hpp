#ifndef BFG_COLOR_HPP
#define BFG_COLOR_HPP

#include <roaring/roaring.hh>

#include "CompactedDBG.hpp"

/** @file src/ColorSet.hpp
* Interface for UnitigColors, the unitig container of k-mer color sets used in ColoredCDBG.
* Code snippets using this interface are provided in snippets/test.cpp.
*/

template<typename Unitig_data_t> class ColoredCDBG;
template<typename Unitig_data_t> class DataAccessor;
template<typename Unitig_data_t> class DataStorage;

template<typename U> using UnitigColorMap = UnitigMap<DataAccessor<U>, DataStorage<U>>;
template<typename U> using const_UnitigColorMap = const_UnitigMap<DataAccessor<U>, DataStorage<U>>;

/** @class UnitigColors
* @brief Represent the k-mer color sets of a unitig. Its template parameter is the type of data
* associated with the unitigs of the graph.
*/
template<typename Unitig_data_t = void>
class UnitigColors {

        typedef Roaring Bitmap;
        typedef Unitig_data_t U;

        template<typename U> friend class ColoredCDBG;
        template<typename U> friend class DataAccessor;
        template<typename U> friend class DataStorage;

    public:

        /** @class UnitigColors_const_iterator
        * @brief See UnitigColors::const_iterator
        */
        template<typename U>
        class UnitigColors_const_iterator : public std::iterator<std::forward_iterator_tag, size_t> {

            template<typename T> friend class UnitigColors;

            private:

                const UnitigColors<U>* cs;

                size_t flag;
                size_t it_setBits;
                size_t color_id;
                size_t cs_sz;

                const Roaring empty_roar;
                Roaring::const_iterator it_roar;

                UnitigColors_const_iterator(const UnitigColors<U>* cs_, const bool beg) :   cs(cs_), it_setBits(0xffffffffffffffff),
                                                                                            color_id(0), it_roar(empty_roar.end()) {

                    flag = cs->setBits & flagMask;

                    if (flag == ptrCompressedBitmap){

                        it_roar = beg ? cs->getConstPtrBitmap()->begin() : cs->getConstPtrBitmap()->end();
                        cs_sz = cs->size();
                    }
                    else cs_sz = (flag == localSingleColor) ? 1 : maxBitVectorIDs;

                    if (!beg) it_setBits = cs_sz;
                }

            public:

                /** Constructor of an empty iterator. The resulting iterator can't be used as it is
                * (it is not associated with any UnitigColors).
                */
                UnitigColors_const_iterator() : cs(nullptr), flag(unoccupied), it_setBits(0),
                                                color_id(0), cs_sz(0), it_roar(empty_roar.end()) {}

                /** Copy assignment operator. After the call to this function, the same iterator exists
                * twice in memory.
                * @param o is the iterator to copy.
                * @return a reference to the iterator copy.
                */
                UnitigColors_const_iterator<U>& operator=(const UnitigColors_const_iterator& o) {

                    cs = o.cs;
                    flag = o.flag;
                    it_setBits = o.it_setBits;
                    color_id = o.color_id;
                    cs_sz = o.cs_sz;
                    it_roar = o.it_roar;

                    return *this;
                }

                /** Indirection operator.
                * @return the current color ID of the iterator.
                */
                size_t operator*() const { return color_id; }

                /** Postfix increment operator, iterate over the next color.
                * @return an iterator on the color ID pointed to before the call to this function.
                */
                UnitigColors_const_iterator<U> operator++(int) {

                    UnitigColors_const_iterator<U> tmp(*this);
                    operator++();
                    return tmp;
                }

                /** Prefix increment operator, iterate over the next color.
                * @return an iterator on the next color ID.
                */
                UnitigColors_const_iterator<U>& operator++() {

                    if (it_setBits == cs_sz) return *this;

                    ++it_setBits;

                    if (flag == ptrCompressedBitmap) {

                        if (it_setBits != 0) ++it_roar;
                        if (it_roar != cs->getConstPtrBitmap()->end()) color_id = *it_roar;
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
                    else if (flag == localSingleColor) color_id = cs->setBits >> 2;

                    return *this;
                }

                /** Equality operator.
                * @return a boolean indicating if two iterators are the same (true) or not (false).
                */
                bool operator==(const UnitigColors_const_iterator& o) const {

                    return  (cs == o.cs) && (flag == o.flag) && (cs_sz == o.cs_sz) &&
                            ((flag == ptrCompressedBitmap) ? (it_roar == o.it_roar) : (it_setBits == o.it_setBits));
                }

                /** Inequality operator.
                * @return a boolean indicating if two iterators are different (true) or not (false).
                */
                bool operator!=(const UnitigColors_const_iterator& o) const { return !operator==(o); }
        };

        /** @typedef const_iterator
        * @brief Iterator for the colors of a unitig. Each color the
        * iterator returns appears on AT LEAST one k-mer of the unitig.
        */
        typedef UnitigColors_const_iterator<U> const_iterator;

        /** Constructor (set up an empty container of k-mer color sets). The UnitigColors is
        * initialized as "unoccupied": the color set is "free" to be used, it is not associated
        * with a unitig.
        */
        UnitigColors();

        /** Copy constructor. After the call to this constructor, the same UnitigColors exists
        * twice in memory.
        * @param o is the color set to copy.
        */
        UnitigColors(const UnitigColors& o); // Copy constructor

        /** Move constructor. After the call to this constructor, the UnitigColors to move is empty
        * (set as "unoccupied") and its content has been transfered (moved) to a new UnitigColors.
        * @param o is the color set to move.
        */
        UnitigColors(UnitigColors&& o); // Move  constructor

        /** Destructor.
        */
        ~UnitigColors();

        /** Copy assignment operator. After the call to this operator, the same UnitigColors exists
        * twice in memory.
        * @param o is the UnitigColors to copy.
        * @return a reference to the current UnitigColors which is a copy of o.
        */
        UnitigColors& operator=(const UnitigColors& o);

        /** Move assignment operator. After the call to this operator, the UnitigColors to move is empty
        * (set as "unoccupied") and its content has been transfered (moved) to another UnitigColors.
        * @param o is the UnitigColors to move.
        * @return a reference to the current UnitigColors having (owning) the content of o.
        */
        UnitigColors& operator=(UnitigColors&& o);

        /** Empty a UnitigColors of its content. If the UnitigColors is associated with a
        * unitig, it is still the case (it is NOT set to "unoccupied").
        */
        void empty();

        /** Add a color in the current UnitigColors to all k-mers of a unitig mapping.
        * @param um is a const_UnitigColorMap representing a unitig mapping for which the color must be added.
        * The color will be added only for the mapping, i.e, unitig[um.dist..um.dist+um.len+k-1]
        * @param color_id is the ID of the color to add.
        */
        void add(const const_UnitigColorMap<U>& um, const size_t color_id);

        /** Check if a color is present on all k-mers of a unitig mapping.
        * @param um is a const_UnitigColorMap representing a unitig mapping. All k-mers of this mapping will be
        * checked for the presence of the color. If true is returned, all k-mers of the mapping have the color.
        * If false is returned, at least one k-mer of the mapping does not have the color.
        * @param color_id is the ID of the color to check the presence.
        * @return a boolean indicating if the color is present on all k-mers of the unitig mapping.
        */
        bool contains(const const_UnitigColorMap<U>& um, const size_t color_id) const;

        /** Get the number of colors of a unitig (sum of the number of colors for each k-mer of the unitig).
        * @return Number of colors (sum of the number of colors for each k-mer of the unitig).
        */
        size_t size() const;

        /** Write a UnitigColors to a stream.
        * @param stream_out is an out stream to which the UnitigColors must be written. It must be
        * opened prior to the call of this function and it won't be closed by this function.
        * @return a boolean indicating if the write was successful.
        */
        bool write(ostream& stream_out) const;

        /** Read a UnitigColors from a stream.
        * @param stream_in is an in stream from which the UnitigColors must be read. It must be
        * opened prior to the call of this function and it won't be closed by this function.
        * @return a boolean indicating if the write was successful.
        */
        bool read(istream& stream_in);

        /** Size of the UnitigColors in bytes.
        * @return Size of the UnitigColors in bytes.
        */
        size_t getSizeInBytes() const;

        /** Optimize the memory of a UnitigColor. Useful in case stretches of overlapping k-mers
        * in a unitig share the same color.
        */
        inline void optimize(){

            if ((setBits & flagMask) == ptrCompressedBitmap) getPtrBitmap()->runOptimize();
        }

        /** Create a constant iterator to the first color of the UnitigColors. Each color the
        * iterator returns appears on AT LEAST one k-mer of the unitig.
        * @return a constant iterator to the first color of the UnitigColors. Each color the
        * iterator returns appears on AT LEAST one k-mer of the unitig.
        */
        const_iterator begin() const {

            const_iterator it(this, true);
            return ++it;
        }

        /** Create a constant iterator to the "past-the-last" color of the UnitigColors.
        * @return a constant iterator to the "past-the-last" color of the UnitigColors.
        */
        const_iterator end() const { return const_iterator(this, false); }

        //void test(const const_UnitigColorMap<U>& um);

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
        void merge(const UnitigColors<U>& cs);

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
