#ifndef BFG_COLOR_HPP
#define BFG_COLOR_HPP

#include "roaring.hh"

#include "CompactedDBG.hpp"
#include "TinyBitmap.hpp"

/** @file src/ColorSet.hpp
* Interface for UnitigColors, the unitig container of k-mer color sets used in ColoredCDBG.
* Code snippets using this interface are provided in snippets/test.cpp.
*/

template<typename Unitig_data_t> class ColoredCDBG;
template<typename Unitig_data_t> class DataAccessor;
template<typename Unitig_data_t> class DataStorage;

/** @class UnitigColors
* @brief Represent the k-mer color sets of a unitig.
*/
class UnitigColors {

        typedef Roaring Bitmap;

        template<typename U> friend class ColoredCDBG;
        template<typename U> friend class DataAccessor;
        template<typename U> friend class DataStorage;

    public:

        class UnitigColors_const_iterator;

        /** @class ColorKmer_ID
        * @brief Represents the position of a k-mer in a unitig and the ID of the color associated
        * with the k-mer at the given position.
        */
        class ColorKmer_ID {

            friend class UnitigColors_const_iterator;
            friend class UnitigColors;

            public:

                /** Get the color ID.
                * @param length_unitig_km is the length of the unitig (in k-mers) associated with the
                * UnitigColors object this ColorKmer_ID object refers to.
                * @return a color ID;
                */
                size_t getColorID(const size_t length_unitig_km) const;

                /** Get the k-mer position (on the forward strand of the unitig).
                * @param length_unitig_km is the length of the unitig (in k-mers) associated with the
                * UnitigColors object this ColorKmer_ID object refers to.
                * @return a k-mer position;
                */
                size_t getKmerPosition(const size_t length_unitig_km) const;

            private:

                size_t ck_id;

                ColorKmer_ID();
                ColorKmer_ID(const size_t id);
                ColorKmer_ID& operator=(const size_t ck_id_);
        };

        /** @class UnitigColors_const_iterator
        * @brief See UnitigColors::const_iterator
        */
        class UnitigColors_const_iterator : public std::iterator<std::forward_iterator_tag, ColorKmer_ID> {

            friend class UnitigColors;

            public:

                /** Constructor of an empty iterator. The resulting iterator can't be used as it is
                * (it is not associated with any UnitigColors).
                */
                UnitigColors_const_iterator();

                /** Copy assignment operator. After the call to this function, the same iterator exists
                * twice in memory.
                * @param o is the iterator to copy.
                * @return a reference to the iterator copy.
                */
                UnitigColors_const_iterator& operator=(const UnitigColors_const_iterator& o);

                /** Indirection operator.
                * @return a ColorKmer_ID object reference representing the position of a k-mer in the unitig
                * and the ID of the color associated with the k-mer at the given position.
                */
                const ColorKmer_ID& operator*() const;

                /** Dereference operator.
                * @return a ColorKmer_ID object pointer representing the position of a k-mer in the unitig
                * and the ID of the color associated with the k-mer at the given position.
                */
                const ColorKmer_ID* operator->() const;

                /** Postfix increment operator: it iterates over the next k-mer of the unitig having the
                * current color or the first k-mer having the next color if all k-mers having the current
                * color have already been visited by this iterator.
                * @return a copy of the iterator before the call to this operator.
                */
                UnitigColors_const_iterator operator++(int);

                /** Prefix increment operator: it iterates over the next k-mer of the unitig having the
                * current color or the first k-mer having the next color if all k-mers having the current
                * color have already been visited by this iterator.
                * @return a reference to the current iterator.
                */
                UnitigColors_const_iterator& operator++();

                /** Color increment operator: it iterates over the first k-mer position of the next color.
                * @return a reference to the current iterator.
                */
                UnitigColors_const_iterator& nextColor(const size_t length_unitig_km);

                /** Equality operator.
                * @return a boolean indicating if two iterators are the same (true) or not (false).
                */
                bool operator==(const UnitigColors_const_iterator& o) const;

                /** Inequality operator.
                * @return a boolean indicating if two iterators are different (true) or not (false).
                */
                bool operator!=(const UnitigColors_const_iterator& o) const;

            private:

                const UnitigColors* cs;

                size_t flag;

                size_t it_setBits;
                size_t cs_sz;

                ColorKmer_ID ck_id;

                const Roaring empty_roar;

                Roaring::const_iterator it_roar;
                TinyBitmap::const_iterator it_t_bmp;

                UnitigColors_const_iterator(const UnitigColors* cs_, const bool beg);
        };

        /** @typedef const_iterator
        * @brief Iterator for the colors of a unitig. The iterator iterates over the colors of a
        * unitig in ascending order. For each color, it iterates over the k-mer positions
        * of the unitig k-mers having this color, in ascending order.
        */
        typedef UnitigColors_const_iterator const_iterator;

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

        /** Check if a UnitigColors is empty (no colors).
        * @return a boolean indicating if the UnitigColors is empty.
        */
        inline bool isEmpty() const { return (size() == 0); }

        /** Add a color in the current UnitigColors to all k-mers of a unitig mapping.
        * @param um is a UnitigMapBase object representing a mapping to a reference unitig for which the
        * color must be added. The color will be added only for the given mapping, i.e, unitig[um.dist..um.dist+um.len+k-1]
        * @param color_id is the ID of the color to add.
        */
        void add(const UnitigMapBase& um, const size_t color_id);

        /** Remove a color in the current UnitigColors for all k-mers of a unitig mapping.
        * @param um is a UnitigMapBase object representing a mapping to a reference unitig for which the
        * color must be removed. The color will be removed only for the given mapping, i.e, unitig[um.dist..um.dist+um.len+k-1]
        * @param color_id is the ID of the color to remove.
        */
        void remove(const UnitigMapBase& um, const size_t color_id);

        /** Check if a color is present on all k-mers of a unitig mapping.
        * @param um is a UnitigMapBase object representing a mapping to a reference unitig. All k-mers of this mapping will be
        * checked for the presence of the color. If true is returned, all k-mers of the mapping have the color.
        * If false is returned, at least one k-mer of the mapping does not have the color.
        * @param color_id is the ID of the color to check the presence.
        * @return a boolean indicating if the color is present on all k-mers of the unitig mapping.
        */
        bool contains(const UnitigMapBase& um, const size_t color_id) const;

        /** Get the number of colors of a unitig (sum of the number of colors for each k-mer of the unitig).
        * @return Number of colors (sum of the number of colors for each k-mer of the unitig).
        */
        size_t size() const;

        /** Get the number of k-mers of a unitig having a given color.
        * @param um is a UnitigMapBase object representing a mapping to a reference unitig
        * @param color_id is the color index
        * @return Number of colors (sum of the number of colors for each k-mer of the unitig).
        */
        size_t size(const UnitigMapBase& um, const size_t color_id) const;

        ColorKmer_ID maximum() const;

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

        /** Create a constant iterator to the first color of the UnitigColors. Each color the
        * iterator returns appears on AT LEAST one k-mer of the unitig.
        * @return a constant iterator to the first color of the UnitigColors. Each color the
        * iterator returns appears on AT LEAST one k-mer of the unitig.
        */
        const_iterator begin() const;

        /** Create a constant iterator to the "past-the-last" color of the UnitigColors.
        * @return a constant iterator to the "past-the-last" color of the UnitigColors.
        */
        const_iterator end() const;

    private:

        inline void releaseMemory(){

            const uintptr_t flag = setBits & flagMask;

            if (flag == ptrBitmap) delete getPtrBitmap();
            else if (flag == localTinyBitmap) t_bmp.empty();

            setBits = localBitVector;
        }

        inline bool isBitmap() const { return ((setBits & flagMask) == ptrBitmap); }
        inline bool isTinyBitmap() const { return ((setBits & flagMask) == localTinyBitmap); }

        void add(const size_t color_id);
        bool contains(const size_t color_km_id) const;

        UnitigColors reverse(const UnitigMapBase& um) const;

        inline Bitmap* getPtrBitmap() const { return reinterpret_cast<Bitmap*>(setBits & pointerMask); }
        inline const Bitmap* getConstPtrBitmap() const { return reinterpret_cast<const Bitmap*>(setBits & pointerMask); }

        static const size_t maxBitVectorIDs; // 64 bits - 2 bits for the color set type = 62

        // asBits and asPointer represent:
        // Flag 0 - A TinyBitmap which can contain up to 65488 uint
        // Flag 1 - A bit vector of 62 bits storing presence/absence of up to 62 integers
        // Flag 2 - A single integer
        // Flag 3 - A pointer to a CRoaring compressed bitmap which can contain up to 2^32 uint

        static const uintptr_t localTinyBitmap; // Flag 0
        static const uintptr_t localBitVector; // Flag 1
        static const uintptr_t localSingleInt; // Flag 2
        static const uintptr_t ptrBitmap; // Flag 3

        static const uintptr_t flagMask; // 0x3
        static const uintptr_t pointerMask; // 0xfffffffffffffffc

        union {

            uintptr_t setBits;
            TinyBitmap t_bmp;
            Bitmap* setPointer;
        };
};

#endif
