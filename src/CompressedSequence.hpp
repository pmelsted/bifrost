#ifndef BFG_COMPRESSED_SEQUENCE_HPP
#define BFG_COMPRESSED_SEQUENCE_HPP

#include <cstring>
#include <string>
#include <stdint.h>

#include "Kmer.hpp"


/* Short description:
 *  - Compress a DNA string by using 2 bits per base instead of 8
 *  - Easily get the DNA string back from the compressed format
 *  - Create a sequence from a kmer
 *  - Get kmers from a sequence
 *  - Get length of a sequence
 *  - Easily get length of matching substring from a given string
 * */
class CompressedSequence {

    public:

        CompressedSequence();
        ~CompressedSequence();
        CompressedSequence(const CompressedSequence& o);
        CompressedSequence& operator=(const CompressedSequence& o);
        explicit CompressedSequence(const char *s);
        explicit CompressedSequence(const string& s);
        explicit CompressedSequence(const Kmer& km);

        const char operator[](size_t index) const;

        void clear();

        size_t size() const;
        Kmer getKmer(size_t offset) const;
        string toString() const;
        void toString(char *s) const;
        void toString(char *s, size_t offset, size_t length) const;
        string toString(size_t offset, size_t length) const;

        //  void setSequence(const CompressedSequence &o, size_t length, size_t offset = 0, bool reversed=false);
        void setSequence(const CompressedSequence& o, size_t start, size_t length, size_t offset = 0, bool reversed = false);
        void setSequence(const char *s, size_t length, size_t offset = 0, bool reversed=false);
        void setSequence(const string& s, size_t length, size_t offset = 0, bool reversed=false);
        void setSequence(const Kmer& km, size_t length, size_t offset = 0, bool reversed=false);

        void reserveLength(size_t new_length);

        CompressedSequence rev() const;
        size_t jump(const char *s, size_t i, int pos, bool reversed) const;

        bool isShort() const;

    private:

        size_t round_to_bytes(const size_t len) const { return (len+3)/4; }
        void _resize_and_copy(size_t new_cap, size_t copy_limit);
        void initShort();
        void setSize(size_t size);

        size_t capacity() const;
        const char *getPointer() const;

        static const uint8_t shortMask = 1;

        union {

            struct {
                uint32_t _length; // size of sequence
                uint32_t _capacity; // capacity of array allocated in bytes
                char *_data; // 0-based 2bit compressed dna string
                char padding[16];
            } asPointer;

            struct {
                uint8_t _size; // 7 bits can index up to 128
                char _arr[31]; // can store 124 nucleotides
            } asBits;
        };
};

#endif // BFG_COMPRESSED_SEQUENCE_HPP
