#ifndef BIFROST_COMPRESSED_SEQUENCE_HPP
#define BIFROST_COMPRESSED_SEQUENCE_HPP

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
        CompressedSequence(const CompressedSequence& o); // Copy constructor
        CompressedSequence(CompressedSequence&& o); // Move constructor

        ~CompressedSequence();

        CompressedSequence& operator=(const CompressedSequence& o); // Copy assignment
        CompressedSequence& operator=(CompressedSequence&& o); // Move assignment

        explicit CompressedSequence(const char *s);
        explicit CompressedSequence(const string& s);
        explicit CompressedSequence(const Kmer& km);

        const char operator[](size_t index) const;

        void clear();

        size_t size() const;

        void toString(char *s, const size_t offset, const size_t length) const;
        string toString(const size_t offset, const size_t length) const;


        inline string toString() const { return toString(0, size()); }
        inline void toString(char *s) const { toString(s, 0, size()); }

        Kmer getKmer(size_t offset) const;
        char getChar(const size_t offset) const;

        //bool compareKmer(const size_t offset, const Kmer& km) const;
        bool compareKmer(const size_t offset, const size_t length, const Kmer& km) const;

        //  void setSequence(const CompressedSequence &o, size_t length, size_t offset = 0, bool reversed=false);
        void setSequence(const CompressedSequence& o, const size_t start, const size_t length, const size_t offset = 0, const bool reversed = false);
        void setSequence(const char *s, const size_t length, const size_t offset = 0, const bool reversed = false);
        void setSequence(const string& s, const size_t length, const size_t offset = 0, const bool reversed=false);
        void setSequence(const Kmer& km, const size_t length, const size_t offset = 0, const bool reversed=false);

        void reserveLength(const size_t new_length);

        CompressedSequence rev() const;

        size_t jump(const char *s, const size_t i, int pos, const bool reversed) const;
        //size_t bw_jump(const char *s, const size_t i, int pos, const bool reversed) const;

        int64_t findKmer(const Kmer& km) const;

        bool isShort() const;

    private:

        inline size_t round_to_bytes(const size_t len) const { return (len+3)/4; }
        void _resize_and_copy(const size_t new_cap, const size_t copy_limit);
        void initShort();
        void setSize(const size_t size);

        size_t capacity() const;
        const unsigned char *getPointer() const;

        static const uint8_t shortMask = 1;

        union {

            struct {
                uint32_t _length; // size of sequence
                uint32_t _capacity; // capacity of array allocated in bytes
                unsigned char *_data; // 0-based 2bit compressed dna string
                unsigned char padding[16];
            } asPointer;

            struct {
                uint8_t _size; // 7 bits can index up to 128
                unsigned char _arr[31]; // can store 124 nucleotides
            } asBits;
        };
};

#endif // BFG_COMPRESSED_SEQUENCE_HPP
