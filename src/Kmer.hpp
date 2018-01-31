#ifndef BFG_KMER_HPP
#define BFG_KMER_HPP

#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 32
#endif

#include <stdio.h>
#include <stdint.h>
#include <cassert>
#include <cstring>
#include <string>
#include <iostream>

#ifndef XXH_NAMESPACE
#define XXH_NAMESPACE BIFROST_HASH_
#endif

extern "C" {
    #include "xxhash.h"
}

class CompressedSequence;

/* Short description:
 *  - Store kmer strings by using 2 bits per base instead of 8
 *  - Easily return reverse complements of kmers, e.g. TTGG -> CCAA
 *  - Easily compare kmers
 *  - Provide hash of kmers
 *  - Get last and next kmer, e.g. ACGT -> CGTT or ACGT -> AACGT
 *  */
class Kmer {

    public:

        friend class CompressedSequence;

        Kmer();
        Kmer(const Kmer& o);
        explicit Kmer(const char *s);

        Kmer& operator=(const Kmer& o);

        void set_empty();
        void set_deleted();
        void set_kmer(const char *s);

        bool operator<(const Kmer& o) const;

        bool operator==(const Kmer& o) const;
        bool operator!=(const Kmer& o) const;

        inline uint64_t hash(const uint64_t seed = 0) const {
            return (uint64_t)XXH64((const void *)bytes, MAX_K/4, seed);
        }

        Kmer twin() const;
        Kmer rep() const;

        Kmer getLink(const size_t index) const;

        Kmer forwardBase(const char b) const;
        Kmer backwardBase(const char b) const;

        void selfForwardBase(const char b);

        std::string getBinary() const;
        char getChar(const size_t offset) const;

        void toString(char *s) const;
        std::string toString() const;

        bool write(std::ostream& stream_out) const;

        // static functions
        static void set_k(unsigned int _k);

        static const unsigned int MAX_K = MAX_KMER_SIZE;
        static unsigned int k;

    private:

        // data fields
        union {

            uint8_t bytes[MAX_K/4];
            uint64_t longs[MAX_K/32];
        };

        //static unsigned int k_bytes;
        //static unsigned int k_longs;
        //static unsigned int k_modmask; // int?

        // By default MAX_K == 64 so the union uses 16 bytes
        // However sizeof(Kmer) == 24
        // Are the 8 extra bytes alignment?

        // private functions
        //void shiftForward(int shift);

        //void shiftBackward(int shift);
};


struct KmerHash {

    size_t operator()(const Kmer& km) const {

        return km.hash();
    }
};




//template<typename T> class MinimizerHashTable;

class Minimizer {

    //template<typename T> friend class MinimizerHashTable;

    public:

        Minimizer();
        Minimizer(const Minimizer& o);
        explicit Minimizer(const char *s);

        Minimizer& operator=(const Minimizer& o);

        void set_empty();
        void set_deleted();

        bool operator<(const Minimizer& o) const;
        bool operator==(const Minimizer& o) const;
        bool operator!=(const Minimizer& o) const;

        void set_minimizer(const char *s);

        inline uint64_t hash(const uint64_t seed = 0) const {

            return (uint64_t)XXH64((const void *)bytes, MAX_G/4, seed);
        }

        Minimizer twin() const;
        Minimizer rep() const;

        Minimizer getLink(const size_t index) const;

        Minimizer forwardBase(const char b) const;
        Minimizer backwardBase(const char b) const;

        std::string getBinary() const;

        void toString(char *s) const;
        std::string toString() const;

        // static functions
        static void set_g(unsigned int _g);

        static const unsigned int MAX_G = MAX_KMER_SIZE;
        static unsigned int g;

    private:

        // data fields
        union {

            uint8_t bytes[MAX_G/4];
            uint64_t longs[MAX_G/32];
        };

        //static unsigned int g_bytes;
        //static unsigned int g_longs;
        //static unsigned int g_modmask; // int?

        // By default MAX_K == 64 so the union uses 16 bytes
        // However sizeof(Kmer) == 24
        // Are the 8 extra bytes alignment?

        // private functions
        //void shiftForward(int shift);

        //void shiftBackward(int shift);
};


struct MinimizerHash {

    size_t operator()(const Minimizer& minz) const {

        return minz.hash();
    }
};

#endif // BFG_KMER_HPP
