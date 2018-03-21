#ifndef BFG_UNITIG_HPP
#define BFG_UNITIG_HPP

#include "Common.hpp"
#include "Kmer.hpp"
#include "CompressedSequence.hpp"
#include "CompressedCoverage.hpp"

/** @file src/Unitig.hpp
* The Unitig interface.
* Code snippets using these interface are provided in snippets/test.cpp.
*/

 /** @class Unitig
* @brief Represent a unitig which is a vertex of the Compacted de Bruijn graph.
* The first template argument T is the type of data associated with the unitigs
* of the graph.
* @var Unitig::data
* Data associated with the unitigs of the graph.
*/
template<typename T>
class Unitig {

    public:

        Unitig() : coveragesum(0) {}
        Unitig(const char* s, bool full = false) : seq(s) { initializeCoverage(full); }

        void initializeCoverage(bool full) {

            const size_t ssz = seq.size(), k = Kmer::k;

            coveragesum = 0;

            ccov = CompressedCoverage();
            ccov.initialize(ssz >= k ? ssz - k + 1 : 0, full);
        }

        void cover(size_t start, size_t end) {

            ccov.cover(start, end);

            if (end < start) swap(start, end);

            __sync_add_and_fetch(&coveragesum, end - start + 1);
        }

        inline size_t numKmers() const { return seq.size( ) -Kmer::k + 1; }
        inline size_t length() const { return seq.size(); }

        /** Return a constant pointer to the data associated with the unitig.
        * @return a constant pointer to the data associated with the unitig.
        */
        inline const T* getData() const { return &data; }
        inline T* getData() { return &data; }

        uint64_t coveragesum;

        CompressedCoverage ccov;
        CompressedSequence seq;

        T data;
};

///@cond NO_DOC
template<>
class Unitig<void> {

    public:

        Unitig() : coveragesum(0) {}
        Unitig(const char* s, bool full = false) : seq(s) { initializeCoverage(full); }

        void initializeCoverage(bool full) {

            const size_t ssz = seq.size(), k = Kmer::k;

            coveragesum = 0;

            ccov = CompressedCoverage();
            ccov.initialize(ssz >= k ? ssz - k + 1 : 0, full);
        }

        void cover(size_t start, size_t end) {

            ccov.cover(start, end);

            if (end < start) swap(start, end);

            __sync_add_and_fetch(&coveragesum,end - start + 1);
        }

        inline size_t numKmers() const { return seq.size( ) - Kmer::k + 1; }
        inline size_t length() const { return seq.size(); }

        inline const void* getData() const { return nullptr; }
        inline void* getData() { return nullptr; }

        uint64_t coveragesum;

        CompressedCoverage ccov;
        CompressedSequence seq;
};
///@endcond

#endif // BFG_CONTIG_HPP
