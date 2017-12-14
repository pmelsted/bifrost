#ifndef BFG_UNITIG_HPP
#define BFG_UNITIG_HPP

#include "Common.hpp"
#include "Kmer.hpp"
#include "CompressedSequence.hpp"
#include "CompressedCoverage.hpp"

/** @file src/Unitig.hpp
* The Unitig interface.
* Code snippets using these interface are provided in snippets.hpp.
*/

/* Short description:
 *  - Use the CompressedSequence class for storing the DNA string
 *  - Use the CompressedCoverage class for storing the kmer coverage
 *  */

 /** @class Unitig
* @brief Represent a unitig which is a vertex of the compacted de Bruijn graph.
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

        /*size_t memory() const {

            size_t m = sizeof(ccov) + sizeof(seq);
            const size_t numkmers = numKmers();
            const size_t seqlength = length();

            if (numkmers > ccov.size_limit) m += ((numkmers + 3) / 4) + 8;
            if (!seq.isShort()) m += ((seqlength + 3) / 4);

            return m;
        }*/

        inline size_t numKmers() const { return seq.size( ) -Kmer::k + 1; }
        inline size_t length() const { return seq.size(); }

        /** Return a constant pointer to the data associated with the unitig.
        * @return a constant pointer to the data associated with the unitig.
        */
        inline const T* getData() const { return &data; }

        inline T* getData() { return &data; }

        /** Set the data associated with the unitig.
        * @param data_ is a constant pointer to data that must be copied to the data
        * associated with the current unitig
        */
        inline void setData(const T* const data_){ data = *data_; }

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

        /*size_t memory() const {

            size_t m = sizeof(ccov) + sizeof(seq);
            const size_t numkmers = numKmers();
            const size_t seqlength = length();

            if (numkmers > ccov.size_limit) m += ((numkmers + 3) / 4) + 8;
            if (!seq.isShort()) m += ((seqlength + 3) / 4);

            return m;
        }*/

        inline size_t numKmers() const { return seq.size( ) - Kmer::k + 1; }
        inline size_t length() const { return seq.size(); }

        inline const void* getData() const { return nullptr; }
        inline void* getData() { return nullptr; }
        inline void setData(const void* const data_){ return; }

        uint64_t coveragesum;

        CompressedCoverage ccov;
        CompressedSequence seq;
};
///@endcond

#endif // BFG_CONTIG_HPP
