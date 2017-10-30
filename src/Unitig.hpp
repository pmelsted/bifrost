#ifndef BFG_UNITIG_HPP
#define BFG_UNITIG_HPP

#include "Common.hpp"
#include "Kmer.hpp"
#include "CompressedSequence.hpp"
#include "CompressedCoverage.hpp"


/* Short description:
 *  - Use the CompressedSequence class for storing the DNA string
 *  - Use the CompressedCoverage class for storing the kmer coverage
 *  */

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

        size_t memory() const {

            size_t m = sizeof(ccov) + sizeof(seq);
            const size_t numkmers = numKmers();
            const size_t seqlength = length();

            if (numkmers > ccov.size_limit) m += ((numkmers + 3) / 4) + 8;
            if (!seq.isShort()) m += ((seqlength + 3) / 4);

            return m;
        }

        inline size_t numKmers() const { return seq.size( ) -Kmer::k + 1; }
        inline size_t length() const { return seq.size(); }

        inline const T* getData() { return &data; }
        inline void setData(const T* const data_){ data = *data_; }
        inline void releaseData(T* const data){ delete data; }

        uint64_t coveragesum;

        CompressedCoverage ccov;
        CompressedSequence seq;

        T data;
};

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

        size_t memory() const {

            size_t m = sizeof(ccov) + sizeof(seq);
            const size_t numkmers = numKmers();
            const size_t seqlength = length();

            if (numkmers > ccov.size_limit) m += ((numkmers + 3) / 4) + 8;
            if (!seq.isShort()) m += ((seqlength + 3) / 4);

            return m;
        }

        inline size_t numKmers() const { return seq.size( ) - Kmer::k + 1; }
        inline size_t length() const { return seq.size(); }

        inline const void* getData() { return nullptr; }
        inline void setData(const void* const data_){ return; }
        inline void releaseData(void* const data){ return; }

        uint64_t coveragesum;

        CompressedCoverage ccov;
        CompressedSequence seq;
};


#endif // BFG_CONTIG_HPP
