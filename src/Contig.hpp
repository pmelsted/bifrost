#ifndef BFG_CONTIG_HPP
#define BFG_CONTIG_HPP

#include "Common.hpp"
#include "Kmer.hpp"
#include "CompressedSequence.hpp"
#include "CompressedCoverage.hpp"


/* Short description:
 *  - Use the CompressedSequence class for storing the DNA string
 *  - Use the CompressedCoverage class for storing the kmer coverage
 *  */
class Contig {
 public:
  Contig() : coveragesum(0) {}
  Contig(const char *s, bool full=false) : seq(s) { initializeCoverage(full); }
  void initializeCoverage(bool full);
  void cover(size_t start, size_t end);

  uint64_t coveragesum;
  CompressedCoverage ccov;
  CompressedSequence seq;

  inline size_t numKmers() const { return seq.size()-Kmer::k+1; }
  inline size_t length() const { return seq.size(); }
  size_t memory() const;
};


#endif // BFG_CONTIG_HPP
