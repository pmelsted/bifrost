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
  ~Contig();
  void initializeCoverage(bool full);

  uint64_t coveragesum;
  CompressedCoverage ccov;
  CompressedSequence seq;

  size_t numKmers() const { return ccov.size(); }
  size_t length() const { return seq.size(); }
};

#endif // BFG_CONTIG_HPP
