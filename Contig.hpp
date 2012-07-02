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
  Contig() {}
  Contig(const char *s, bool full=false) : seq(s) { allocateCov(full); }
  ~Contig();
  void allocateCov(bool full);

  //uint8_t *cov;
  uint32_t covlength;
  CompressedCoverage *covp;
  CompressedSequence seq;
  // TODO: do we store the links here?
};

#endif // BFG_CONTIG_HPP
