#ifndef BFG_CONTIG_HPP
#define BFG_CONTIG_HPP

#include "Common.hpp"
#include "Kmer.hpp"
#include "CompressedSequence.hpp"


/* Short description: 
 *  - Use the CompressedSequence class for storing the DNA string 
 *  */
class Contig {
public:
  Contig() : cov(0) {}
  Contig(const char *s) : seq(s) { allocateCov(); }
  ~Contig();

  uint8_t *cov;
  uint32_t covlength;
  CompressedSequence seq;
  // TODO: do we store the links here?
};

#endif // BFG_CONTIG_HPP
