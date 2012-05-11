#ifndef BFG_CONTIG_HPP
#define BFG_CONTIG_HPP

#include "Common.hpp"
#include "Kmer.hpp"
#include "CompressedSequence.hpp"


class Contig {
public:
  Contig() : cov(0) {}
  Contig(const char *s) : seq(s),cov(0) {}
  ~Contig();

  uint32_t *cov;
  CompressedSequence seq;
  // TODO: do we store the links here?
};

#endif // BFG_CONTIG_H
