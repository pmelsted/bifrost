#ifndef BFG_FINDCONTIG_HPP
#define BFG_FINDCONTIG_HPP

#include "Kmer.hpp"
#include "BloomFilter.hpp"
#include "Common.hpp"

struct FindContig {
  Kmer end;
  size_t dist;
  bool selfloop;
  string s;
  FindContig(Kmer km, size_t i, bool b, string _s) : end(km), dist(i), selfloop(b), s(_s) {}
};

bool isNeighbor(Kmer a, Kmer b);
FindContig find_contig_forward(BloomFilter &bf, Kmer km);

#endif // BFG_FINDCONTIG_HPP
