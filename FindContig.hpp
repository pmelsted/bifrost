#ifndef BFG_FINDCONTIG_HPP
#define BFG_FINDCONTIG_HPP

#include "Kmer.hpp"
#include "BlockedBloomFilter.hpp"
#include "Common.hpp"

struct FindContig {
  Kmer end;
  size_t dist;
  int selfloop; // 0 for no selfloop, 1 for regular, 2 for reversed
  string s;
  size_t deg;
  FindContig(Kmer km, size_t i, int b, string _s, size_t d) : end(km), dist(i), selfloop(b), s(_s), deg(d) {}
};

bool isNeighbor(Kmer a, Kmer b);
FindContig find_contig_forward(BlockedBloomFilter &bf, Kmer km);

#endif // BFG_FINDCONTIG_HPP
