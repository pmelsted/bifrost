#ifndef BFG_CONTIGMETHODS_HPP
#define BFG_CONTIGMETHODS_HPP

#include <string>
#include "Common.hpp"
#include "Kmer.hpp"
#include "KmerMapper.hpp"
#include "BloomFilter.hpp"
#include "FindContig.hpp"

/* Structs for Contig and Kmer information */

struct NewContig {
  string seq;
  size_t start;
  size_t end;
  size_t read_index;
  NewContig(string s, size_t i, size_t j, size_t ri) : seq(s), start(i), end(j), read_index(ri) {}
};

struct CheckContig {
  ContigRef cr;
  size_t dist;
  bool repequal;
  CheckContig(ContigRef ref, size_t i, bool eq) : cr(ref), dist(i), repequal(eq) {}
};


struct MakeContig {
  string seq;
  int selfloop; // 0 for no selfloop, 1 for regular, 2 for reversed
  size_t pos;
  MakeContig(string s, int l, size_t p) : seq(s), selfloop(l), pos(p) {}
};

/* Methods for Contig and Kmer information and mapping */

void getMappingInfo(const bool repequal, const int32_t pos, const size_t dist, size_t &kmernum, int32_t &cmppos);
CheckContig check_contig(BloomFilter &bf, KmerMapper &mapper, Kmer km);
MakeContig make_contig(BloomFilter &bf, KmerMapper &mapper, Kmer km);

#endif // BFG_CONTIGMETHODS_HPP
