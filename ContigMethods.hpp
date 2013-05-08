#ifndef BFG_CONTIGMETHODS_HPP
#define BFG_CONTIGMETHODS_HPP

#include <string>
#include "Common.hpp"
#include "Kmer.hpp"
#include "KmerMapper.hpp"
#include "BloomFilter.hpp"
#include "FindContig.hpp"

/* Structs for Contig and Kmer information */

/*
struct NewContig {
  string seq;
  size_t start;
  size_t end;
  size_t read_index;
  int selfloop; // 0 for no selfloop, 1 for regular, 2 for reversed
  NewContig(string s, size_t i, size_t j, size_t ri, int sl) 
       : seq(s), start(i), end(j), read_index(ri), selfloop(sl) {}
       };*/

/*

struct CheckContig {
  ContigRef cr;
  size_t dist;
  bool repequal;
  CheckContig(ContigRef ref, size_t i, bool eq) : cr(ref), dist(i), repequal(eq) {}
};
*/

struct ContigMap {
  /**
     xxxxxxxxxxxxxxxxxxxxxxxxxxxxx contig
    |dist ->         xxxxxxxxxyyyyyyy  read on forward strand                
                     | len ->|
          yyyyXXXXXXXX                 read on reverse strand
    | dist -> | len  | 

   */
  Kmer head; // reference to start of contig
  size_t dist; // 0-based distance from start of contig
  size_t len;  // length of match, >= 1
  size_t size; // length of the contig in k-mers
  bool isEmpty; // true if proper match found
  bool isShort; // true if the contig is short
  bool strand; // true for forward strand
  ContigMap(Kmer ref, size_t i, size_t l,  size_t sz, bool eq, bool sh) : dist(i), strand(eq), size(sz), len(l), isShort(sh), head(ref), isEmpty(false){}
  ContigMap(size_t l = 1) : isEmpty(true), len(l) {}
};


struct NewContig {
  Kmer km;
  string read;
  size_t pos;
  NewContig(Kmer o, string &s, size_t p) : km(o), read(s), pos(p) {}
};

/*
struct MakeContig {
  string seq;
  int selfloop; // 0 for no selfloop, 1 for regular, 2 for reversed
  size_t pos;
  bool empty;
  MakeContig(string s, int l, size_t p, bool e) : seq(s), selfloop(l), pos(p), empty(e) {}
  }; */

/* Methods for Contig and Kmer information and mapping */

// void getMappingInfo(const bool repequal, const int32_t pos, const size_t dist, size_t &kmernum, int32_t &cmppos);
// CheckContig check_contig(BloomFilter &bf, KmerMapper &mapper, Kmer km);
// MakeContig make_contig(BloomFilter &bf, KmerMapper &mapper, Kmer km);

#endif // BFG_CONTIGMETHODS_HPP
