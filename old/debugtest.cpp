#include <iostream>
#include "ContigMapper.hpp"
#include <cstdio>


int main() {
  int k = 31;
  Kmer::set_k(k);
  const char* s = "AGCTGGTTGCGCGTCTCCTGCGGCAGCTGCGCAACAGCTTCAAAGTAGTAGCAAATCTGCGCCAGCAAACGGCTGATGTTAATCGAGTTAGCCGAGTTTAACCCTAGCGCCACTTTCAGTTCTTCATCATCAAACGCCTGCTTCACCAGCGCCTGACAGGCATCGAAATCGCCGTCGATGGCAACAGTTTCGATATTGCCGCCCAATGTACAGAACAGTTTTTCTTGCAGTGGACTGATTTTGCCTCGTGGATAGAGGATAACCACTTTCACATTCGGTAAACCGTAGAAAGCATGAGCCACTGCCGCTCCGGTATCACCGGAGGTCGCGGTCAGAATGGTCACTGGCTTATCACCCGCAATATGGGTCAGCATTTGTGCCATAAAGCGACCGCCGAAATCTTTAAATGCCAGCGTTGGCCCGTGGAACAATTCCA";

  const char* read = s;

  ContigMapper cm;
  BloomFilter bf;
  FILE* f = fopen("../test/tiny/test.bf","r");
  bf.ReadBloomFilter(f);
  fclose(f);

  bf.count();
  cm.mapBloomFilter(&bf);
  size_t pos = 50;
  Kmer km(read + pos);
  cout << bool(bf.contains(km)) << endl;

  ContigMap cc = cm.findContig(km,read,pos);
  
  if (cc.isEmpty) {
    string s;
    cm.findContigSequence(km,s);
    //cout << s << endl;
    cm.addContig(km, read, pos);
  }
  cout << "find contig" << endl;

  cc = cm.findContig(km,read,pos);
  cout << cc.isEmpty <<  cc.dist << " " << cc.len << " " << cc.size << endl;;
  return 0;
}
