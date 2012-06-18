#include <cstdio>
#include <ctime>
#include <iostream>
#include <map>

#include "../KmerMapper.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  int k = 4, n = 8; 
  srand(time(NULL));

  Kmer::set_k(k);

  KmerMapper mapper1, mapper2;

  
  char s1[] = "ACGGTTT";
  char s2[] = "TTTCCCC";

  Kmer km1(s1+3), km2(s2+3);
  mapper1.addContig(s1);
  mapper1.addContig(s2);
  ContigRef cr1 = mapper1.getContig(0);
  cr1.ref.contig->cov[0] += 1;
  ContigRef cr2 = mapper1.getContig(1);
  cr2.ref.contig->cov[0] += 2;
  assert(cr1.ref.contig->seq.toString() == (string) s1);
  assert(cr2.ref.contig->seq.toString() == (string) s2);
  cr1 = mapper1.find(km1);
  cr2 = mapper1.find(km2);

  ContigRef joined = mapper1.joinContigs(cr1, cr2);
  Contig newc = *(mapper1.getContig(joined).ref.contig);
  assert(newc.seq.toString() == "ACGGTTTCCCC");
  assert(newc.covlength == 8);
  assert(newc.cov[0] == 1);
  assert(newc.cov[1] == 0);
  assert(newc.cov[4] == 2);

  
  char s3[] = "AAGGCCC";
  char s4[] = "ATATGGG";
  mapper2.addContig(s3);
  mapper2.addContig(s4);

  
  cr1 = mapper2.getContig(0);
  cr1.ref.contig->cov[1] += 5;
  cr2 = mapper2.getContig(1);
  cr2.ref.contig->cov[3] += 3;
  assert(cr1.ref.contig->seq.toString() == (string) s3);
  assert(cr2.ref.contig->seq.toString() == (string) s4);
  km1 = Kmer(s3+3);
  km2 = Kmer(s4+3);
  cr1 = mapper2.find(km1);
  cr2 = mapper2.find(km2);
  joined = mapper2.joinContigs(cr1, cr2);
  newc = *(mapper2.getContig(joined).ref.contig);
  
  assert(newc.cov[0] == 0);
  assert(newc.cov[1] == 5);
  assert(newc.cov[2] == 0);
  assert(newc.cov[4] == 3);
  assert(newc.covlength == 8);


  assert(mapper2.getContig(joined).ref.contig->seq.toString() == "AAGGCCCATAT");
  cout << &argv[0][2] << " completed successfully" << endl;
}
