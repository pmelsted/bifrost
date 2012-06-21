#include <cstdio>
#include <ctime>
#include <iostream>
#include <map>

#define private public 
#include "../KmerMapper.hpp"
#undef private

using namespace std;

int main(int argc, char *argv[]) {
  int k = 4, n = 8; 
  srand(time(NULL));

  Kmer::set_k(k);

  KmerMapper mapper1, mapper2, mapper3, mapper4, mapper5, mapper6, mapper7;

  
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


  char s5[] = "AAAATCCCC";
  mapper3.addContig(s5);
  cr1 = mapper3.getContig(0);
  Contig *cn = cr1.ref.contig;
  for (size_t j=0; j<6; ++j) {
    cn->cov[j] = 2;
  }
  cn->cov[3] = 1;
  mapper3.splitAndJoinContigs();
  assert(mapper3.contigCount() == 3);
  assert(mapper3.contigs[0].isEmpty());
  assert(!mapper3.contigs[1].isEmpty());
  assert(!mapper3.contigs[2].isEmpty());
  assert(mapper3.contigs[1].ref.contig->seq.toString() == "AAAATC");
  assert(mapper3.contigs[2].ref.contig->seq.toString() == "TCCCC");
  for (size_t j=0; j<=2; ++j) {
    assert(mapper3.contigs[1].ref.contig->cov[j] == 2);
  }
  for (size_t j=0; j<=1; ++j) {
    assert(mapper3.contigs[2].ref.contig->cov[j] == 2);
  }
  
  char s6[] = "AAAATCCCC";
  mapper4.addContig(s6);
  cr1 = mapper4.getContig(0);
  cn = cr1.ref.contig;
  for (size_t j=0; j<=5; ++j) {
    cn->cov[j] = 1;
  }
  for (size_t j=2; j<=4; ++j) {
    cn->cov[j] = 2;
  }
  mapper4.splitAndJoinContigs();
  assert(mapper4.contigCount() == 2);
  assert(mapper4.contigs[0].isEmpty());
  assert(!mapper4.contigs[1].isEmpty());
  assert(mapper4.contigs[1].ref.contig->seq.toString() == "AATCCC");
  for (size_t j=0; j<=2; ++j) {
    assert(mapper4.contigs[1].ref.contig->cov[j] == 2);
  }
  
  char s7[] = "AAAATCCCC";
  mapper5.addContig(s7);
  cr1 = mapper5.getContig(0);
  cn = cr1.ref.contig;
  for (size_t j=0; j<=5; ++j) {
    cn->cov[j] = 2;
  }
  cn->cov[0] = 1;
  cn->cov[2] = 1;
  cn->cov[5] = 1;
  mapper5.splitAndJoinContigs();
  assert(mapper5.contigCount() == 3);
  assert(mapper5.contigs[0].isEmpty());
  assert(!mapper5.contigs[1].isEmpty());
  assert(!mapper5.contigs[2].isEmpty());
  assert(mapper5.contigs[1].ref.contig->seq.toString() == "AAAT");
  assert(mapper5.contigs[2].ref.contig->seq.toString() == "ATCCC");
  assert(mapper5.contigs[1].ref.contig->cov[0] == 2);
  for (size_t j=0; j<=1; ++j) {
    assert(mapper5.contigs[2].ref.contig->cov[j] == 2);
  }

  char s8[] = "ACGAAAG";
  char s9[] = "AAGCTTA";
  mapper6.addContig(s8);
  mapper6.addContig(s9);
  cr1 = mapper6.getContig(0);
  cr2 = mapper6.getContig(1);
  for (size_t j=0; j<=3; ++j) {
    cr1.ref.contig->cov[j] = 2;
    cr2.ref.contig->cov[j] = 2;
  }
  mapper6.splitAndJoinContigs();
  assert(mapper6.contigs.size() == 3);

  assert(mapper6.contigs[2].ref.contig->seq.toString() == "ACGAAAGCTTA");
  Kmer k1(s8);
  ContigRef _c1 = mapper6.find(k1);
  assert(_c1.ref.idpos.id == 2);
  assert(_c1.ref.idpos.pos == 0);
  Kmer k2(s9);
  ContigRef _c2 = mapper6.find(k2);
  cout << _c2.ref.idpos.pos << endl;
  assert(_c2.ref.idpos.id == 2);
  assert(_c2.ref.idpos.pos == 4);
  

  char s10[] = "AGTCAGTTAAC";
  char s11[] = "AACGTAGG";
  mapper7.addContig(s10);
  mapper7.addContig(s11);
  cr1 = mapper7.getContig(0);
  cr2 = mapper7.getContig(1);
  for (size_t j=0; j<=10; ++j) {
    cr1.ref.contig->cov[j] = 2;
  }
  for (size_t j=0; j<=7; ++j) {
    cr2.ref.contig->cov[j] = 2;
  }
  mapper7.splitAndJoinContigs();
  assert(mapper7.contigs.size() == 3);

  assert(mapper7.contigs[2].ref.contig->seq.toString() == "AGTCAGTTAACGTAGG");

  

  
  
  cout << &argv[0][2] << " completed successfully" << endl;
}
