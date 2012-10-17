#include <cstdio>
#include <ctime>
#include <iostream>
#include <map>

#define private public 
#include "../KmerMapper.hpp"
#undef private
#include "../ContigMethods.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  int k = 4, n = 8; 
  srand(time(NULL));

  Kmer::set_k(k);
  Kmer km_del;
  km_del.set_deleted();

  KmerMapper mapper1, mapper2, mapper3, mapper4, mapper5, mapper6, mapper7;

  
  char s1[] = "ACGGTTT";
  char s2[] = "TTTCCCC";

  Kmer km1(s1+3), km2(s2+3);
  mapper1.addContig(s1);
  mapper1.addContig(s2);
  ContigRef cr1 = mapper1.getContig(0);
  cr1.ref.contig->cover(0,3);
  cr1.ref.contig->cover(0,3);
  ContigRef cr2 = mapper1.getContig(1);
  cr2.ref.contig->cover(0,3);
  cr2.ref.contig->cover(0,3);
  assert(cr1.ref.contig->seq.toString() == (string) s1);
  assert(cr2.ref.contig->seq.toString() == (string) s2);
  cr1 = mapper1.find(km1);
  cr2 = mapper1.find(km2);

  int id = mapper1.joinContigs(cr1, cr2, 1, 1);
  assert(id == 2);
  Contig newc = *(mapper1.getContig(2).ref.contig);
  assert(newc.seq.toString() == "ACGGTTTCCCC");
  assert(newc.numKmers() == 8);

  char s3[] = "AAGGCCC";
  char s4[] = "ATATGGG";
  mapper2.addContig(s3);
  mapper2.addContig(s4);

  
  cr1 = mapper2.getContig(0);
  cr1.ref.contig->cover(0,3);
  cr1.ref.contig->cover(0,3);
  cr2 = mapper2.getContig(1);
  cr2.ref.contig->cover(0,3);
  cr2.ref.contig->cover(0,3);
  assert(cr1.ref.contig->seq.toString() == (string) s3);
  assert(cr2.ref.contig->seq.toString() == (string) s4);
  km1 = Kmer(s3+3);
  km2 = Kmer(s4+3);
  cr1 = mapper2.find(km1);
  cr2 = mapper2.find(km2);
  id = mapper2.joinContigs(cr1, cr2, 1, -1);
  assert(id == 2);
  newc = *(mapper2.getContig(5).ref.contig);
  
  assert(newc.numKmers() == 8);
  assert(mapper2.getContig(id).ref.contig->seq.toString() == "AAGGCCCATAT");

  
  /* Test the deletion of kmers */
  KmerMapper dMapper;
  char dContig[] = "ACGGTTTCCCC";
  dMapper.addContig(dContig);
  dMapper.map.set_deleted_key(km_del);
  Kmer dFirst(dContig);
  assert(dMapper.stride == 4);
  int numkmers = strlen(dContig) - 4 + 1;
  for(int i=0; i < (numkmers - 1); ++i){
    Kmer km(dContig+i);
    Kmer rep = km.rep();
    if ((i % 4) == 0) {
      assert(!dMapper.find(km).isEmpty());
      dMapper.map.erase(km.rep());
    }
    assert(dMapper.find(km).isEmpty());
  }
  Kmer km(dContig + (numkmers -1));
  Kmer rep = km.rep();
  assert(!dMapper.find(km).isEmpty());
  dMapper.map.erase(km.rep());
  assert(dMapper.find(km).isEmpty());

  /* Test the splitContigs method */
  KmerMapper splitMapper;
  char splitContigString[] = "ACACTAGAGTAAAA";
  splitMapper.addContig(splitContigString);
  ContigRef splitRef = splitMapper.getContig(0);
  Contig *splitContig = splitRef.ref.contig;
  splitContig->cover(0,10);
  splitContig->cover(0,0);
  splitContig->cover(5,5);
  splitContig->cover(10,10);
  
  splitMapper.map.set_deleted_key(km_del);
  vector<pair<int, int> > spv = splitContig->ccov.splittingVector();
  assert(spv[0].first == 0); 
  assert(spv[0].second == 1); 
  assert(spv[1].first == 5); 
  assert(spv[1].second == 6); 
  assert(spv[2].first == 10); 
  assert(spv[2].second == 11); 
  pair<size_t, size_t> splitpair = splitMapper.splitContigs();
  assert(splitpair.first == 2);
  assert(splitpair.second == 0);
  assert(splitMapper.getContig(1).ref.contig->seq.toString() == "ACAC");
  assert(splitMapper.getContig(2).ref.contig->seq.toString() == "AGAG");
  assert(splitMapper.getContig(3).ref.contig->seq.toString() == "AAAA");

  
  /* Test the splitAndJoinContigs method */  
  char s5[] = "AAAATCCCC";
  mapper3.addContig(s5);
  cr1 = mapper3.getContig(0);
  Contig *cn = cr1.ref.contig;
  cn->cover(0,5);
  cn->cover(0,2);
  cn->cover(4,5);
  mapper3.splitAndJoinContigs();
  assert(mapper3.contigCount() == 3);
  assert(mapper3.contigs[0].isEmpty());
  assert(!mapper3.contigs[1].isEmpty());
  assert(!mapper3.contigs[2].isEmpty());
  assert(mapper3.contigs[1].ref.contig->seq.toString() == "AAAATC");
  assert(mapper3.contigs[2].ref.contig->seq.toString() == "TCCCC");
  assert(mapper3.contigs[1].ref.contig->coveragesum >= 6);
  assert(mapper3.contigs[2].ref.contig->coveragesum >= 4);
  
  char s6[] = "AAAATCCCC";
  mapper4.addContig(s6);
  cr1 = mapper4.getContig(0);
  cn = cr1.ref.contig;
  cn->cover(0,5);
  cn->cover(2,4);
  mapper4.splitAndJoinContigs();
  assert(mapper4.contigCount() == 2);
  assert(mapper4.contigs[0].isEmpty());
  assert(!mapper4.contigs[1].isEmpty());
  assert(mapper4.contigs[1].ref.contig->seq.toString() == "AATCCC");
  
  assert(mapper3.contigs[1].ref.contig->coveragesum >= 6);
  
  char s7[] = "AAAATCCCC";
  mapper5.addContig(s7);
  cr1 = mapper5.getContig(0);
  cn = cr1.ref.contig;
  cn->cover(0,5);
  cn->cover(1,1);
  cn->cover(3,4);
  mapper5.splitAndJoinContigs();
  assert(mapper5.contigCount() == 3);
  assert(mapper5.contigs[0].isEmpty());
  assert(!mapper5.contigs[1].isEmpty());
  assert(!mapper5.contigs[2].isEmpty());
  assert(mapper5.contigs[1].ref.contig->seq.toString() == "AAAT");
  assert(mapper5.contigs[2].ref.contig->seq.toString() == "ATCCC");

  char s8[] = "ACGAAAG";
  char s9[] = "AAGCTTA";
  mapper6.addContig(s8);
  mapper6.addContig(s9);
  
  cr1 = mapper6.getContig(0);
  cr1.ref.contig->cover(0,3);
  cr1.ref.contig->cover(0,3);

  cr2 = mapper6.getContig(1);
  cr2.ref.contig->cover(0,3);
  cr2.ref.contig->cover(0,3);
  
  mapper6.splitAndJoinContigs();
  assert(mapper6.contigs.size() == 2 ); // Two candidates for AAGX


  KmerMapper revself;


  
  cout << &argv[0][2] << " completed successfully" << endl;
}
