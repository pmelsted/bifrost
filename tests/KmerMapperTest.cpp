#include <cstdio>
#include <ctime>
#include <iostream>
#include <map>

#include "../KmerMapper.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  /*
	if (argc < 3) {
		cout << "usage: KmerMapperTest <k> <n>\n where k is kmer size and n is contig size" << endl;
		return 1;
	}
	unsigned int k = atoi(argv[1]), n = atoi(argv[2]);

  char letters[] = {'A', 'C', 'G', 'T'};
  char *s = new char[n+1];
  for(int i=0;i<n;i++) 
    s[i] = letters[rand() & 3];
  s[n] = '\0';
  */
  int k = 4, n = 8; 
  srand(time(NULL));

  Kmer::set_k(k);
  KmerMapper mapper1, mapper2;

  
  char s1[] = "ACGGTTT";
  char s2[] = "TTTCCCC";

  Kmer km1(s1+4), km2(s2+4);
  mapper1.addContig(s1);
  mapper1.addContig(s2);
  ContigRef cr1 = mapper1.getContig(0);
  ContigRef cr2 = mapper1.getContig(1);
  assert(cr1.ref.contig->seq.toString() == (string) s1);
  assert(cr2.ref.contig->seq.toString() == (string) s2);
  cr1 = mapper1.find(km1);
  cr2 = mapper1.find(km2);
  ContigRef joined = mapper1.joinContigs(cr1,cr2);
  assert(mapper1.getContig(joined).ref.contig->seq.toString() == "ACGGTTTCCCC");

  
  char s3[] = "GGGGAAA";

  mapper2.addContig(s2);
  mapper2.addContig(s3);

  cr1 = mapper2.getContig(0);
  cr2 = mapper2.getContig(1);
  joined = mapper1.joinContigs(cr1, cr2);

  assert(mapper2.getContig(joined).ref.contig->seq.toString() == "ACGGTTTCCCC");

  cout << &argv[0][2] << " completed successfully" << endl;
}
