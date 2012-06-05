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
  KmerMapper mapper;

  
  char s1[] = "ACGGTTTT";
  char s2[] = "TTTTCCCC";

  Kmer km1(s1+4), km2(s2+4);
  mapper.addContig(s1);
  mapper.addContig(s2);
  ContigRef cr1 = mapper.getContig(0);
  ContigRef cr2 = mapper.getContig(1);
  assert(cr1.ref.contig->seq.toString() == (string)  s1);
  assert(cr2.ref.contig->seq.toString() == (string)  s2);
  cr1 = mapper.find(km1);
  cr2 = mapper.find(km2);
  ContigRef joined = mapper.joinContigs(cr1,cr2);
  assert(mapper.getContig(joined).ref.contig->seq.toString() == "ACGGTTTTCCCC");
  cout << &argv[0][2] << " completed successfully" << endl;
}
