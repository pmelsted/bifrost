#include <cstdio>
#include <ctime>
#include <iostream>
#include <map>

#include "../KmerMapper.hpp"

using namespace std;

int main(int argc, char *argv[]) {
	if (argc < 2) {
		cout << "usage: KmerMapperTest <k> <n>\n where k is kmer size and n is contig size" << endl;
		return 1;
	}
	unsigned int k = atoi(argv[1]), n = atoi(argv[2]);
  srand(time(NULL));

	// take bits as parameter
  Kmer::set_k(k);
  KmerMapper mapper;
  char letters[] = {'A', 'C', 'G', 'T'};

  char *s = new char[n+1];
  for(int i=0;i<n;i++) 
    s[i] = letters[rand() & 3];
  s[n] = '\0';
  cout << s << endl;

  Kmer km(s), km2;
  mapper.addContig(s);
  map<Kmer, ContigRef> vf;
  ContigRef a = vf[km2];
  assert(a.isEmpty());

}
