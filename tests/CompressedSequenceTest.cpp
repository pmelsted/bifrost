#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <ctime>
#include <map>

using namespace std;

#include "../Kmer.hpp"
#include "../CompressedSequence.hpp"


int main(int argc, char *argv[]) {
  if (argc < 2) {
  cout << "usage: CompressedSequenceTest <n>\n and a random sequence of length 2^n will be tested" << endl;
      return 1;
  }
    
	srand(time(NULL));

    // Read k from argument and set kmer size
  unsigned int k = atoi(argv[1]), i; 
	unsigned int LIM = 1 << k;

	char *s = new char[LIM];
	char *out = new char[LIM];
	char *rev = new char[LIM];


  char letters[] = {'A', 'C', 'G', 'T'};
  map<int, int> baseKey;
	baseKey['A'] = 0; baseKey['C'] = 1; baseKey['G'] = 2; baseKey['T'] = 3;
	
	for(i=0;i<LIM;i++)
		s[i] = letters[rand() & 3];

	CompressedSequence C1, C2;
	C1 = CompressedSequence(s);
	C1.toString(out);
	assert(strcmp(s, out) == 0);
	C2 = C1.rev();
	C2.toString(rev);
	for(i=0;i<LIM;i++) {
		assert(rev[LIM-i-1] == letters[3-baseKey[s[i]]]);
	}

  string m = "CCCG";
  string m2;
  CompressedSequence C3("TTTTTTTT");
  C3.setSequence(m.c_str(), 4, 4, false);
  assert(C3.toString() == (string) "TTTTCCCG");
  C3.setSequence(m.c_str(), 4, 4, true);
  assert(C3.toString() == (string) "TTTTCGGG");

  cout << &argv[0][2] << " completed successfully" << endl;
  return 0;
}
