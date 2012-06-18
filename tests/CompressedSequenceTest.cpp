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

  // Read n from argument
  unsigned int n = atoi(argv[1]), i; 
  unsigned int LIM = 1 << n;

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

  Kmer::set_k(4);
  Kmer km = C3.getKmer(4);
  char *kmrstr = new char[5];
  km.toString(kmrstr);
  assert((string) kmrstr == (string) "CGGG");


  string T = "TTTT";
  assert(C3.jump(T.c_str(), 0, 0, false) == 4);
  T = "TTCG";
  assert(C3.jump(T.c_str(), 0, 2, false) == 4);
  T = "ACCCG";
  assert(C3.jump(T.c_str(), 1, 7, true) == 4);
  T = string(10, 'A');
  T.append("GAAA");
  assert(C3.jump(T.c_str(), 10, 4, true) == 4);
  T = "TT";
  assert(C3.jump(T.c_str(), 0, 0, false) == 2);
  T = "CCC";
  assert(C3.jump(T.c_str(), 0, 7, true) == 3);
  T = "AAAA";
  assert(C3.jump(T.c_str(), 0, 3, true) == 4);
  T = "NNNN";

  for (size_t m=0; m<8; m++) {
    for (size_t h=0; m<4; m++) {
      assert(C3.jump(T.c_str(), h, m, true) == 0);
      assert(C3.jump(T.c_str(), h, m, false) == 0);
    }
  }


  CompressedSequence C4("ACGT");
  

  cout << &argv[0][2] << " completed successfully" << endl;
  return 0;
}
