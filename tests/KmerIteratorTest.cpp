#include <cstdio>
#include <ctime>
#include <iostream>
#include <map>

#include "../KmerIterator.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  int k = 31;
  Kmer::set_k(k);

  string s = "GGAGAAGTGCGTAATCTGGCCCAGCGCAGCGCCCAGGCGGCTCGTGAANTTAAAAGCCTGATTGAAGACTCGGTGGGGGAAGTGGATGTTGGCTCTACGC";

  KmerIterator iter(s.c_str());

  Kmer km, rep;
  km = iter->first;
  char kmrstr[200];
  km.toString(kmrstr);
  printf("%u : %d : %s\n", 0, iter->second, kmrstr);
  for(size_t j=1; j < s.size() - k; j++) {
    iter.raise(km, rep);
    km.toString(kmrstr);
    printf("%u : %d : %s\n", j, iter->second, kmrstr);
  }
  
  cout << &argv[0][2] << " completed successfully" << endl;
}
