#include <cstdio>
#include <ctime>
#include <iostream>
#include <map>

#define private public 
#include "../KmerMapper.hpp"
#undef private
#include "../BloomFilter.hpp"
#include "../ContigMethods.hpp"

using namespace std;

int main(int argc, char **argv) {
  if (argc < 2) {
    cout << "Missing BloomFilter file argument" << endl;
    return 1;
  }
  int k = 31;

  Kmer::set_k(k);
  Kmer km_del;
  km_del.set_deleted();

  KmerMapper m;

  FILE *f = fopen(argv[1], "r");
  if (f == NULL) {
    cout << "BloomFilter file: " << argv[1] << " not found" << endl;
    return 1;
  }
  BloomFilter BF;
  BF.ReadBloomFilter(f);
  fclose(f);

  /* Custom tests */
  Kmer km1("ATATATACTATATAGTATATATACTATATAGTATATATAC");
  //string seq = "GTATATATACTATATAGTATATATACTATATAGTATATATA";
  MakeContig mc = make_contig(BF, m, km1);
  cout << mc.seq << endl;


  
  cout << &argv[0][2] << " completed successfully" << endl;
}
