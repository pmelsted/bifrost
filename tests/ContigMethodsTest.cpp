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
  //Kmer km1("ATATATACTATATAGTATATATACTATATAGTATATATAC");
  char s[] = "ATATATACTATATAGTATATATACTATATAGTATATATAC";
  int len = strlen(s);
  for(int i=0; i+31 <= len; ++i) {
    Kmer km(s + i);
    MakeContig mc = make_contig(BF, m, km);
    printf("km: %s contig: %s\n", km.toString().c_str(), mc.seq.c_str());
    /*
    FindContig fc1 =  find_contig_forward(BF, km);
    printf("sl: %d km: %s seq: %s\n", fc1.selfloop, km.toString().c_str(), fc1.s.c_str());
    if (fc1.selfloop == 2) {
        Kmer end =  find_contig_forward(BF, km.twin()).end.twin();
        printf("end=%s\n", end.toString().c_str());
        FindContig fc2 =  find_contig_forward(BF, end);
        printf("Contig after end trix seq: %s\n", fc2.s.c_str());
    }
    */
  }
  //string seq = "GTATATATACTATATAGTATATATACTATATAGTATATATA";


  
  cout << &argv[0][2] << " completed successfully" << endl;
}
