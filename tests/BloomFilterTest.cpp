#include <cstdio>
#include <ctime>
#include <iostream>

#include "../BloomFilter.hpp"
#include "../Kmer.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "usage: BloomFilterTest <bits>\n where bits is the bit count per variable" << endl;
    return 1;
  }
  unsigned int bits = atoi(argv[1]), wrong = 0, counter = 0;
  cout << "bits=" << bits << endl;
  int limit = 1000000;

  // take bits as parameter
  BloomFilter BF (limit, (size_t) bits, (uint32_t) time(NULL));



  // insert ints 0,...,1M
  for (int j=0; j<limit;j++)
    BF.insert(j);

  // assert 0,...,1M
  for (int j=0; j<limit;j++)
    assert(BF.contains(j));
  
  for (int j=limit; j < 2*limit ;j++) {
    if (BF.contains(j))
      wrong++;
    counter++;
  }

  Kmer::set_k(31);
  Kmer km("ACGTACGTACGTACGTACGTACGTACGTACG");
  printf("sizeof(Kmer) == %d\n", (int) sizeof(Kmer));
  printf("sizeof(km) == %d\n", (int) sizeof(km));
  printf("pointer cast == %p\n", (const void*) &km);

  cout << "Count == " << BF.count() << endl;
  cout << "Contains(km) == " << BF.contains(km) << endl;
  cout << "Insert(km) == " << BF.insert(km) << endl;
  BF.count();
  cout << "Contains(km) == " << BF.contains(km) << endl;
  assert(BF.contains(km));
  
  BF.count();
  Kmer km2("ACGTACGTACGTACGTACGTACGTACGTAGG"); // AGG vs ACG in the end
  cout << "Contains(km2) == " << BF.contains(km2) << endl;
  cout << "Insert(km2) == " << BF.insert(km2) << endl;
  cout << "Contains(km2) == " << BF.contains(km2) << endl;
  assert(BF.contains(km2));

  BF.count();
  
  FILE *fp = fopen("testBloom.bf", "wb");
  BF.WriteBloomFilter(fp);
  fclose(fp);



  // compute false positive rate
  printf("False positive ratio: %.6f\n", wrong / (0.0 + counter));
  cout << &argv[0][2] << " completed successfully" << endl;
}
