#include <cstdio>
#include <ctime>
#include <iostream>

#include "../BlockedBloomFilter.hpp"

using namespace std;

int main(int argc, char *argv[]) {
	if (argc < 2) {
		cout << "usage: BlockedBloomFilter <bits>\n where bits is the bit count per variable" << endl;
		return 1;
	}
	unsigned int bits = atoi(argv[1]), limit = 10000000, wrong = 0, counter = 0;

	// take bits as parameter
	BlockedBloomFilter BF(limit, (size_t) bits, (uint32_t) time(NULL));


	// insert ints 0,...,1M
	for (int j=0; j<=limit/10;j++)
		BF.insert(j);

	// assert 0,...,1M
	for (int j=0; j<=limit/10;j++)
		assert(BF.contains(j));
	
	for (int j=1+limit/10; j<=limit ;j++) {
		if (BF.contains(j))
			wrong++;
		counter++;
	}

	// compute false positive rate
	printf("False positive ratio: %.6f\n", wrong / (0.0 + counter));
  cout << &argv[0][2] << " completed successfully" << endl;
}
