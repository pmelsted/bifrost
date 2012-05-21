#include <cstdio>
#include <ctime>

#include "../BlockedBloomFilter.hpp"

int main(void) {
	//BF.insert();
  // take bits as parameter
	BlockedBloomFilter BF(10, (size_t) 4, (uint32_t) time(NULL));
	BF.insert(30);
	for (int j=0; j<100;j++)
		cout << BF.contains(j) <<  endl;


	// insert ints 0,...,1M

	// assert 0,...,1M

	// check 1M+1,...,10M

	// compute false positive rate
}
