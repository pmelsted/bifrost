#include <cstdio>
#include <ctime>

#include "../BlockedBloomFilter.hpp"

int main(void) {
	//BF.insert();
	BlockedBloomFilter BF(10, (size_t) 4, (uint32_t) time(NULL));
	BF.insert(30);
	for (int j=0; j<100;j++)
		cout << BF.contains(j) <<  endl;
}
