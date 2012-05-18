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
		cout << " required string in argument " << endl;
        return 1;
    }

    // Read k from argument and set kmer size
    unsigned int k = strlen(argv[1]); 
	char *s = argv[1];
	std::cout << s << std::endl;
	
	char *tmp1 = (char *) malloc(k * sizeof(char));;
	char *tmp2 = (char *) malloc(k * sizeof(char));;


    Kmer::set_k(k);

    //char letters[] = {'A', 'C', 'G', 'T'};
    //map<int, int> baseKey = {{'A',0}, {'C',1}, {'G',2}, {'T',3}};

    Kmer K1, K2;
	//CompressedSequence C1, C2;

	K1 = Kmer(s);
	K1.toString(tmp1);
	
	//C1 = CompressedSequence(s);
	//C1.toString(tmp2);

	cout << "From Kmer:" << tmp1 << endl;
	cout << "From CompressedSequence:" << tmp2 << endl;
    return 0;
}
