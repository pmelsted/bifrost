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
	
	//char *tmp1 = (char *) malloc(k * sizeof(char));;
	char *tmp2 = (char *) malloc(k * sizeof(char));;
	//char *tmp3 = (char *) malloc(k * sizeof(char));;
	char *tmp3 = new char[k];



    //char letters[] = {'A', 'C', 'G', 'T'};

	CompressedSequence C1, C2;
	C1 = CompressedSequence(s);
	C1.toString(tmp2);
	C2 = C1.rev();
	C2.toString(tmp3);

	cout << "From CompressedSequence: " << tmp2 << endl;
	
	cout << "From rev: " << tmp3 << endl;

    return 0;
}
