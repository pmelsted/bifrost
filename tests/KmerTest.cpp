#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "../Kmer.hpp"

using namespace std;

void twinString(char *from, char *to) {
    while (*from) {
        switch(from[0]) {
            case 'A': 
                to[0] = 'T';
                break;
            case 'C':
                to[0] = 'G';
                break;
            case 'G':
                to[0] = 'C';
                break;
            case 'T':
                to[0] = 'A';
                break;
        }
        from++;
        to++;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "usage: KmerTest <k>\n Where k is kmer size" << endl;
        return 0;
    }

    // Read k from argument and set kmer size
    unsigned int k = atoi(argv[1]), i;
    Kmer::set_k(k);

    // Allocate strings to use
    char *real = (char *) malloc((k + 1) * sizeof(char));
    char *fromkmer = (char *) malloc((k + 1) * sizeof(char));
    char *realtwin = (char *) malloc((k + 1) * sizeof(char));

    for (i=0; i<k; i++) 
        real[i] = 'A';
    real[k] = '\0';

    Kmer KM(real);

    // Check the kmer
    KM.toString(fromkmer);
    if (strcmp(real, fromkmer) != 0) {
        cout << "Was expecting: " << real << " from toString() but got: " << fromkmer << endl;
        return 0;
    }

    // Check the twin kmer
    Kmer TW = KM.twin();
    TW.toString(fromkmer);
    twinString(real, realtwin);
    if (strcmp(fromkmer, realtwin) != 0) {
        cout << "Was expecting: " << realtwin << " from toString() but got: " << fromkmer << endl;
        return 0;
    }

    cout << "All tests completed successfully" << endl;
    return 0;
}
