#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../Kmer.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "usage: KmerTest <k>\n Where k is kmer size" << endl;
        return 0;
    }

    // Read k from argument
    unsigned int k = atoi(argv[1]), i;
    Kmer::set_k(k);

    // Allocate memory for the test string and initialize it
    char *s = (char *) malloc((k + 1) * sizeof(char));
    for(i=0; i<k; i++) 
        s[i] = 'A';
    s[k] = '\0';
    printf("Test string: %s\n", s);
    const char *test = "AAAAGG";

    // Make a kmer
    Kmer KM(s);
    //Kmer TW = KM.twin();

    // Allocate memory for a temporary string
    char *t = (char *) malloc((k + 1) * sizeof(char));
    char *tp = t;
    KM.toString(t);
    printf("KM.toString() gives: %s\n", tp);

    return 0;
}
