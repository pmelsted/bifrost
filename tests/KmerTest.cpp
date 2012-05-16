#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>

#include "../Kmer.hpp"

using namespace std;

void twinString(char *from, char *to, unsigned int length) {
    while (*from) {
        switch(from[0]) {
            case 'A': 
                to[--length] = 'T';
                break;
            case 'C':
                to[--length] = 'G';
                break;
            case 'G':
                to[--length] = 'C';
                break;
            case 'T':
                to[--length] = 'A';
                break;
        }
        from++;
    }
}

void makeKmerString(char *s, unsigned int i) {
    while (*s) {
        switch(i % 4) {
            case 0: 
                s[0] = 'A';
                break;
            case 1:
                s[0] = 'C';
                break;
            case 2:
                s[0] = 'G';
                break;
            case 3:
                s[0] = 'T';
                break;
            printf("s[0]=%c\n", s[0]);
        }
        s++;
        i /= 4;
    }
}


int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "usage: KmerTest <k>\n Where k is kmer size" << endl;
        return 0;
    }

    // Read k from argument and set kmer size
    unsigned int k = atoi(argv[1]), index = 0, limit = (int) pow(4, k);
    int j;
    Kmer::set_k(k);

    // Allocate strings to use
    char *real = (char *) malloc((k + 1) * sizeof(char));
    char *fromkmer = (char *) malloc((k + 1) * sizeof(char));
    char *realtwin = (char *) malloc((k + 1) * sizeof(char));
    char letters[] = {'A', 'C', 'G', 'T'};
    for(index=0; index <k; index++)
        real[index] = 'A';
    real[k] = '\0';
    index = 0;
    Kmer K, TWIN, FW, BACK;
    while (index < limit) {
        makeKmerString(real, index++);

        K = Kmer(real);

        // Check the kmer
        K.toString(fromkmer);
        if (strcmp(real, fromkmer) != 0) {
            cout << "Was expecting the base string to be: " << real << " but got: " << fromkmer << endl;
            return 0;
        }

        // Check the twin kmer
        TWIN = K.twin();
        TWIN.toString(fromkmer);
        twinString(real, realtwin, k);
        if (strcmp(fromkmer, realtwin) != 0) {
            cout << "Was expecting twin to be: " << realtwin << " but got: " << fromkmer << endl;
            cout << "The base string was: " << real << endl;
            return 0;
        }
        
        for(j=0; j<4; j++) { 
            // Check the forward base
            FW = K.forwardBase(letters[j]);
            FW.toString(fromkmer);
            if ((strncmp(&real[1], fromkmer, k-1) != 0) || (letters[j] != fromkmer[k-1])) {
                printf("Was expecting the forward base to be: %s%c  but got: %s \n", &real[1], letters[j], fromkmer);
                return 0;
            }
        }
        
        for(j=0; j<4; j++) { 
            // Check the backward base
            BACK = K.backwardBase(letters[j]);
            BACK.toString(fromkmer);
            if ((strncmp(real, &fromkmer[1], k-1) != 0) || (letters[j] != fromkmer[0])) {
                real[k-1] = '\0';
                printf("Was expecting the backward base to be: %c%s  but got: %s \n", letters[j], real, fromkmer);
                return 0;
            }
        }
    }

    cout << "All tests completed successfully" << endl;
    return 0;
}
