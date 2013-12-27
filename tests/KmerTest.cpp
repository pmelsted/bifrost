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


// maps a number uniquely to a string
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


// tests all different kmers for the given kmer size
int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "usage: KmerTest <k>\n Where k is kmer size" << endl;
    return 1;
  }

  // Read k from argument and set kmer size
  unsigned int k = atoi(argv[1]), index = 0, limit = (unsigned int) ( 1 << (2*k));
  int j;
  Kmer::set_k(k);

  // Allocate strings to use
  char *real = (char *) malloc((k + 1) * sizeof(char));
  char *fromkmer = (char *) malloc((k + 1) * sizeof(char));
  char *realtwin = (char *) malloc((k + 1) * sizeof(char));
  char *tmp = (char *) malloc((k + 1) * sizeof(char));
  char *last = (char *) malloc((k + 1) * sizeof(char));
  char letters[] = {'A', 'C', 'G', 'T'};
  for(index=0; index <k; index++) {
    real[index] = 'A';
  }
  real[k] = '\0';
  strcpy(last,real);
  index = 0;
  Kmer K, Kp, TWIN, FW, BACK;
  while (index < limit) {
    makeKmerString(real, index);
    K = Kmer(real);

    if (true && k < 8) { // this has been checked for k = 1,..,7
      if (index > 0) {
        // Verify the operators
        for (unsigned int i2 = 0; i2 < index; i2++) {
          makeKmerString(last,i2);
          Kp = Kmer(last);
          
          if (Kp < K) {
            if (strcmp(last, real) >= 0) {
              cout << "Kmer with string: " << last << " is less than kmer with string: " << real << endl;
              return 1;
            }
          } else {
            if (strcmp(last, real) < 0) {
              cout << "Kmer with string: " << last << " is greater or equal than kmer with string: " << real << endl;
              return 1;
            }
          }
          
          
          if (Kp == K) {
            if (strcmp(last, real) != 0) {
              cout << "Kmer with string: " << last << " is equal to kmer with string: " << real << endl;
              return 1;
            }
          }
        }
      }
    }

    // Verify toString from this kmer
    K.toString(fromkmer);
    if (strcmp(real, fromkmer) != 0) {
      cout << "Was expecting the base string to be: " << real << " but got: " << fromkmer << endl;
      return 1;
    }

    // Verify toString from this kmer's twin
    TWIN = K.twin();
    TWIN.toString(fromkmer);
    twinString(real, realtwin, k);
    if (strcmp(fromkmer, realtwin) != 0) {
      cout << "Was expecting twin to be: " << realtwin << " but got: " << fromkmer << endl;
      cout << "The base string was: " << real << endl;
      return 1;
    }
    Kmer tw2(realtwin);
    if (TWIN != tw2) {
      cout << TWIN.toString() << endl;
      cout << TWIN.getBinary() << endl;
      cout << tw2.toString() << endl;
      cout << tw2.getBinary() << endl;
    }

    assert(TWIN == tw2);
    assert(memcmp(&TWIN,&tw2,sizeof(Kmer)) == 0);
    Kmer tw3 = TWIN.twin();
    assert(tw3 == K);
    assert(memcmp(&K,&tw3,sizeof(Kmer)) == 0);

      
    for(j=0; j<4; j++) { 
      // Verify toString from this kmer's forward bases
      FW = K.forwardBase(letters[j]);
      FW.toString(fromkmer);
      if ((strncmp(&real[1], fromkmer, k-1) != 0) || (letters[j] != fromkmer[k-1])) {
        printf("Was expecting the forward base to be: %s%c  but got: %s \n", &real[1], letters[j], fromkmer);
        return 1;
      }
      strncpy(tmp,real+1,k-1);
      tmp[k-1] = letters[j];
      Kmer fw2(tmp);
      if (fw2 != FW) {
        printf("index: %d, j: %d\n, tmp %s, real %s",index,j,tmp,real);
        return 1;
      }
      assert(memcmp(&fw2,&FW,sizeof(Kmer)) == 0);
      assert(FW.backwardBase(real[0]) == K);
    }
      
    for(j=0; j<4; j++) { 
      // Verify toString from this kmer's backward bases
      BACK = K.backwardBase(letters[j]);
      BACK.toString(fromkmer);
      if ((strncmp(real, &fromkmer[1], k-1) != 0) || (letters[j] != fromkmer[0])) {
        real[k-1] = '\0';
        printf("Was expecting the backward base to be: %c%s  but got: %s \n", letters[j], real, fromkmer);
        return 1;
      }
      strncpy(tmp+1,real,k-1);
      tmp[0] = letters[j];
      Kmer bw2(tmp);
      assert(bw2 == BACK);
      assert(memcmp(&bw2,&BACK,sizeof(Kmer)) == 0);
      assert(BACK.forwardBase(real[k-1]) == K);
    }
    index++;
  }

  cout << &argv[0][2] << " completed successfully" << endl;
  return 0;
}
