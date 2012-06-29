#include "Contig.hpp"


// use:  delete c;
// pre:  c is a pointer to a contig
// post: the memory location which c points to has been freed
Contig::~Contig() {
  if (cov != NULL) {
    delete[] cov;
    delete[] covp;
    covp = NULL;
    cov = NULL;
  }
}


// use:  c = allocateCov();
// pre:  cov is a NULL pointer and seq is not NULL
// post: cov has space to store coverage for all kmers in this contig 
void Contig::allocateCov() {
  if (seq.size() >= Kmer::k) {
    covlength = seq.size() - Kmer::k + 1;
    covp = new CompressedCoverage(covlength);
    cov = new uint8_t[covlength];
    for(int i=0;i < covlength;++i) {
      cov[i] = 0;
    }
  }
} 
