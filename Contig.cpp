#include "Contig.hpp"


// use:  delete c;
// pre:  c is a pointer to a contig
// post: the memory location which c points to has been freed
Contig::~Contig() {
  if (cov != NULL) {
    delete[] cov;
    cov = NULL;
  }
}

// use:  c = allocateCov();
// pre:  cov is a NULL pointer and seq is not NULL
// post: cov is a array with zeros of length : seq.size() - Kmer::k +1
void Contig::allocateCov() {
  covlength = seq.size()-Kmer::k+1;
  cov = new uint8_t[covlength];
  memset(cov, 0, covlength);
}
