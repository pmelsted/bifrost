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
