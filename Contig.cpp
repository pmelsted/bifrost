#include "Contig.hpp"

Contig::~Contig() {
  if (cov != NULL) {
    delete[] cov;
    cov = NULL;
  }
}
