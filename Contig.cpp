#include "Contig.hpp"


Contig::~Contig() {
}


// use:  c = allocateCov(full);
// pre:  cov is a NULL pointer and seq is not NULL
// post: cov has space to store coverage for all kmers in this contig 
void Contig::initializeCoverage(bool full) {
  size_t ssz = seq.size(), k = Kmer::k;
  assert(ssz >= k);
  coveragesum = ssz - k + 1;
  ccov.initialize(coveragesum, full);
} 
