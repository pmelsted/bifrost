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

// use:   c.cover(start,end)
// pre:   0 <= start <= end < c.numKmers();
// post:  contig is covered from start to end and coveragesum is updated, threadsafe
void Contig::cover(size_t start, size_t end) {
  assert(0 <= start);
  assert(start <= end);
  ccov.cover(start,end);
  __sync_add_and_fetch(&coveragesum,end-start+1);
}
