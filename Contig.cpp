#include "Contig.hpp"


// use:  c = allocateCov(full);
// pre:  cov is a NULL pointer and seq is not NULL
// post: cov has space to store coverage for all kmers in this contig 
void Contig::initializeCoverage(bool full) {
  size_t ssz = seq.size(), k = Kmer::k;
  assert(ssz >= k);
  coveragesum = 0;
  ccov.initialize(ssz - k + 1, full);
} 

// use:   c.cover(start,end)
// pre:   0 <= start , end < c.numKmers();
// post:  contig is covered from start to end (inclusive) and coveragesum is updated, threadsafe
void Contig::cover(size_t start, size_t end) {
  ccov.cover(start,end);
  if( end < start) {
    swap(start,end);
  }
  __sync_add_and_fetch(&coveragesum,end-start+1);
}


// use:  i = c.memory();
// post: i is the total memory used by c in bytes 
size_t Contig::memory() const {
  size_t m = sizeof(ccov) + sizeof(seq);
  size_t numkmers = numKmers();
  if (numkmers > ccov.size_limit) {
    m += ((numkmers + 3) / 4) + 8;
  }
  size_t seqlength = length();
  if (!seq.isShort()) {
    m += ((seqlength + 3) / 4);
  }
  return m;
}
