#include <vector>
#include "FindContig.hpp"

// use:  r = isNeighbor(a,b)
// pre:
// post: r is true if a[1:k-1]+c == b for some c
bool isNeighbor(Kmer a, Kmer b) {
  for (size_t i = 0; i < 4; ++i) {
    if (b == a.forwardBase(alpha[i])) {
      return true;
    }
  }
  return false;
}

// use:  fc = find_contig_forward(bf,km,s);
// pre:  
// post: km is contained in a contig c with respect to the
//       bloom filter graph bf and fc.end is the forward endpoint (wrt km direction)
//       (if we see the original kmer again we stop to prevent infinite loop)
//       and c contains fc.dist kmers until the end (including km)
//       the sequence of the contig is stored in fc.s
FindContig find_contig_forward(BloomFilter &bf, Kmer km) {
  assert(bf.contains(km.rep()));

  int j;
  bool selfloop = false;
  size_t i, dist = 1;
  
  Kmer first = km, end = km;
  string s = km.toString();
  
  while (true) {
    assert(bf.contains(end.rep()));
    size_t fw_count = 0;
    j = -1;
    for (i = 0; i < 4; ++i) {
      Kmer fw_rep = end.forwardBase(alpha[i]).rep();
      if (bf.contains(fw_rep)) {
        j = i;
        ++fw_count;
        if (fw_count > 1) {
          break;
        }
      }
    }

    if (fw_count != 1) {
      break;
    }
    
    Kmer fw = end.forwardBase(alpha[j]);
    assert(0 <= j && j < 4);
    assert(bf.contains(fw.rep()));
    if (first == fw) {
      selfloop = true;
      break;
    }

    size_t bw_count = 0;
    for (i = 0; i < 4; ++i) {
      Kmer bw_rep = fw.backwardBase(alpha[i]).rep();
      if (bf.contains(bw_rep)) {
        ++bw_count;
        if (bw_count > 1) {
          break;
        }
      }
    }

    assert(bw_count >= 1);
    if (bw_count != 1) {
      break;
    }
    
    end = fw;
    ++dist;
    s += alpha[j];
  }
  return FindContig(end, dist, selfloop, s);
}
