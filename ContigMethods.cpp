#include "ContigMethods.hpp"

static const char beta[4] = {'T','G','A','C'}; // c -> beta[(c & 7) >> 1] maps: 'A' <-> 'T', 'C' <-> 'G'

// use:  getMappingInfo(repequal, pos, dist, kmernum, cmppos)
// pre:  Originally we have a kmer, call it km1. We go forward from this kmer dist times and get another kmer, call it km2.
//       repequal is true <==> km2 == km2.rep() 
//       km2.rep() maps to position pos in some contig, call it c.
// post: if (pos >= 0) == repequal: 
//            c.seq.getKmer(kmernum) == km1
//            kmernum + k == cmppos
//        else: 
//            c.seq.getKmer(kmernum) == km1.twin()
//            kmernum - 1 == cmppos
void getMappingInfo(const bool repequal, const int32_t pos, const size_t dist, size_t &kmernum, int32_t &cmppos) {
  size_t k = Kmer::k; 
  if (pos >= 0) {
    if (repequal) {
      cmppos = pos - dist + k;
      kmernum = cmppos - k;
    } else {
      cmppos = pos - 1 + dist;
      kmernum = cmppos +1;
    }
  } else {
    if (repequal) {
      cmppos = -pos + dist -k;
      kmernum = cmppos +1;
    } else {
      cmppos = -pos + 1 - dist; // Equivalent: (-pos +1 -k) - dist + k
      kmernum = cmppos - k;
    }
  }
}


// use:  cc = check_contig(bf, km, mapper);
// post: if km does not map to a contig: cc.cr.isEmpty() == true and cc.dist == 0
//       else: km is in a contig which cc.cr maps to and cc.dist is the distance 
//             from km to the mapping location 
//             (cc.eq == true):  km has the same direction as the contig
//             else:  km has the opposite direction to the contig
CheckContig check_contig(BlockedBloomFilter &bf, KmerMapper &mapper, Kmer km) {
  ContigRef cr = mapper.find(km);
  if (!cr.isEmpty()) {
    return CheckContig(cr, 0, km == km.rep());
  }
  int i, j;
  size_t dist = 1;
  bool found = false;
  Kmer end = km;
  while (dist < mapper.stride) {
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
    cr = mapper.find(fw);
    end = fw;
    if (!cr.isEmpty()) {
      found = true;
      break;
    }
    ++dist;
  }
  if (found) {
    return CheckContig(cr, dist, end == end.rep());
  } else {
    return CheckContig(ContigRef(), 0, 0);
  }
}

// use:  mc = make_contig(bf, mapper, km, s);
// pre:  km is not contained in a mapped contig in mapper 
// post: Finds the forward and backward limits of the contig
//       which contains km  according to the bloom filter bf and puts it into mc.seq
//       mc.pos is the position where km maps into this contig
MakeContig make_contig(BlockedBloomFilter &bf, KmerMapper &mapper, Kmer km) {
  /** There are several cases here:
   *
   * Case 0: Regular contig, no self-loops
   * Case 1: Self-looping contig:  firstkm -> ... -> lastkm -> firstkm -> ... ->lastkm
   * Case 2: Hairpinned contigs:
   *  a) firstkm -> ... -> lastkm -> twin(lastkm) -> ... -> twin(firstkm)
   *  b) twin(lastkm) -> ... -> twin(firstkm) -> firstkm -> ... -> lastkm
   *  c) firstkm -> ... -> lastkm -> twin(lastkm) -> ... -> twin(firstkm) -> firstkm -> ... -> lastkm -> ... (can repeat infinitely)
   *
   **/

  size_t k = Kmer::k;
  string seq;
  FindContig fc_fw = find_contig_forward(bf, km);
  int selfloop = fc_fw.selfloop;
  bool empty = false;

  if (selfloop == 0) {
    // Case 0, Case 1 or Case 2b 
    // No reverse self-loop on forward strand or two connections from km
  } else if (selfloop == 1) {
    // Case 1
    // We don't want to grow the contig backwards, it would duplicate kmers
    return MakeContig(fc_fw.s, selfloop, 0,false); 
  } else if (selfloop == 2) {
    // Case 2a or Case 2c
    // Reverse self-loop found on forward strand
    // Maybe we don't have all the contig yet, because km might not be equal to 
    // firstkm as described in Case 2a and Case 2c above
  } 

  FindContig fc_bw = find_contig_forward(bf, km.twin());
  ContigRef cr_tw_end = mapper.find(fc_bw.end);
  assert(cr_tw_end.isEmpty());

  // isolated contig, TODO: set min size limit 
  if (fc_fw.deg == 0 && fc_bw.deg == 0) {
    size_t seqlen = fc_fw.s.size() + fc_bw.s.size() - k + 1;
    if (seqlen < 2*k) {
      empty=true;
    }
  }

  if (fc_bw.selfloop == 0) {
    // (selfloop == 0) => Case 0
    // (selfloop == 2) => Case 2a
  } else if (fc_bw.selfloop == 1) { 
    // Case 1
    // Since selfloop != 1, there are two connections from km, into the loop and out of it
    assert(selfloop == 0);
    assert(fc_fw.dist == 1);  
  } else if (fc_bw.selfloop == 2) {
    // Reverse self-loop found on backward strand
    // (selfloop == 0) => Case 2a
    // (selfloop == 2) => Case 2c
    selfloop = fc_bw.selfloop;
  }

  // After: seq == twin(fc_bw.s)[:-k] + fc_fw.s
  if (fc_bw.dist > 1) {
    seq.reserve(fc_bw.s.size() + fc_fw.s.size() - k);
    for (size_t j = fc_bw.s.size() - 1; j >= k; --j) {
      seq.push_back(beta[(fc_bw.s[j] & 7) >> 1]);
    }
    seq += fc_fw.s;
  } else {
    seq = fc_fw.s;
  }

  return MakeContig(seq, selfloop, fc_bw.dist - 1,empty);
}
