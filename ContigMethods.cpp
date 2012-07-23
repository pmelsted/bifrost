#include "ContigMethods.hpp"

static const char beta[4] = {'T','G','A','C'}; // c -> beta[(c & 7) >> 1] maps: 'A' <-> 'T', 'C' <-> 'G'

// use:  getMappingInfo(repequal, pos, dist, kmernum, cmppos)
// pre:  
// post: cmppos is the first character after the kmer-match at position pos
void getMappingInfo(const bool repequal, const int32_t pos, const size_t dist, size_t &kmernum, int32_t &cmppos) {
  size_t k = Kmer::k; 
  // Now we find the right location of the kmer inside the contig
  // to increase coverage 
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
      cmppos = -pos + 1 - dist; // Original: (-pos +1 -k) - dist + k
      kmernum = cmppos - k;
    }
  }
}


// use:  cc = check_contig_(bf,km,mapper);
// pre:  
// post: if km does not map to a contig: cc.cr.isEmpty() == true and cc.dist == 0
//       else: km is in a contig which cc.cr maps to and cc.dist is the distance 
//             from km to the mapping location 
//             (cc.eq == true):  km has the same direction as the contig
//             else:  km has the opposite direction to the contig
CheckContig check_contig(BloomFilter &bf, KmerMapper &mapper, Kmer km) {
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
MakeContig make_contig(BloomFilter &bf, KmerMapper &mapper, Kmer km) {
  size_t k = Kmer::k;
  string seq;
  FindContig fc_fw = find_contig_forward(bf, km);

  if (fc_fw.selfloop == 1) {
    return MakeContig(fc_fw.s, 0); 
  } else if (fc_fw.selfloop == 2) {
    FindContig fc_bw = find_contig_forward(bf, km.twin());
    Kmer realfirsttwin(fc_bw.s.substr(fc_bw.s.size() - k, k).c_str());
    Kmer realfirst = realfirsttwin.twin();
    fc_fw = find_contig_forward(bf, realfirst);
    return MakeContig(fc_fw.s, fc_bw.s.size() - k); 
  }

  FindContig fc_bw = find_contig_forward(bf, km.twin());
  ContigRef cr_tw_end = mapper.find(fc_bw.end);
  assert(cr_tw_end.isEmpty());

  if (fc_bw.dist > 1) {
    seq.reserve(fc_bw.s.size() + fc_fw.s.size() - k);
    // copy reverse part of fc_bw.s not including k
    for (size_t j = fc_bw.s.size() - 1; j >= k; --j) {
      seq.push_back(beta[(fc_bw.s[j] & 7) >> 1]);
    }
    seq += fc_fw.s; // append fc_fw.s
  } else {
    seq = fc_fw.s;
  }

  return MakeContig(seq, fc_bw.dist - 1);
}
