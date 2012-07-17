#include "ContigMethods.hpp"

static const char alpha[4] = {'A','C','G','T'};
static const char beta[4] = {'T','G','A','C'}; // c -> beta[(c & 7) >> 1] maps: 'A' <-> 'T', 'C' <-> 'G'

// use:  getMappingInfo(repequal, pos, dist, k, kmernum, cmppos)
// pre:  
// post: cmppos is the first character after the kmer-match at position pos
void getMappingInfo(const bool repequal, const int32_t pos, const size_t dist, const size_t k, size_t &kmernum, int32_t &cmppos) {
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
    //return make_pair(cr, make_pair(dist, end == end.rep())); 
  } else {
    return CheckContig(ContigRef(), 0, 0);
    //return make_pair(ContigRef(), make_pair(0,0));
  }
}


// use:  mc = make_contig(bf, mapper, km, s);
// pre:  km is not contained in a mapped contig in mapper 
// post: Finds the forward and backward limits of the contig
//       which contains km  according to the bloom filter bf and puts it into mc.seq
//       mc.pos is the position where km maps into this contig
MakeContig make_contig(BloomFilter &bf, KmerMapper &mapper, Kmer km) {
  size_t k  = Kmer::k;
  string seq, seq_fw(k, 0), seq_bw(k, 0);
  FindContig fc_fw = find_contig_forward(bf, km, &seq_fw);

  if (fc_fw.selfloop) {
    return MakeContig(seq_fw, 0); 
  }

  FindContig fc_bw = find_contig_forward(bf, km.twin(), &seq_bw);
  ContigRef cr_tw_end = mapper.find(fc_bw.end);
  assert(cr_tw_end.isEmpty());

  if (fc_bw.dist > 1) {
    seq.reserve(seq_bw.size() + seq_fw.size() - k);
    // copy reverse part of seq_bw not including k
    for (size_t j = seq_bw.size() - 1; j >= k; --j) {
      seq.push_back(beta[(seq_bw[j] & 7) >> 1]);
    }
    seq += seq_fw; // append seq_fw
  } else {
    seq = seq_fw;
  }

  return MakeContig(seq, fc_bw.dist - 1);
}


// use:  fc = find_contig_forward(bf,km,s);
// pre:  
// post: km is contained in a contig c with respect to the
//       bloom filter graph bf and fc.end is the forward endpoint (wrt km direction)
//       and c contains fc.dist kmers until the end (including km)
//       if s is not NULL the sequence of the contig is stored in s
FindContig find_contig_forward(BloomFilter &bf, Kmer km, string* s) {
  int j;
  bool selfloop = false;
  size_t i,dist = 1;
  vector<char> v;

  Kmer first = km, end = km;

  assert(bf.contains(km.rep()));
  if (s != NULL) {
    char t[Kmer::MAX_K+1];
    km.toString(t);
    for (i = 0; i < Kmer::k; ++i) {
      v.push_back(t[i]);
    }
  }
  
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
    if (s != NULL) {
      v.push_back(alpha[j]);
    }
  }

  if (s != NULL) {
    s->clear();
    s->reserve(v.size());
    s->insert(s->begin(), v.begin(), v.end());
  }
  return FindContig(end, dist, selfloop);
}

