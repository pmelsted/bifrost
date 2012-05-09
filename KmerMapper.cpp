#include "KmerMapper.hpp"
#include <cmath>


KmerMapper::~KmerMapper() {
  vector<ContigRef>::iterator it,it_end;
  it_end = contigs.end();
  for (it = contigs.begin(); it != it_end; ++it) {
    if (it->isContig && it->ref.contig != NULL) {
      delete it->ref.contig;
      it->ref.contig = NULL;
    }
  }
}

ContigRef KmerMapper::addContig(const string &s) {
  addContig(s.c_str());
}

ContigRef KmerMapper::addContig(const char *s) {
  // check that it doesn't map, our responsibility or not?
  ContigRef cr;
  cr.ref.contig = new Contig(s);
  // add to the map
  
  

  return cr;
}

ContigRef KmerMapper::find(const Kmer km) {
  Kmer rep = km.rep();
  int dir = (rep == km) ? 1 : -1;
  iterator it = map.find(km.rep());
  if (it == map.end()) {
    return ContigRef();
  }
  
  ContigRef cr = it->second;
  assert(!cr.isContig); // cannot be a pointer
  
  ContigRef a = find_rep(cr);
  it->second = a; // shorten the tree
  
  a.ref.idpos.pos *= dir; // modify the copy
  return a;
}

// use:  r = m.joinContigs(a,b);
// pre:  a and b are not contig pointers
// post: r is a contigref that points to a newly created contig
//       formed by joining a+b with proper direction, sequences
//       pointed to by a and b have been forwarded to the new contig r

ContigRef KmerMapper::joinContigs(ContigRef a, ContigRef b) {
  //join a to b
  a = find_rep(a);
  uint32_t a_id = a.ref.idpos.id;
  b = find_rep(b);
  uint32_t b_id = b.ref.idpos.id;
  assert(contigs[a_id].isContig && contigs[b_id].isContig);

  // fix this mess
  CompressedSequence &sa = (contigs[a_id].ref)->contig.seq;
  CompressedSequence &sb = (contigs[b_id].ref)->contig.seq;

  int direction = 0;
  if(sa.getKmer(sa.size()-Kmer::k) == sb.getKmer(0)) {
    direction = 1;
  }
  if(sa.getKmer(sa.size()-Kmer::k) == sb.getKmer(sb.size()-Kmer::k).twin()) {
    direction = -1;
  }

  assert(direction != 0); // what if we want b+a?

  Contig *joined = new Contig(0); // allocate new contig
  joined->seq.reserveLength(sa.size()+sb.size());
  joined->seq.setSequence(sa,0,sa.size());
  //TODO: fix bug for reversed
  joined->seq.setSequence(sb,sa.size(),sb.size() + 0000000); //copy sequences

  ContigRef cr;
  cr.ref.contig = joined;
  uint32_t id = (uint32_t) contigs.size();
  contigs.push_back(cr); // add to set
  
  // invalidated old contigs
  delete contigs[a_id].ref.contig;
  delete contigs[b_id].ref.contig;
  
  contigs[a_id] = ContigRef(id,0);
  contigs[b_id] = ContigRef(id,sa.size()); 

  return ContigRef(id,0); // points to newly created contig
}

/* Will we need this at all?
ContigRef KmerMapper::extendContig(ContigRef a, const string &s) {
  
}

ContigRef KmerMapper::extendContig(ContigRef a, const char *s) {
  
}
*/


// use:  r = m.find_rep(a);
// pre:  
// post: if a.isContig is false, then r.isContig is false and
//       m.contig[a.ref.idpos.id].isContig is true, otherwise r == a
// Finds the contigref just before the contig pointer
ContigRef KmerMapper::find_rep(ContigRef a) {
  uint32_t id;
  ContigRef b = a;
  if (a.isContig) {
    return a;
  }
  int32_t pos = a.ref.idpos.pos;

  while (true) {
    id = a.ref.idpos.id;
    b = contigs[id];
    if (b.isContig) {
      break;
    }

    pos += abs(a.ref.idpos.pos);
    int sign = (pos >= 0) ? 1 : -1;
    sign *= (a.ref.idpos.pos >= 0) ? 1 : -1;
    pos *= sign;
    a = ContigRef(id,pos);    
  }
  return a;
}
