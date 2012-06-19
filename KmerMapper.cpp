#include "KmerMapper.hpp"
#include <cmath>
#include <iostream>


static const char alpha[4] = {'A','C','G','T'};

// use:  reverse(s);
// pre:  s != NULL
// post: s has been reversed
void reverse(uint8_t *s, int len) {
  for (int i=0;i<len/2;i++) {
    s[i]^=s[len-i-1];                                                                                                                        
    s[len-i-1]^=s[i];
    s[i]^=s[len-i-1];
  }
}

// use:  delete m;
// pre:  m is a pointer to a KmerMapper
// post: the memory that this KmerMapper had allocated has been freed 
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


// same as addContig(const char *s) but with string
size_t KmerMapper::addContig(const string &s) {
  return addContig(s.c_str());
}


// use:  id = mapper.addContig(s);
// pre:  s is a string of 'A','C','G' and 'T's
// post: a contig with s as string has been added to mapper
//       the reps of the first and last kmer in this contig now map to this contig
//       the reps of the kmers between the first and last that do not overlap each other
//       also map to this contig
//       id is the id of the contig
size_t KmerMapper::addContig(const char *s) {
  // check that it doesn't map, our responsibility or not?
  ContigRef cr;
  cr.ref.contig = new Contig(s);
   
  uint32_t id = (uint32_t) contigs.size();
  contigs.push_back(cr);

  size_t len = cr.ref.contig->seq.size()-Kmer::k+1;
  bool last = false;
  size_t pos;
  int32_t ipos;
  for (pos = 0; pos < len; pos += stride) {
    if (pos == len-1) {
      last = true;
    }
    Kmer km(s+pos);
    Kmer rep = km.rep();

    ipos  = (km == rep) ? (int32_t) pos : -((int32_t)(pos+Kmer::k-1));
    map.insert(make_pair(rep,ContigRef(id,ipos)));    
  }
  if (!last) {
    pos = len-1;
    Kmer km(s+pos);
    Kmer rep = km.rep();
    ipos  = (km == rep) ? (int32_t) pos : -((int32_t)(pos+Kmer::k-1));
    map.insert(make_pair(rep,ContigRef(id,ipos)));  
  }
  return id;
}


// use:  cr = mapper.find(km);
// pre:  
// post: If the rep of kmer km maps to a contig, cr is the contigref that maps the rep
//       to a contig, else cr is an empty contigref 
ContigRef KmerMapper::find(const Kmer km) {
  iterator it = map.find(km.rep());
  if (it == map.end()) {
    return ContigRef();
  }
  
  ContigRef cr = it->second;
  assert(!cr.isContig); // cannot be a pointer
  
  ContigRef a = find_rep(cr);
  it->second = a; // shorten the tree for future reference

  return a;
}

// use:  r = m.joinContigs(a, b);
// pre:  a and b are not contig pointers
//       the last Kmer::k-1 bases in the last kmer in the contig that a refers to
//       are the same as the first Kmer::k-1 bases in the first kmer or the twin of 
//       the last kmer in the contig that b refers to
// post: r is a contigref that points to a newly created contig
//       formed by joining a+b with proper direction, sequences
//       pointed to by a and b have been forwarded to the new contig r
ContigRef KmerMapper::joinContigs(ContigRef a, ContigRef b) {
  //join a to b
  size_t k = Kmer::k;
  a = find_rep(a);
  uint32_t a_id = a.ref.idpos.id;
  b = find_rep(b);
  uint32_t b_id = b.ref.idpos.id;
  assert(contigs[a_id].isContig && contigs[b_id].isContig);

  Contig *ca = contigs[a_id].ref.contig;
  Contig *cb = contigs[b_id].ref.contig;
  CompressedSequence &sa = ca->seq;
  CompressedSequence &sb = cb->seq;

  int direction = 0;
  Kmer aLast = sa.getKmer(sa.size() - k);
  Kmer bFirst = sb.getKmer(0);
  Kmer bLast = sb.getKmer(sb.size() - k).twin();
  Kmer next;
  for(int i=0; i<4; ++i) {
    next = aLast.forwardBase(alpha[i]);
    if (next == bFirst) {
      direction = 1;
      break;
    }
    if (next == bLast) {
      direction = -1;
      break;
    }
  }

  assert(direction != 0);

  Contig *joined = new Contig(0); // allocate new contig
  joined->seq.reserveLength(sa.size() + sb.size() - k + 1);
  joined->seq.setSequence(sa, 0, sa.size(), 0, false); // copy all from a, keep orientation of a
  joined->seq.setSequence(sb, k - 1, sb.size() - k + 1, sa.size(), direction == -1); // copy from b, reverse if neccessary

  joined->allocateCov();
 
  for(unsigned int i=0; i <= sa.size() - k; ++i) {
    joined->cov[i] = ca->cov[i]; 
  } 

  if (direction == -1) {
    reverse(cb->cov, cb->covlength);
  }

  for(unsigned int i=0; i <= sb.size() - k; ++i) {
    joined->cov[sa.size() - k + 1 + i] = cb->cov[i]; 
  } 

  ContigRef cr;
  cr.ref.contig = joined;
  uint32_t id = (uint32_t) contigs.size();
  contigs.push_back(cr); // add to contigs set, is now at position id
  
  // invalidated old contigs
  delete contigs[a_id].ref.contig;
  delete contigs[b_id].ref.contig;
  
  contigs[a_id] = ContigRef(id, 0);
  contigs[b_id] = ContigRef(id, sa.size()); 

  // TODO: fix stride issues, release k-mers, might improve memory

  return ContigRef(id, 0); // points to newly created contig
}


// use:  contig = mapper.getContig(_id);
// pre:  
// post: contig is the contig with id _id
ContigRef KmerMapper::getContig(const size_t id) const {
  ContigRef a = contigs[id];
  if (!a.isContig) {
    a = find_rep(a);
    a = contigs[a.ref.idpos.id];
  }
  return a;
}


// use:  contig = mapper.getContig(cr);
// pre:  
// post: contig is the contig that cr maps to 
ContigRef KmerMapper::getContig(const ContigRef ref) const {
  if (ref.isContig) {
    return ref;
  } else {
    return getContig(ref.ref.idpos.id);
  }
}


// use:  mapper.printContig(_id);
// pre:  _id is in mapper.contigs
// post: details about the contig whose id is _id has been printed to cout
void KmerMapper::printContig(const size_t id) {
  if (id >= contigs.size()) {
    cerr << "invalid reference " << id << endl;
  } else {
    ContigRef a = contigs[id];
    if (a.isContig) {
      string s = a.ref.contig->seq.toString();
      //cout << "contig " << id << ": length "  << s.size() << endl;
      cout << s << endl;
      //cout << "kmers mapping: " << endl;
      const char *t = s.c_str();
      char tmp[Kmer::MAX_K+1];
      for (int i = 0; i < s.length()-Kmer::k+1; i++) {
        Kmer km(t+i);
        if (!find(km).isEmpty()) {
          km.rep().toString(tmp);
          ContigRef km_rep = find(km);
          //cout << string(i,' ') << tmp << " -> (" << km_rep.ref.idpos.id << ", " << km_rep.ref.idpos.pos << ")"  << endl;
        }
      }
    } else {
      ContigRef rep = find_rep(a);
      //cout << "-> (" << rep.ref.idpos.id << ", " << rep.ref.idpos.pos << ")" << endl;
    }
  }
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
ContigRef KmerMapper::find_rep(ContigRef a) const {
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


// use:  mapper.splitAndJoinContigs();
// pre:  No contig is longer than 8192 
// post: The contigs in mapper have been splitted and joined
void KmerMapper::splitAndJoinContigs() {
  Kmer km, rep, end, km_del;
  km_del.set_deleted();
  map.set_deleted_key(km_del);

  size_t firstchar, lastchar, contigcount = contigs.size();
  size_t lengthbefore, k = Kmer::k;
  uint32_t covlength;
  char cstr[8192];
  char *p;
  uint8_t *covp;
  Contig *now;
  ContigRef cr;
  
  for(size_t contigid = 0; contigid < contigcount; ++contigid) {
    cr = contigs[contigid];
    if (!cr.isContig) {
      continue;
    }

    p = &cstr[0];
    now = cr.ref.contig;
    firstchar = 0;
    lastchar = now->seq.size() - 1;
    assert(lastchar < 8192);
    covlength = now->covlength;
    assert(covlength + k - 2 == lastchar);
    end = Kmer(&cstr[covlength-1]);

    strcpy(cstr, now->seq.toString().c_str());
    cstr[lastchar + 1] = 0;

    // Trim the contig if either end is only covered once 
    while (1 + lastchar - firstchar >= k && now->cov[firstchar] == 1) {
      if (firstchar % stride == 0) {
        km = Kmer(p);
        rep = km.rep();
        map.erase(rep);
        iterator it = map.find(rep);
        assert(it == map.end());
      }
      ++p;
      ++firstchar;
    }

    lengthbefore = covlength;
    while (1 + lastchar - firstchar >= k && now->cov[covlength-1] == 1) {
      assert(lastchar - k +1 == covlength -1);
      if ( (covlength -1) % stride == 0) {
        km = Kmer(&cstr[covlength-1]);
        rep = km.rep();
        map.erase(rep);
        iterator it = map.find(rep);
        assert(it == map.end());
      }

      cstr[lastchar] = 0;
      --lastchar;
      --covlength;

    }

    if (firstchar > 0 || covlength != lengthbefore) {
      if ((lengthbefore -1) % stride != 0) {
        rep = end.rep();
        map.erase(rep);
        iterator it = map.find(rep);
        assert(it == map.end());
      }
      now->seq.clear();
      now->seq.setSequence(p, lastchar - firstchar + 1, 0, false);
      covp = now->cov;
      now->cov = new uint8_t[covlength - firstchar];
      for(size_t q=0; q < covlength - firstchar; ++q) {
        now->cov[q] = covp[firstchar+q];
      }
      delete[] covp;

      // TODO: Remap the contig
    }

    assert(strncmp(now->seq.toString().c_str(), p, lastchar-firstchar+1) == 0);

    if (1 + lastchar - firstchar < k) {
      // TODO: Delete the contig contigs[contigid].ref.contig 
      contigs[contigid] = ContigRef();
    }
  }
}

void KmerMapper::printContigs() {
  size_t contigcount = contigs.size();
  ContigRef cr; 

  for(size_t contigid = 0; contigid < contigcount; ++contigid) {
    cr = contigs[contigid];
    if (!cr.isContig) {
      continue;
    }
    cout << cr.ref.contig->seq.toString() << endl;
  }
}
