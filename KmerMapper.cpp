#include "KmerMapper.hpp"
#include <cmath>
#include <iostream>


static const char alpha[4] = {'A','C','G','T'};

// use:  reverse(s);
// pre:  s != NULL
// post: s has been reversed
void reverse(uint8_t *s, int len) {
  for (int i=0; i<len/2; i++) {
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
  mapContig(id, cr.ref.contig->seq.size()-Kmer::k+1, s);
  return id;
}


// use:  mapper.mapContig(id, len, s);
// pre:  the contig whose string sequence is s has not been mapped before
// post: the contig has been mapped
void KmerMapper::mapContig(uint32_t id, size_t len, const char *s) {
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
    map.insert(make_pair(rep, ContigRef(id, ipos)));    
  }
  if (!last) {
    pos = len-1;
    Kmer km(s+pos);
    Kmer rep = km.rep();
    ipos  = (km == rep) ? (int32_t) pos : -((int32_t)(pos+Kmer::k-1));
    map.insert(make_pair(rep, ContigRef(id, ipos)));  
  }
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

  int a_direction = 0, b_direction = 0;
  Kmer aFirst = sa.getKmer(0).twin();
  Kmer aLast = sa.getKmer(sa.size() - k);
  Kmer bFirst = sb.getKmer(0);
  Kmer bLast = sb.getKmer(sb.size() - k).twin();

  if (isNeighbor(aLast,bFirst)) {
    a_direction = 1;
    b_direction = 1;
  } else if (isNeighbor(aLast,bLast)) {
    a_direction = 1;
    b_direction = -1;
  } else if (isNeighbor(aFirst, bLast)) {
    a_direction = -1;
    b_direction = -1;
  } else if (isNeighbor(aFirst, bFirst)) {
    a_direction = -1;
    b_direction = 1;
  } else {
    cerr << "sa:" << endl << sa.toString() << endl << endl;
    cerr << "sb:" << endl << sb.toString() << endl << endl;
    assert(0);
  }

  Contig *joined = new Contig(); // allocate new contig
  joined->seq.reserveLength(sa.size() + sb.size() - k + 1);
  joined->seq.setSequence(sa, 0, sa.size(), 0, a_direction == -1); // copy all from a, keep orientation of a
  joined->seq.setSequence(sb, k - 1, sb.size() - k + 1, sa.size(), b_direction == -1); // copy from b, reverse if neccessary

  joined->allocateCov(true); // true because the joined contig will have full coverage

  // reverse coverages if neccessary
  
  /*
  if (a_direction == -1) {
    reverse(ca->cov, ca->covlength);
  }
  if (b_direction == -1) {
    reverse(cb->cov, cb->covlength);
  }
  
  for(unsigned int i=0; i <= sa.size() - k; ++i) {
    joined->cov[i] = ca->cov[i]; 
  }
  
  for(unsigned int i=0; i <= sb.size() - k; ++i) {
    joined->cov[sa.size() - k + 1 + i] = cb->cov[i]; 
  } 
  */

  ContigRef cr;
  cr.ref.contig = joined;
  uint32_t id = (uint32_t) contigs.size();
  contigs.push_back(cr); // add to contigs set, is now at position id

  size_t sa_size = sa.size();
  size_t sb_size = sb.size();
  
  // invalidated old contigs
  delete contigs[a_id].ref.contig;
  delete contigs[b_id].ref.contig;
  
  if (a_direction == 1) {
    contigs[a_id] = ContigRef(id, 0);
  } else {
    contigs[a_id] = ContigRef(id, -sa_size + k);
  }
  
  if (b_direction == 1) {
    contigs[b_id] = ContigRef(id, sa_size - k + 1);
  } else {
    contigs[b_id] = ContigRef(id, 1 - sa_size - sb_size);
  }


  
  
  // TODO: fix stride issues, release k-mers, might improve memory
  assert(!contigs[a_id].isContig);
  assert(!contigs[b_id].isContig);
  assert(contigs[id].isContig);
  assert(!contigs[id].isEmpty());

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
      cout << "contig " << id << ": length "  << s.size() << endl;
      cout << s << endl;
      cout << "kmers mapping: " << endl;
      const char *t = s.c_str();
      char tmp[Kmer::MAX_K+1];
      for (size_t i = 0; i < s.length()-Kmer::k+1; i++) {
        Kmer km(t+i);
        if (!find(km).isEmpty()) {
          km.rep().toString(tmp);
          ContigRef km_rep = find(km);
          cout << string(i,' ') << tmp << " -> (" << km_rep.ref.idpos.id << ", " << km_rep.ref.idpos.pos << ")"  << endl;
        }
      }
    } else {
      ContigRef rep = find_rep(a);
      cout << "-> (" << rep.ref.idpos.id << ", " << rep.ref.idpos.pos << ")" << endl;
    }
  }
}


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

    pos += abs(b.ref.idpos.pos);
    int sign = (pos >= 0) ? 1 : -1;
    sign *= (b.ref.idpos.pos >= 0) ? 1 : -1;
    pos *= sign;
    a = b;
  }
  return ContigRef(id, pos);
}


// use:  d = mapper.splitAndJoinContigs();
// pre:  
// post: The contigs in mapper have been splitted and joined
//       d is the increase of contigs after split and join
int KmerMapper::splitAndJoinContigs() {
  Kmer km, rep, end, km_del;
  km_del.set_deleted();
  map.set_deleted_key(km_del);

  int splitted = splitContigs();
  int joined = joinContigs();
  return splitted - joined;
}


// use:  joined = mapper.joinContigs()
// pre:  
// post: contigs that really should be connected have been connected 
//       joined is the number of contigs joined
int KmerMapper::joinContigs() {
  size_t joined = 0;
  
  Contig *c;
  ContigRef cr;

  for(size_t contigid = 0; contigid < contigs.size(); ++contigid) {
    cr = contigs[contigid];
    if (!cr.isContig || cr.isEmpty()) {
      continue;
    }
    c = cr.ref.contig;
    
    ContigRef found;
    Kmer start_twin = c->seq.getKmer(0).twin();
    Kmer end = c->seq.getKmer(c->covlength-1);
    if (checkContigForward(c, end, found)) {
      joinContigs(ContigRef(contigid, 0), found); // this -> found
      ++joined;
    } else if (checkContigForward(c, start_twin, found)) {
      joinContigs(found, ContigRef(contigid, 0)); // found -> this
      ++joined;
    }
  }
  return joined;
}

bool KmerMapper::checkContigForward(Contig* c, Kmer km, ContigRef &found) {
  ContigRef b,cand;
  Kmer fw_km;
  size_t fw_count=0,bw_count=0;
  for(size_t i = 0; i < 4; ++i) {
    Kmer fw = km.forwardBase(alpha[i]);
    ContigRef b = find(fw);
    if (!b.isEmpty()) {
      fw_count += 1;
      cand = b;
      fw_km = fw;
    }
  }

  if (fw_count == 1 && getContig(cand).ref.contig != c) { // one fw-neighbor and no self-loop
    for (size_t i = 0; i < 4; i++) {
      Kmer bw = fw_km.backwardBase(alpha[i]);
      b = find(bw);
      if (!b.isEmpty()) {
        bw_count += 1;
      }
    }
    if (bw_count == 1) {
      found = cand;
      return true;
    }
  }
  return false;
}


// use:  count = mapper.splitContigs()
// pre:  
// post: all contigs with 1 coverage somewhere have been split on those locations
//       count is the number of contigs splitted
int KmerMapper::splitContigs() {
  size_t splitted = 0;
  size_t k = Kmer::k, contigcount = contigs.size();
  size_t covlength, seqlength, cstr_len = 2*k+1;
  size_t nextid = contigcount;
  char *cstr = (char*) malloc(cstr_len);
  bool okay;

  Contig *c;
  ContigRef cr;

  for(size_t contigid = 0; contigid < contigcount; ++contigid) {
    cr = contigs[contigid];
    
    if (!cr.isContig) {
      continue;
    }

    c = cr.ref.contig;
    covlength = c->covlength;

    // Check if this contig has 1 coverage somewhere
    /*
    okay = true;
    for(size_t index=0; index < covlength; ++index) {
      if (c->cov[index] <= 1) {
        okay = false;
        break;
      }
    }
    */
    okay = c->covp->isFull();
    
    if (okay) {
      // The new CompressedSequence class
      //assert(c->covp->isFull());
      continue;
    }
    
    // The new CompressedSequence class
    //assert(!c->covp->isFull());
    
    seqlength = c->seq.size();

    if (seqlength >= cstr_len) {
      cstr_len = 2*seqlength+1;
      cstr = (char*) realloc(cstr, cstr_len);
    }
    
    strcpy(cstr, c->seq.toString().c_str());
    cstr[seqlength] = 0;

    vector<pair<int, int> > v = c->covp->getSplittingVector();
    /*
    size_t a = 0, b = 0;
    vector<pair<int, int> > v;
    // The new CompressedSequence class

    // put [start,end] of covered subintervals of cstr into v
    while (b != covlength) {
      while (a < covlength &&  c->cov[a] <= 1) {
        a++;
      }
      if (a == covlength) {
        break;
      }
      b = a;
      while (b < covlength && c->cov[b] > 1) {
       b++; 
      }
      v.push_back(make_pair(a,b));
      a = b;
    }
    */
    // The new CompressedSequence class
    // assert(v == splittingVector)

    // unmap the contig
    for(size_t index = 0; index < covlength; ++index) {
      if (index % stride == 0) {
        Kmer km = Kmer(&cstr[index]);
        Kmer rep = km.rep();
        map.erase(rep);
      }
    }
    if ((covlength - 1) % stride != 0) {
      Kmer km = Kmer(&cstr[covlength-1]);
      Kmer rep = km.rep();
      map.erase(rep);
    }
    

    // add the subcontigs to contigs and map them
    splitted += v.size() - 1;
    for(size_t index = 0; index < v.size(); ++index) {
      size_t a = v[index].first, b = v[index].second;
      string s(&cstr[a], (b - a) + k - 1);
      ContigRef newcr;
      Contig *newc = new Contig(s.c_str(), true);
      newcr.ref.contig = newc;
      /*
      size_t i = 0;
      while (a < b) {
        newc->cov[i++] = c->cov[a++];
      }
      */
      contigs.push_back(newcr);
      mapContig(nextid++, newc->covlength, s.c_str());
    }

    delete c;
    contigs[contigid] = ContigRef();
  }
  free(cstr);
  return splitted;
}


// use:  mapper.printContigs()
// pre:  
// post: All the contigs in mapper have been printed, line by line, to stdout
void KmerMapper::printContigs() {
  size_t contigcount = contigs.size();

  for(size_t contigid = 0; contigid < contigcount; ++contigid) {
    ContigRef cr = contigs[contigid];
    if (cr.isContig && !cr.isEmpty()) {
      cout << cr.ref.contig->seq.toString() << endl;
    }
  }
}
