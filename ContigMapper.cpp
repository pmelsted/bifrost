#include "ContigMapper.hpp"
#include "CompressedSequence.hpp"
#include <string>
#include <iterator>
#include <algorithm>


// for debugging
#include <iostream>

size_t stringMatch(const string &a, const string& b, size_t pos) {
  return distance(a.begin(),mismatch(a.begin(), a.end(), b.begin() + pos).first);
}

// use: delete cm
// pre:  
// post: all memory allocated has been released
ContigMapper::~ContigMapper() {
  // we do not own bf pointer
  // long contigs could have pointers, but those should already be released
  hmap_long_contig_t::iterator it = lContigs.begin(),it_end=lContigs.end();
  for (; it != it_end; ++it) {
    delete it->second; // the contig pointer
  }
}

// user: i = cm.contigCount()
// pre:  
// post: i is the number of contigs in the mapper
size_t ContigMapper::contigCount() const {
  return lContigs.size() + sContigs.size();
}

// use:  cm.mapBloomFilter(bf)
// pre:  bf != null
// post: uses the bloom filter bf to map reads
void ContigMapper::mapBloomFilter(const BloomFilter* bf) {
  this->bf = bf;
}

// use:  cm.mapRead(km,pos,cc)
// pre:  cc is a reference to a current contig in cm, km maps to cc
// post: the coverage information in cc has been updated
void ContigMapper::mapRead(const ContigMap& cc) {
  if (cc.isEmpty) { // nothing maps, move on
    return;
  } else if (cc.isShort) { 
    // find short contig
    hmap_short_contig_t::iterator it = sContigs.find(cc.head);
    CompressedCoverage& cov = it->second; // reference to the info
    // increase coverage
    cov.cover(cc.dist, cc.dist+ ((cc.repequal) ? cc.len -1 : -cc.len+1));
  } else {
    // find long contig
    hmap_long_contig_t::iterator it = lContigs.find(cc.head);
    Contig* cont = it->second;
    if (cc.repequal) {
      cont->cover(cc.dist, cc.dist + cc.len -1);
    } else {
      cont->cover(cc.dist -cc.len+2, cc.dist);
    }
  }
}

// use: b = cm.addContig(km,read)
// pre:
// post: either contig string containsin has been added and b == true
//       or it was present and the coverage information was updated, b == false
//       NOT Threadsafe!
bool ContigMapper::addContig(Kmer km, const string& read, size_t pos) {
  // find the contig string to add
  string s;
  findContigSequence(km,s);
  size_t k = Kmer::k;
  
  // head is the front 
  const char* c = s.c_str();
  size_t len = s.size();

  Kmer headrep = Kmer(c).rep();
  
  ContigMap cc = this->find(headrep);
  bool found = !cc.isEmpty;

  if (!found) {
    // proper new contig
    if (s.size() < limit) {
      // create a short contig
      sContigs.insert(make_pair(headrep, CompressedCoverage(s.size()-k+1)));
    } else {
      lContigs.insert(make_pair(headrep, new Contig(c)));
      
      for (size_t i = k; i < len-k+1; i += k) {
        shortcuts.insert(make_pair(Kmer(c+i),make_pair(headrep,i)));
      }
      
    }
  }
  
  
  // map the read
  cc = findContig(km,read, pos);
  mapRead(cc);

  return found;
}


// use:  s = cm.findContigSequence(km)
// pre:  km is in the bloom filter
// post: s is the contig containing the kmer km
//       and the rep of the first k-mer in s is smaller (wrt. < operator)
//       rep of the last kmer
void ContigMapper::findContigSequence(Kmer km, string& s) {
  string fw_s;
  Kmer end = km;
  char c;
  while (fwBfStep(end,end,c)) {
    fw_s.push_back(c);
  }
  string bw_s;
  Kmer front = km;
  while (bwBfStep(front,front,c)) {
    bw_s.push_back(c);
  }
  reverse(bw_s.begin(), bw_s.end());
  
  size_t k = Kmer::k;
  s.reserve(k + fw_s.size()+bw_s.size());
  s.append(bw_s);

  char tmp[Kmer::MAX_K];
  km.toString(tmp);
  s.append(tmp);

  s.append(fw_s);

  const char *t = s.c_str();
  Kmer head(t);
  Kmer tail(t+s.size()-k);
  if (tail.rep() < head.rep()) { // reverse complement the string
    s = CompressedSequence(s).rev().toString();
  }
}


// use:  cc = cm.findContig(km)
// pre:  s[pos,pos+k-1] is the kmer km
// post: cc contains either the reference to the contig position
//       or empty if none found
ContigMap ContigMapper::findContig(Kmer km, const string& s, size_t pos) const {
  assert(bf != NULL);
  size_t k = Kmer::k;
  
  // dist is the distance to the beginning.
  size_t dist = 0;
  Kmer end = km;

  // need to check if we find it right away
  ContigMap cc;
  cc = this->find(end);

  char c;
  string fw_s;
  // check <k steps ahead in fw direction
  while (cc.isEmpty && dist < k && fwBfStep(end, end, c)) {
    fw_s.push_back(c);
    cc = this->find(end);
    ++dist;
  }

  
  if (!cc.isEmpty) {
    //TODO: check this logic
    if (cc.isShort) {
      // figure out length of match.
      size_t len = 1 + stringMatch(fw_s, s, pos+k);
      return ContigMap(cc.head, dist, cc.size, len, end == cc.head, true);
    } else {
      size_t len = 1 + stringMatch(fw_s, s, pos+k);

      if (len == 1+dist) {
        // end is now the matching k-mer
        const CompressedSequence& seq = lContigs.find(cc.head)->second->seq;
        //len += seq.jump(s.c_str(), pos+dist+k, cc.dist+k, cc.repequal) -1;
        //cout << "end kmer:     " << end.toString() << endl;
        //cout << "read match:   "  <<s.substr(pos+dist) << endl;
        //cout << "contig match: " << seq.toString().substr(cc.dist) << endl;
        //cout << "contig:       " << seq.toString() << endl;
        //cout << "len before" << len << endl;
        //len += stringMatch(seq.toString().substr(cc.dist+k),s, pos+dist+k);
        len += seq.jump(s.c_str(), pos+dist+k, cc.dist+k,!cc.repequal); // TODO: check arithmetic here!
        //cout << "len after" << len << endl;

      } else {
        // read does not match the contig string
        
      }
      size_t km_dist = cc.dist + ((cc.repequal) ? -dist : dist);
      return ContigMap(cc.head, km_dist,len,  cc.size, cc.repequal, false);
    }
  } 

  //TODO check logic in backwards matching strings. is it really 1?
  Kmer front = km;
  while (cc.isEmpty && dist < k && bwBfStep(front,front,c)) {
    cc = this->find(front);
    ++dist;
  }

  // 

  if (!cc.isEmpty) {
    if (cc.isShort) {
      return ContigMap(cc.head, dist, cc.size, 1, front == cc.head, true);
    } else {
      dist = cc.dist + (cc.repequal) ? dist : (-dist);
      return ContigMap(cc.head, dist, cc.size, 1, front == cc.head, false);
    }
  }


  // ok, no match but how many steps can we jump ahead
  return ContigMap(1+stringMatch(fw_s,s,pos+k));
}

// use:  b = cm.bwBfStep(km,front,c)
// pre:  km is in the bloom filter
// post: b is true if km is inside a contig, in that 
//       case end is the bw link and c is the nucleotide used for the link.
//       if b is false, front and c are not updated
//       if km is an isolated self link (e.g. 'AAA') i.e. end == km then returns false
bool ContigMapper::bwBfStep(Kmer km, Kmer &front, char &c) const {
  size_t i,j;
  size_t bw_count = 0;
  //size_t k = Kmer::k;

  // check bw direction
  j = -1;
  for (i = 0; i < 4; ++i) {
    Kmer bw_rep = front.backwardBase(alpha[i]).rep();
    if (bf->contains(bw_rep)) {
      j = i;
      ++bw_count;
      if (bw_count > 1) {
        break;
      }
    }
  }

  if (bw_count != 1) {
    return false;
  }

  // only one k-mer in the bw link

  Kmer bw = front.backwardBase(alpha[j]);
  size_t fw_count = 0;
  for (i = 0; i < 4; ++i) {
    Kmer fw_rep = bw.forwardBase(alpha[i]).rep();
    if (bf->contains(fw_rep)) {
      ++fw_count;
      if (fw_count > 1) {
        break;
      }
    }
  }

  assert(fw_count >= 1);
  if (fw_count != 1) {
    return false;
  }

  if (bw != km) {
    // exactly one k-mer in bw, character used is c
    front = bw;
    c = alpha[j];
    return true;
  } else {
    return false;
  }
}

// use:  b = cm.fwBfStep(km,end,c)
// pre:  km is in the bloom filter
// post: b is true if km is inside a contig, in that
//       case end is the fw link and c is the nucleotide used for the link.
//       if b is false, end and c are not updated
//       if km is an isolated self link (e.g. 'AAA') i.e. end == km then returns false
bool ContigMapper::fwBfStep(Kmer km, Kmer &end, char &c) const {
  size_t i,j;
  size_t fw_count = 0;
  //size_t k = Kmer::k;

  // check fw direction
  j = -1;
  for (i = 0; i < 4; ++i) {
    Kmer fw_rep = end.forwardBase(alpha[i]).rep();
    if (bf->contains(fw_rep)) {
      j = i;
      ++fw_count;
      if (fw_count > 1) {
        break;
      }
    }
  }
  
  if (fw_count != 1) {
    return false;
  }
  // only one k-mer in fw link

  Kmer fw = end.forwardBase(alpha[j]);

  // check bw from fw link
  size_t bw_count = 0;
  for (i = 0; i < 4; ++i) {
    Kmer bw_rep = fw.backwardBase(alpha[i]).rep();
    if (bf->contains(bw_rep)) {
      ++bw_count;
      if (bw_count > 1) {
        break;
      }
    }
  }
  
  assert(bw_count >= 1);
  if (bw_count != 1) {
    return false;
  }

  if (fw != km) {
    // exactly one k-mer fw, character used is c
    end = fw;
    c = alpha[j];
    return true;
  } else {
    return false;
  }
    
}

// use:  cc = cm.find(km)
// pre:
// post: cc is not empty if there is some info about km
//       in the contig map.
ContigMap ContigMapper::find(Kmer km) const {
  hmap_short_contig_t::const_iterator sit;
  hmap_long_contig_t::const_iterator lit;
  hmap_shortcut_t::const_iterator sc_it;

  Kmer rep = km.rep();
  
  if ((sit = sContigs.find(rep)) != sContigs.end()) {
    return ContigMap(rep, 0, sit->second.size(), 1, km == rep,true);
  } else if ((lit = lContigs.find(rep)) != lContigs.end()) {
    return ContigMap(rep, 0, lit->second->length(), 1, km == rep, false);
  } else {
    if ((sc_it = shortcuts.find(km)) != shortcuts.end()) {
      lit = lContigs.find(sc_it->second.first);
      return ContigMap(lit->first,sc_it->second.second, lit->second->ccov.size(), 1, true,false); // found on fw strand
    } else if ((sc_it = shortcuts.find(rep)) != shortcuts.end()) {
      lit = lContigs.find(sc_it->second.first);
      return ContigMap(lit->first,sc_it->second.second, lit->second->ccov.size(), 1, false,false); // found on rev strand
    }
  }
  return ContigMap();
  
}


// use:  <<splitted, deleted>, joined> = mapper.splitAndJoinAllContigs();
// post: The contigs in mapper have been splitted and joined
//       splitted: the number of splitted contigs 
//       deleted: the number of deleted contigs 
//       joined: the number of joined contigs 
pair<pair<size_t, size_t>, size_t> ContigMapper::splitAndJoinAllContigs() {
  Kmer km, rep, end, km_del;

  // Set the deleted key so we can unmap contigs in splitAllContigs 
  km_del.set_deleted();
  
  lContigs.set_deleted_key(km_del);
  sContigs.set_deleted_key(km_del);
  shortcuts.set_deleted_key(km_del);


  pair<size_t, size_t> splitpair = splitAllContigs();
  size_t joined = joinAllContigs();
  return make_pair(splitpair, joined);
}


// use:  joined = mapper.joinAllContigs()
// post: all contigs that could be connected have been connected 
//       joined is the number of joined contigs
size_t ContigMapper::joinAllContigs() {
  size_t joined = 0;

  /*
  for(size_t contigid = 0; contigid < contigs.size(); ++contigid) {
    cr = contigs[contigid];
    if (!cr.isContig || cr.isEmpty()) {
      continue;
    }
    c = cr.ref.contig;
    
    ContigRef found;
    Kmer start_twin = c->seq.getKmer(0).twin();
    Kmer end = c->seq.getKmer(c->numKmers()-1);
    int dir = 0;
    if ((dir = checkContigForward(c, end, found)) != 0) {
      joinTwoContigs(ContigRef(contigid, 0), found, 1, dir); // this -> found
    } else if ((dir = checkContigForward(c, start_twin, found)) != 0) {
      joinTwoContigs(found, ContigRef(contigid, 0), -dir, 1); // found -> this, -dir because we used twin(first)
    }
    joined += (dir != 0); // increase joined by 1 if (dir == ±1) else increase by 0
  }
  */

  return joined;
}


// use:  splitted, deleted = mapper.splitAllContigs()
// post: All contigs with 1 coverage somewhere have been split where the coverage is 1
//       splitted is the number of contigs splitted
//       deleted is the number of contigs deleted
//       Now every contig in mapper has coverage >= 2 everywhere
pair<size_t, size_t> ContigMapper::splitAllContigs() {
  size_t k = Kmer::k;

  size_t l_contigcount = lContigs.size();
  size_t s_contigcount = sContigs.size();


  // for each short-contig
  typedef vector<pair<int,int> > split_vector_t;
  vector<string> split_contigs; 
  for (hmap_short_contig_t::iterator it = sContigs.begin(); it != sContigs.end(); ) {
    // check if we should split it up
    if (! it->second.isFull()) {
      string s;
      findContigSequence(it->first,s);
      split_vector_t sp = it->second.splittingVector();

      // remember small contigs
      for (split_vector_t::iterator sit = sp.begin(); sit != sp.end(); ++sit) {
	size_t pos = sit->first;
	size_t len = sit->second - pos;
	split_contigs.push_back(s.substr(pos,len+k));	
      }
      
      // erase the split contig
      sContigs.erase(it++); // note: post-increment
    } else {
      ++it;
    }
  }
  
  // insert short contigs
  for (vector<string>::iterator it = split_contigs.begin(); it != split_contigs.end(); ++it) {
    if (it->size() >= k) {
      const char *s = it->c_str();
      Kmer head = Kmer(s).rep();
      Kmer tail = Kmer((s + it->size()-k)).rep();
      if (tail < head) {
	swap(head,tail);
      } 
      // insert contigs into the long contigs, so that the sequence is stored!
      Contig *cont = new Contig(s,true);
      cont->coveragesum = 2 * (it->size()-k); // fake sum, TODO: keep track of this!
      lContigs.insert(make_pair(head,cont));
    }
  }

  split_contigs.clear();



  // long contigs
  vector<pair<string, uint64_t> > long_split_contigs;
  for (hmap_long_contig_t::iterator it = lContigs.begin(); it != lContigs.end(); ) {
    if (! it->second->ccov.isFull()) {
      const string &s = it->second->seq.toString();
      CompressedCoverage &ccov = it->second->ccov;
      pair<size_t, size_t> lowpair = ccov.lowCoverageInfo();
      size_t lowcount = lowpair.first;
      size_t lowsum = lowpair.second;
      size_t totalcoverage = it->second->coveragesum - lowsum;

      // remember pieces
      split_vector_t sp = it->second->ccov.splittingVector();
      for (split_vector_t::iterator sit = sp.begin(); sit != sp.end(); ++sit) {
	size_t pos = sit->first;
	size_t len = sit->second - pos;
	long_split_contigs.push_back(make_pair(s.substr(pos,len+k),(totalcoverage * len)/(ccov.size() - lowcount)));
      }

      // remove shortcuts 
      const char *c = s.c_str();
      for (size_t i = k; i < s.size()-k+1; i += k) {
	shortcuts.erase(Kmer(c+i));
      }

      // erase the split contig
      delete it->second;
      it->second = NULL;
      lContigs.erase(it++); // note: post-increment
    } else {
      ++it;
    }
  }

  // insert the pieces back
  for (vector<pair<string, uint64_t> >::iterator it = long_split_contigs.begin(); it != long_split_contigs.end(); ++it) {
    if (it->first.size() >= k) {
      const char *s = it->first.c_str();
      size_t len = it->first.size()-k+1;
      Kmer head = Kmer(s).rep();
      Kmer tail = Kmer((s+len-1)).rep();
      string tmp;

      if (tail < head) {
	swap(head,tail);
	CompressedSequence c(s);
	tmp = c.rev().toString(); // TODO: create better utility methods!
	s = tmp.c_str();
      }

      // insert new contig
      Contig *cont = new Contig(s,true);
      cont->coveragesum = it->second;
      lContigs.insert(make_pair(head, cont));
      
      // create new shortcuts
      for (size_t i = k; i < len; i+= k) {
	shortcuts.insert(make_pair(Kmer(s+i),make_pair(head,i)));
      }
    }
  }

  long_split_contigs.clear();

  return make_pair(0,0); // TODO: return actual numbers
}


// use:  count2 = mapper.writeContigs(count1, contigfilename, graphfilename);
// pre:  the program has permissions to open contigfilename and graphfilename 
// post: all the contigs have been written to the file: contigfilename
//       the De Brujin graph has been written to the file: graphfilename
//       count2 is the number of real contigs and we assert that count1 == count2 
size_t ContigMapper::writeContigs(int count1, string contigfilename, string graphfilename) {
  /* This is the schema for the outputfiles: 
    --- graphfile:
    contigcount kmersize                    (only in the first line of the file)
    id_length_ratio
    bw1 bw2 bw3 bw4                         (at most 4) // Backward maps
    fw1 fw2 fw3 fw4                         (at most 4) // Forward maps
    ...

    --- contigfile:
    >contigID
    sequence
    ...
  */
  size_t count2 = 0;
  /*
  ofstream contigfile, graphfile;
  contigfile.open(contigfilename.c_str());
  graphfile.open(graphfilename.c_str());
  assert(!contigfile.fail() && !graphfile.fail());

  size_t k = Kmer::k;
  int count2 = 0;
  size_t contigcount = contigs.size();
 
  graphfile << count1 << " " << k << "\n"; 

  for(size_t id = 0; id < contigcount; ++id) {
    ContigRef cr = contigs[id];
    if (!cr.isContig || cr.isEmpty()) {
      continue;
    }
    ++count2;

    Contig *c = cr.ref.contig;
   
    size_t length = c->length(), numkmers = length - k + 1, coveragesum = c->coveragesum;
    float ratio = coveragesum / (0.0 + numkmers);
    
    Kmer first = c->seq.getKmer(0), last = c->seq.getKmer(length - k);
    
    contigfile << ">contig" << id << "\n" << c->seq.toString() << "\n";
    graphfile << id << "_" <<  length << "_" << ratio << "\n";

    
    for (size_t i=0; i<4; ++i) {
      Kmer bw = first.backwardBase(alpha[i]);
      ContigRef prevcr = find(bw);
      if (!prevcr.isEmpty()) {
        Contig *oc = getContig(prevcr).ref.contig;
        size_t oid = prevcr.ref.idpos.id;
        Kmer oFirst = oc->seq.getKmer(0);
        Kmer oLast = oc->seq.getKmer(oc->length() - k);
        if (oid == id) {
          if (isNeighbor(oLast, first)) {
            cerr << "Self-looped, id: " << id << ", seq: " << c->seq.toString() << endl; 
          } else {
            cerr << "Hairpinned, id: " << id << ", seq: " << c->seq.toString() << endl; 
          }
        }
        if (isNeighbor(oLast, first) || isNeighbor(oFirst.twin(), first)) { 
          graphfile << oid << " ";
        } else {
          assert(isNeighbor(oLast.twin(), first) || isNeighbor(oFirst, first));
        }
      }
    }

    graphfile << "\n";

    for (size_t i=0; i<4; ++i) {
      Kmer fw = last.forwardBase(alpha[i]);
      ContigRef fwcr = find(fw);
      if (!fwcr.isEmpty()) {
        Contig *oc = getContig(fwcr).ref.contig;
        size_t oid = fwcr.ref.idpos.id;
        Kmer oFirst = oc->seq.getKmer(0);
        Kmer oLast = oc->seq.getKmer(oc->length() - k);
        if (oid == id) {
          if (isNeighbor(last, oFirst)) {
            cerr << "Self-looped, id: " << id << ", seq: " << c->seq.toString() << endl; 
          } else {
            cerr << "Hairpinned, id: " << id << ", seq: " << c->seq.toString() << endl; 
          }
        }
        if (isNeighbor(last, oFirst) || isNeighbor(last, oLast.twin())) {
          graphfile << oid << " ";
        } else {
          assert(isNeighbor(last, oFirst.twin()) || isNeighbor(last, oLast));
        }
      }
    }

    graphfile << "\n";
  }

  // Flush and close
  contigfile.close(); 
  graphfile.close();  
  */
  return count2;
}
