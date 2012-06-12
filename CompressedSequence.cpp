#include "Common.hpp"
#include "Kmer.hpp"
#include "CompressedSequence.hpp"


static const char bases[256] = {
  'A','C','G','T','N','N','N','N',  'N','N','N','N','N','N','N','N', 
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N', 
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N', 
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N', 
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N', 
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N', 
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N', 
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N', 
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N', 
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N', 
  'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N'
};

// bits['A'] = bits['a'] = 0 
// bits['C'] = bits['c'] = 1
// bits['G'] = bits['g'] = 2 
// bits['T'] = bits['t'] = 3
static const uint8_t bits[256] = {
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};


// use:  delete c;
// pre:  c is a pointer to a CompressedSequence
// post: the memory which the CompressedSequence had allocated has been freed
CompressedSequence::~CompressedSequence() {
  if (_capacity > 0 && _data != NULL) {
    delete[] _data;
    _data = NULL;
  }
}


// use:  _cs = CompressedSequence(cs);
// pre:   
// post: the DNA string in _cs and is the same as in cs
CompressedSequence::CompressedSequence(const CompressedSequence& o) : _length(0),_capacity(0),_data(0) {
  setSequence(o,0,o._length);
}


// use:  _cs = cs;
// pre:   
// post: the DNA string in _cs is the same as in cs
CompressedSequence& CompressedSequence::operator=(const CompressedSequence& o) {
  setSequence(o,0,o._length);
  return *this;
}
  

// use:  cs = CompressedSequence(s);
// pre:  s has only the characters 'A','C','G' and 'T' and can have any length
// post: the DNA string in cs is now the same as s
CompressedSequence::CompressedSequence(const char *s) : _length(0),_capacity(0),_data(0) {
  if (s != NULL) {
    setSequence(s,strlen(s));
  }
}


// same as above except with a string not a char array 
CompressedSequence::CompressedSequence(const string& s) : _length(0),_capacity(0),_data(0) {
  setSequence(s.c_str(), s.size());
}


// use:  cs = CompressedSequence(km);
// pre:   
// post: the DNA string in cs is now the same as the DNA string in km
CompressedSequence::CompressedSequence(const Kmer &km) : _length(0),_capacity(0),_data(0) {
  setSequence(km, Kmer::k);
}


// use:  c = cs[index];
// pre:  0 <= index < cs.size()
// post: c is character nr. index in the DNA string inside cs
const char CompressedSequence::operator[](size_t index) const {
  size_t i = index / 4;
  size_t j = index % 4;
  size_t idx = ((_data[i]) >> (2*j)) & 0x03;
  return bases[idx];
}


//void CompressedSequence::setSequence(const CompressedSequence &o, size_t length, size_t offset, bool reversed) {
//  setSequence(o,0,length,offset,reversed);
//}


// use:  a.setSequence(b, start, length, offset, reversed);
// pre:  start+length <= b._length, offset <= a._length
// post: copies compressed sequence from b to a (reverse complement if reversed == true)
//       the string copied from b is from [start,...,start+length-1]
//       the positions in a that are updated are [offset,...,offset+length-1]
//       capacity of a might be updated to fit the new string.
void CompressedSequence::setSequence(const CompressedSequence &o, size_t start, size_t length, size_t offset, bool reversed) {
  assert(length + start <= o._length);
  if (round_to_bytes(length+offset) > _capacity) {
    _resize_and_copy(round_to_bytes(length+offset), _length);
  } 
  
  size_t w_index = offset;
  size_t r_index = reversed ? length-start-1 : start;
  size_t wi,wj,ri,rj;
  
  for (size_t i = 0; i < length; i++) {
    wi = w_index / 4;
    wj = w_index % 4;
    ri = r_index / 4;
    rj = r_index % 4;
    _data[wi] &= ~(0x03<<(2*wj)); // clear bits
    uint8_t nucl = (o._data[ri] >> (2*rj)) & 0x03; // nucleotide stored in o
    if (reversed) {
      nucl = 3-nucl; // reverse sequence
    }
    _data[wi] |= (nucl << (2*wj));

    w_index++;
    if (reversed) {
      r_index--;
    } else {
      r_index++;
    }
  }
  if (offset + length > _length) 
    _length = offset + length;
}



// use:  cs.reserveLength(new_length);
// pre:  
// post: The DNA string in cs has space for at least new_length bases 
void CompressedSequence::reserveLength(size_t new_length) {
  if (round_to_bytes(new_length) > _capacity) {
    _resize_and_copy(round_to_bytes(new_length), _length);
  }
}


// use:  cs._resize_and_copy(new_cap, copy_limit);
// pre:  
// post: The DNA string in cs has space for at least new_length bases
//       the first copy_limit characters of cs are the same as before this method
void CompressedSequence::_resize_and_copy(size_t new_cap, size_t copy_limit) {
  if (new_cap > _capacity) {
    char *new_data = new char[new_cap];
    size_t bytes = round_to_bytes(copy_limit); 
    memcpy(new_data,_data,bytes);
    delete[] _data;
    _data = new_data;
    _capacity = new_cap;
  }
}


// use:  a.setSequence(s, length, offset, reversed);
// pre:  length <= strlen(s), offset <= a._length
// post: copies substring s[0,...,length-1] to a (reverse complement if reversed == true)
//       the positions in a that are updated are [offset,...,offset+length-1]
//       capacity of a might be updated to fit the new string.
void CompressedSequence::setSequence(const char *s, size_t length, size_t offset, bool reversed) {
  if(round_to_bytes(length+offset) > _capacity) {
    _resize_and_copy(round_to_bytes(length+offset), _length);
  }

  for (size_t index = offset; index < offset+length; index++) {
    size_t i = index / 4;
    size_t j = index % 4;
    _data[i] &= ~(0x03 << (2*j)); // set bits to 0, default
    uint8_t c = reversed ? bases[0x03-bits[(uint8_t)*(s+length+offset-index-1)]] : *(s+index-offset);
    _data[i] |= (bits[c] << (2*j));
  }
  if (offset + length > _length) 
    _length = offset + length;
}


// use:  cs.setSequence(s, start, length, offset, reversed);
// pre:  0 <= start + length < o.size()
// post: If reversed is false then: cs[offset,...,offset+length-1] = s[0,...,start+length-1]
//       else: cs[offset,...,offset+length-1] is the reverse complement of s[0,...,start+length-1]
void CompressedSequence::setSequence(const string &s, size_t length, size_t offset, bool reversed) {
  setSequence(s.c_str(),length,offset,reversed);
}


// use:  cs.setSequence(km, length, offset, reversed);
// pre:  0 <= length < cs._length,
//       length <= Kmer::k
// post: If reversed is false then: cs[offset,...,offset+length-1] 
//         is the first length characters from the DNA string in km
//       else: cs[offset,...,offset+length-1] is the first length characters from the
//         reverse complement of the DNA string in km
void CompressedSequence::setSequence(const Kmer &km, size_t length, size_t offset, bool reversed) {
  char s[Kmer::MAX_K+1];
  km.toString(&s[0]);
  setSequence(s, length, offset, reversed);
}


// use:  s = cs.toString();
// pre:   
// post: s is the DNA string from cs 
string CompressedSequence::toString() const {
  return toString(0,_length);
}


// use:  s = cs.toString(offset, length);
// pre:  offset + length <= cs.size(),
// post: s is the DNA string from c[offset,...,offset+length-1]
string CompressedSequence::toString(size_t offset, size_t length) const {
  assert(offset+length <= _length);
  string s(length,0);
  size_t i,j,idx;
  for (size_t index = offset; index < offset+length; index++) {
    i = index / 4;
    j = index % 4;
    idx = ((_data[i]) >> (2*j)) & 0x03;
    s[index-offset] = bases[idx];
  }
  return s;
}


// use:  cs.toString(s);
// pre:  s has space for cs.size() characters 
// post: s is the same as the DNA string from cs
void CompressedSequence::toString(char *s) const {
  toString(s,0,_length);
}


// use:  s = cs.toString(offset, length);
// pre:  offset + length <= cs.size()
//       s has space for length characters 
// post: s is the same as cs[offset,...,offset+length-1]
void CompressedSequence::toString(char *s, size_t offset, size_t length) const {
  assert(offset+length <= _length);
  size_t i,j,idx;
  for (size_t index = offset; index < offset+length; index++) {
    i = index / 4;
    j = index % 4;
    idx = ((_data[i]) >> (2*j)) & 0x03;
    s[index-offset] = bases[idx];
  }
  s[length] = 0; // 0-terminated string
}


// use:  km = cs.getKmer(offset);
// pre:  offset + Kmer::k <= cs._length
// post: The DNA string in km is cs[offset,...,offset+Kmer::k-1]
Kmer CompressedSequence::getKmer(size_t offset) const {
  char s[Kmer::MAX_K+1];
  toString(&s[0], offset, Kmer::k);
  return Kmer(s);
}


// use:  _cs = cs.rev();
// pre:   
// post: _cs is the reverse complement CompressedSequence with respect to cs,
//       i.e. if the DNA string in cs is 'GTCA'
//          then the DNA string in _cs is 'TGAC'
CompressedSequence CompressedSequence::rev() const {
  CompressedSequence r;
  r.setSequence(*this, 0, _length, 0, true);
  return r;
}


// use:  j = cs.jump(s,i,pos,reversed)
// pre:  0 <= i < s.length, -1 <= pos < cs._length if reversed true, 0 <= pos <= cs._length if reversed false
// post: if reversed == false: s[i...i+j-1] == cs._data[pos...pos+j-1], 0 <= j <= min(s.length-i, cs._length-pos)
//       if reversed == true : reverse_complement(s[i...i+j-1]) == cs._data[pos-j+1...pos], 0 <= j <= min(s.length-i, pos+1)
size_t CompressedSequence::jump(const char *s, size_t i, int pos, bool reversed) const {
  assert(i>=0);
  assert(pos >= -1);
  assert(0 <= _length - pos); // this prevents -1 <= _length from giving false
  size_t j = 0;
  int dir = (reversed) ? -1 : 1;
  int limit = (reversed) ? -1 : _length; // index limit, lower or upper bound
  size_t a,b,idx;
  for (int index = pos; s[i+j] != 0 && index != limit; index+= dir, j++) {
    assert(index >= 0);
    assert(index < _length);
    a = index / 4;
    b = index % 4;
    idx = ((_data[a]) >> (2*b)) & 0x03;
    if ((!reversed && s[i+j] != bases[idx]) || (reversed && s[i+j] != bases[3-idx])) {
      break;
    }
  }
  return j;
}
