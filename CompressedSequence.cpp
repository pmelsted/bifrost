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


CompressedSequence::~CompressedSequence() {
  if (_capacity > 0 && _data != NULL) {
    delete[] _data;
    _data = NULL;
  }
}


CompressedSequence::CompressedSequence(const CompressedSequence& o) : _length(0),_capacity(0),_data(0) {
  setSequence(o,0,o._length);
}

CompressedSequence& CompressedSequence::operator=(const CompressedSequence& o) {
  setSequence(o,0,o._length);
  return *this;
}
  
CompressedSequence::CompressedSequence(const char *s) : _length(0),_capacity(0),_data(0) {
  if (s != NULL) {
    setSequence(s,strlen(s));
  }
}

CompressedSequence::CompressedSequence(const string& s) : _length(0),_capacity(0),_data(0) {
  setSequence(s.c_str(), s.size());
}

CompressedSequence::CompressedSequence(const Kmer &km) : _length(0),_capacity(0),_data(0) {
  setSequence(km, Kmer::k);
}

const char CompressedSequence::operator[](size_t index) const {
  size_t i = index / 4;
  size_t j = index % 4;
  size_t idx = ((_data[i]) >> (2*j)) & 0x03;
  return bases[idx];
}

//void CompressedSequence::setSequence(const CompressedSequence &o, size_t length, size_t offset, bool reversed) {
//  setSequence(o,0,length,offset,reversed);
//}


// use:  a.setSequence(b, start, length, offset, reversed)
// pre:  start+length <= o._length, offset <= a._length
// post: copies compressed sequence from b to a (reverse complement if reversed == true)
//       the string copied from b is from [start,...,start+length-1]
//       the positions in a that are update are [offset,...,offset+length-1]
//       capacity of a might be updated to fit the new string.
void CompressedSequence::setSequence(const CompressedSequence &o, size_t start, size_t length, size_t offset, bool reversed) {
  assert( o._length + start <= length);
  if (round_to_bytes(length+offset) > _capacity) {
    _resize_and_copy(round_to_bytes(length+offset), offset);
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
  _length = offset+length;
}



void CompressedSequence::reserveLength(size_t new_length) {
  if (round_to_bytes(new_length) > _capacity) {
    _resize_and_copy(round_to_bytes(new_length), _length);
  }
}


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

void CompressedSequence::setSequence(const char *s, size_t length, size_t offset, bool reversed) {
  if(round_to_bytes(length+offset) > _capacity) {
    _resize_and_copy(round_to_bytes(length+offset), offset);
  }

  for (size_t index = offset; index < offset+length; index++) {
    size_t i = index / 4;
    size_t j = index % 4;
    _data[i] &= ~(0x03 << (2*j)); // set bits to 0, default
    uint8_t c = reversed ? bases[0x03-bits[(uint8_t)*(s+length+offset-index-1)]] : *(s+index);
    _data[i] |= (bits[c] << (2*j));
  }
  _length = offset+length;
}

void CompressedSequence::setSequence(const string &s, size_t length, size_t offset, bool reversed) {
  setSequence(s.c_str(),length,offset,reversed);
}

void CompressedSequence::setSequence(const Kmer &km, size_t length, size_t offset, bool reversed) {
  char s[Kmer::MAX_K+1];
  km.toString(&s[0]);
  setSequence(s,length, offset, reversed);
}

string CompressedSequence::toString() const {
  return toString(0,_length);
}

string CompressedSequence::toString(size_t offset, size_t length) const {
  assert(offset+length <= _length);
  string s(length,0);
  size_t i,j,idx;
  for (size_t index = offset; index < offset+length; index++) {
    i = index / 4;
    j = index % 4;
    idx = ((_data[i]) >> (2*j)) & 0x03;
    s[index] = bases[idx];
  }
  return s;
}

void CompressedSequence::toString(char *s) const {
  toString(s,0,_length);
}

void CompressedSequence::toString(char *s, size_t offset, size_t length) const {
  assert(offset+length <= _length);
  size_t i,j,idx;
  for (size_t index = offset; index < offset+length; index++) {
    i = index / 4;
    j = index % 4;
    idx = ((_data[i]) >> (2*j)) & 0x03;
    s[index] = bases[idx];
  }
  s[length] = 0; // 0-terminated string
}

Kmer CompressedSequence::getKmer(size_t offset) const {
  char s[Kmer::MAX_K+1];	
  toString(&s[0]);
  return Kmer(s);
}

CompressedSequence CompressedSequence::rev() const {
	CompressedSequence r;
	r.setSequence(*this, 0, _length, 0, true);
	return r;
}
