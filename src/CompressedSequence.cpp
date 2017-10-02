#include <iostream>

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


CompressedSequence::CompressedSequence() {
  initShort();
}

void CompressedSequence::initShort() {
  asBits._size = 1; // short and size 0
  memset(&asBits._arr[0],0,31); // clear other bits
}

// use:  delete c;
// pre:  c is a pointer to a CompressedSequence
// post: the memory which the CompressedSequence had allocated has been freed
CompressedSequence::~CompressedSequence() {
  clear();
}
// use:  _cs = CompressedSequence(cs);
// pre:
// post: the DNA string in _cs and is the same as in cs
CompressedSequence::CompressedSequence(const CompressedSequence& o) {
  if (o.isShort()) {
    asBits._size = o.asBits._size;
    memcpy(asBits._arr, o.asBits._arr, 31);
  } else {
    setSequence(o,0,o.size()); // copy sequence and pointers etc.
  }
}


// use:  _cs = cs;
// pre:
// post: the DNA string in _cs is the same as in cs
CompressedSequence& CompressedSequence::operator=(const CompressedSequence& o) {
  if (o.isShort()) {
    asBits._size = o.asBits._size;
    memcpy(asBits._arr, o.asBits._arr,31); // plain vanilla copy
  } else {
    setSequence(o,0,o.size()); // copy sequence and pointers etc.
  }
  return *this;
}


// use:  s
// pre:  s has only the characters 'A','C','G' and 'T' and can have any length
// post: the DNA string in cs is now the same as s
CompressedSequence::CompressedSequence(const char *s) {
  initShort();
  if (s != NULL) {
    setSequence(s,strlen(s));
  }
}


// same as with char *s but with string
CompressedSequence::CompressedSequence(const string& s) {
  initShort();
  setSequence(s.c_str(), s.size());
}


// use:  cs = CompressedSequence(km);
// pre:
// post: the DNA string in cs is now the same as the DNA string in km
CompressedSequence::CompressedSequence(const Kmer& km) {
  initShort();
  setSequence(km, Kmer::k);
}


// use:  c = cs[index];
// pre:  0 <= index < cs.size()
// post: c is character nr. index in the DNA string inside cs
const char CompressedSequence::operator[](size_t index) const {
  const char *data = getPointer();
  size_t i = index / 4;
  size_t j = index % 4;
  size_t idx = ((data[i]) >> (2*j)) & 0x03;
  return bases[idx];
}


//void CompressedSequence::setSequence(const CompressedSequence &o, size_t length, size_t offset, bool reversed) {
//  setSequence(o,0,length,offset,reversed);
//}

bool CompressedSequence::isShort() const {
  return ((asBits._size & shortMask) == 1);
}

const char *CompressedSequence::getPointer() const {
  if (isShort()) {
    return &(asBits._arr[0]);
  } else {
    return asPointer._data;
  }
}

size_t CompressedSequence::capacity() const {
  if (isShort()) {
    return 31; // 31 bytes
  } else {
    return asPointer._capacity;
  }
}

void CompressedSequence::setSize(size_t size) {
  if (isShort()) {
    asBits._size = ((0x7F & size) << 1) | 1; // set short flag
  } else {
    asPointer._length = (0x7FFFFFFF & size) << 1;
  }
}

size_t CompressedSequence::size() const {
  if (isShort()) {
    return (asBits._size >> 1);
  } else {
    return (asPointer._length >> 1);
  }
}

// use:  a.setSequence(b, start, length, offset, reversed);
// pre:  start+length <= b._length, offset <= a._length
// post: copies compressed sequence from b to a (reverse complement if reversed == true)
//       the string copied from b is from [start,...,start+length-1]
//          (reverse complement of [o._length-1-start-length,...,o._length-1-start] if reversed == true)
//       the positions in a that are updated are [offset,...,offset+length-1]
//       capacity of a might be updated to fit the new string.
void CompressedSequence::setSequence(const CompressedSequence& o, size_t start, size_t length, size_t offset, bool reversed) {
  assert(length + start <= o.size());

  if (round_to_bytes(length+offset) > capacity()) {
    _resize_and_copy(round_to_bytes(length+offset),size());
  }

  char *data = const_cast<char *>(getPointer());
  const char *odata = o.getPointer();

  size_t w_index = offset;
  size_t r_index = reversed ? o.size()-start-1 : start;
  size_t wi,wj,ri,rj;

  for (size_t i = 0; i < length; i++) {
    wi = w_index / 4;
    wj = w_index % 4;
    ri = r_index / 4;
    rj = r_index % 4;
    data[wi] &= ~(0x03<<(2*wj)); // clear bits
    uint8_t nucl = (odata[ri] >> (2*rj)) & 0x03; // nucleotide stored in o
    if (reversed) {
      nucl = 3-nucl; // reverse sequence
    }
    data[wi] |= (nucl << (2*wj));

    w_index++;
    if (reversed) {
      r_index--;
    } else {
      r_index++;
    }
  }
  // new length?
  if (offset + length > size()) {
    setSize(offset+length);
  }
}



// use:  cs.reserveLength(new_length);
// pre:
// post: The DNA string in cs has space for at least new_length bases
void CompressedSequence::reserveLength(size_t new_length) {
  if (round_to_bytes(new_length) > capacity() ) {
    _resize_and_copy(round_to_bytes(new_length), size());
  }
}


// use:  cs._resize_and_copy(new_cap, copy_limit);
// pre:
// post: The DNA string in cs has space for at least new_length bases
//       the first copy_limit characters of cs are the same as before this method
void CompressedSequence::_resize_and_copy(size_t new_cap, size_t copy_limit) {
  if (new_cap <= capacity()) {
    return;
  }

  char *new_data = new char[new_cap]; // allocate new storage
  size_t bytes = round_to_bytes(copy_limit);
  memcpy(new_data, getPointer(), bytes); // copy old data

  if (isShort()) {
    size_t sz = size();
    asBits._size = 0; // this is now a long sequence.
    setSize(sz);
    asPointer._data = new_data;
    asPointer._capacity = new_cap;
  } else {
    delete[] asPointer._data;
    asPointer._data = new_data;
    asPointer._capacity = new_cap;
  }
}


// use:  a.setSequence(s, length, offset, reversed);
// pre:  length <= strlen(s), offset <= a._length
// post: copies substring s[0,...,length-1] to a (reverse complement if reversed == true)
//       the positions in a that are updated are [offset,...,offset+length-1]
//       capacity of a might be updated to fit the new string.
void CompressedSequence::setSequence(const char *s, size_t length, size_t offset, bool reversed) {
  if(round_to_bytes(length+offset) > capacity()) {
    _resize_and_copy(round_to_bytes(length+offset), size());
  }
  char *data = const_cast<char *>(getPointer());


  for (size_t index = offset; index < offset+length; index++) {
    size_t i = index / 4;
    size_t j = index % 4;
    data[i] &= ~(0x03 << (2*j)); // set bits to 0, default
    uint8_t c = reversed ? bases[0x03-bits[(uint8_t)*(s+length+offset-index-1)]] : *(s+index-offset);
    data[i] |= (bits[c] << (2*j));
  }

  if (offset + length > size()) {
    setSize(offset + length);
  }
}


// use:  cs.setSequence(s, start, length, offset, reversed);
// pre:  0 <= start + length < o.size()
// post: If reversed is false then: cs[offset,...,offset+length-1] = s[0,...,start+length-1]
//       else: cs[offset,...,offset+length-1] is the reverse complement of s[0,...,start+length-1]
void CompressedSequence::setSequence(const string& s, size_t length, size_t offset, bool reversed) {
  setSequence(s.c_str(),length,offset,reversed);
}


// use:  cs.setSequence(km, length, offset, reversed);
// pre:  0 <= length < cs._length,
//       length <= Kmer::k
// post: If reversed is false then: cs[offset,...,offset+length-1]
//         is the first length characters from the DNA string in km
//       else: cs[offset,...,offset+length-1] is the first length characters from the
//         reverse complement of the DNA string in km
void CompressedSequence::setSequence(const Kmer& km, size_t length, size_t offset, bool reversed) {
  char s[Kmer::MAX_K+1];
  km.toString(&s[0]);
  setSequence(s, length, offset, reversed);
}


// use:  s = cs.toString(offset, length);
// pre:  offset + length <= cs.size(),
// post: s is the DNA string from c[offset,...,offset+length-1]
string CompressedSequence::toString(const size_t offset, const size_t length) const {

    const char *data = getPointer();

    assert(offset+length <= size());

    string s(length, 0);

    size_t i,j,idx;

    for (size_t index = offset; index < offset+length; index++) {

        i = index / 4;
        j = index % 4;
        idx = ((data[i]) >> (2*j)) & 0x03;
        s[index-offset] = bases[idx];
    }

    return s;

    /*char v, tmp = data[offset / 4] >> (2 * (offset % 4));

    for (size_t i = offset; i < offset + length; i++, tmp >>= 2) {

        if (i % 4 == 0) tmp = data[i / 4];

        v = tmp & 0x3;
        s[i-offset] = 0x40 | (v + 1) | (0x1 << (v-1)*2);
    }

    return s;*/
}


// use:  s = cs.toString(offset, length);
// pre:  offset + length <= cs.size()
//       s has space for length characters
// post: s is the same as cs[offset,...,offset+length-1]
void CompressedSequence::toString(char *s, const size_t offset, const size_t length) const {

    const char *data = getPointer();

    assert(offset+length <= size());

    size_t i,j,idx;

    for (size_t index = offset; index < offset+length; index++) {

        i = index / 4;
        j = index % 4;
        idx = ((data[i]) >> (2*j)) & 0x03;
        s[index-offset] = bases[idx];
    }

    s[length] = 0; // 0-terminated string

    /*char v, tmp = data[offset / 4] >> (2 * (offset % 4));

    for (size_t i = offset; i < offset + length; i++, tmp >>= 2, s++) {

        if (i % 4 == 0) tmp = data[i / 4];

        v = tmp & 0x3;
        *s = 0x40 | (v + 1) | (0x1 << (v-1)*2);
    }

    *s = '\0';*/
}

Kmer CompressedSequence::getKmer(const size_t offset) const {

    /*Kmer km;

    const char* data = getPointer();

    const size_t len = offset + Kmer::k;

    size_t j = 0;

    uint64_t tmp = (uint64_t)(data[offset / 4] >> (2 * (offset % 4)));

    for (size_t i = offset; i != len; i++, j++, tmp >>= 2) {

        if (i % 4 == 0) tmp = (uint64_t)(data[i / 4]);

        km.longs[j/32] = (km.longs[j/32] << 2) | (tmp & 0x3);
    }

    km.longs[j/32] <<= 2 * (32 - (j % 32));

    return km;*/

    Kmer km;

    const char* data = getPointer();

    const size_t len = offset + Kmer::k;
    const size_t nlongs = (Kmer::k + 31) / 32;

    size_t j = 0;

    uint64_t tmp_data = (uint64_t)(data[offset / 4] >> (2 * (offset % 4)));

    for (size_t i = offset; j < nlongs; j++){

        uint64_t tmp_km = 0;
        const size_t end = len < i + 32 ? len : i + 32;

        for (; i != end; i++, tmp_data >>= 2){

            if (i % 4 == 0) tmp_data = (uint64_t) data[i / 4];
            tmp_km = (tmp_km << 2) | (tmp_data & 0x3);
        }

        km.longs[j] = tmp_km;
    }

    km.longs[j-1] <<= 2 * (32 - (Kmer::k % 32));

    return km;
}

bool CompressedSequence::compareKmer(const size_t offset, const Kmer& km) const {

    /*const char* data = getPointer();

    const size_t len = offset + Kmer::k;

    uint64_t tmp_km, tmp_data = (uint64_t)(data[offset / 4] >> (2 * (offset % 4)));

    for (size_t i = offset, j = 0; i != len; i++, j++, tmp_data >>= 2, tmp_km <<= 2) {

        if (i % 4 == 0) tmp_data = (uint64_t)(data[i / 4]);
        if (j % 32 == 0) tmp_km = km.longs[j / 32];

        if ((tmp_data & 0x3) != (tmp_km >> 62)) return false;
    }

    return true;*/

    const char* data = getPointer();

    const size_t len = offset + Kmer::k;
    const size_t nlongs = (Kmer::k + 31) / 32;

    uint64_t tmp_data = (uint64_t)(data[offset / 4] >> (2 * (offset % 4)));

    for (size_t i = offset, j = 0; j < nlongs; j++){

        uint64_t tmp_km = km.longs[j];
        const size_t end = len < i + 32 ? len : i + 32;

        for (; i != end; i++, tmp_km <<= 2, tmp_data >>= 2){

            if (i % 4 == 0) tmp_data = (uint64_t) data[i / 4];
            if ((tmp_data & 0x3) != (tmp_km >> 62)) return false;
        }
    }

    return true;
}

int64_t CompressedSequence::findKmer(const Kmer& km) const {

    int k = Kmer::k;
    size_t sz = size();

    if (sz >= k){

        Kmer km_cs = getKmer(0);
        if (km_cs == km) return 0;

        if (sz > k){

            size_t i = k;
            const char* data = getPointer();
            char tmp = data[i/4] >> (2 * (i % 4));

            for (; i < sz; i++, tmp >>= 2){

                if (i%4 == 0) tmp = data[i/4];
                km_cs.selfForwardBase(bases[tmp & 0x3]);

                if (km_cs == km) return i-k+1;
            }
        }
    }

    return -1;
}

// use:  _cs = cs.rev();
// pre:
// post: _cs is the reverse complement CompressedSequence with respect to cs,
//       i.e. if the DNA string in cs is 'GTCA'
//          then the DNA string in _cs is 'TGAC'
CompressedSequence CompressedSequence::rev() const {
  CompressedSequence r;
  r.setSequence(*this, 0, size(), 0, true);
  return r;
}


// use:  j = cs.jump(s,i,pos,reversed)
// pre:  0 <= i < s.length, -1 <= pos < cs._length if reversed true, 0 <= pos <= cs._length if reversed false
// post: if reversed == false
//         s[i...i+j-1] == cs._data[pos...pos+j-1], 0 <= j <= min(s.length-i, cs._length-pos)
//       else
//         reverse_complement(s[i...i+j-1]) == cs._data[pos-j+1...pos], 0 <= j <= min(s.length-i, pos+1)
size_t CompressedSequence::jump(const char *s, size_t i, int pos, bool reversed) const {

    assert(i >= 0);
    assert(i < strlen(s));
    assert(pos >= -1);
    assert(0 <= size() - pos); // this prevents -1 <= _length from giving false

    const char* data = getPointer();

    size_t i_cpy = i;

    if (reversed){

        if (pos == -1) return 0;

        int idx_div = pos / 4;
        int idx_mod = 2 * (pos % 4);

        for (; (s[i_cpy] != '\0') && (pos != -1); pos--, i_cpy++, idx_mod -= 2) {

            if (idx_mod == -2){
                idx_div--;
                idx_mod = 6;
            }

            if (s[i_cpy] != bases[3-((data[idx_div] >> idx_mod) & 0x03)]) break;
        }
    }
    else {

        size_t cs_size = size();

        if (pos == cs_size) return 0;

        char tmp = data[pos/4] >> (2 * (pos % 4));

        for (; (s[i_cpy] != '\0') && (pos != cs_size); pos++, i_cpy++, tmp >>= 2) {

            if (pos%4 == 0) tmp = data[pos/4];
            if (s[i_cpy] != bases[tmp & 0x3]) break;
        }
    }

    return i_cpy - i;
}

size_t CompressedSequence::bw_jump(const char *s, size_t i, int pos, bool reversed) const {

    assert(i >= 0);
    assert(i < strlen(s));
    assert(pos >= -1);
    assert(0 <= size() - pos); // this prevents -1 <= _length from giving false

    const char* data = getPointer();

    size_t i_cpy = i;

    if (reversed){

        size_t cs_size = size();

        if (pos == cs_size) return 0;

        char tmp = data[pos/4] >> (2 * (pos % 4));

        for (; (i_cpy != -1) && (pos != cs_size); pos++, i_cpy--, tmp >>= 2) {

            if (pos%4 == 0) tmp = data[pos/4];
            if (s[i_cpy] != bases[3 - (tmp & 0x3)]) break;
        }
    }
    else {

        if (pos == -1) return 0;

        int idx_div = pos / 4;
        int idx_mod = 2 * (pos % 4);

        for (; (i_cpy != -1) && (pos != -1); pos--, i_cpy--, idx_mod -= 2) {

            if (idx_mod == -2){
                idx_div--;
                idx_mod = 6;
            }

            if (s[i_cpy] != bases[(data[idx_div] >> idx_mod) & 0x03]) break;
        }
    }

    return i - i_cpy;
}

void CompressedSequence::clear() {
  if (!isShort()) { // release memory if needed
    if (asPointer._capacity > 0 && asPointer._data != NULL) {
      delete[] asPointer._data;
      asPointer._data = NULL;
    }
  }
}
