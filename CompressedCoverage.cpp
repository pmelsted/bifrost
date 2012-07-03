#include <string>
#include <cstdio>
#include <stdint.h>
#include <assert.h>
#include <sstream>
#include "CompressedCoverage.hpp"

using namespace std;

size_t round_to_bytes(const size_t len)  { return (len+3)/4; }
 

CompressedCoverage::CompressedCoverage(size_t sz, bool full) {
  if (sz > 0) {
    initialize(sz, full);
  } 
}


void CompressedCoverage::initialize(size_t sz, bool full) {
  if (sz <= size_limit) {
    asBits = intptr_t(0); // zero out
    asBits |= tagMask;  // set 0-bit to 1
    asBits |= (sizeMask & (sz << 2)); // set bits 2-6 to size;
    if (full) {
      asBits |= fullMask;
    }
  } else {
    if (!full) {
      uint8_t* ptr = new uint8_t[8+round_to_bytes(sz)];
      *(reinterpret_cast<uint32_t*>(ptr)) = (uint32_t) sz; // first 4 bytes store size
      *(reinterpret_cast<uint32_t*>(ptr+4)) = (uint32_t) sz; // next  4 bytes store number of uncovered bases
      memset(ptr+8, 0, round_to_bytes(sz)); // 0 out array allocated
      asPointer = ptr; // last bit is 0
    } else {
      asBits = fullMask;
      asBits |= (sz << 32);
    }
  }
  assert(sz == size());
}


CompressedCoverage::~CompressedCoverage() {
  releasePointer();
}


void CompressedCoverage::releasePointer() {
  if ((asBits & tagMask) == 0 && (asBits & fullMask) != fullMask) {
    // release pointer
    uint8_t* ptr = getPointer();
    size_t sz = size();
    intptr_t *change = &asBits;
    intptr_t oldval = *change;
    intptr_t newval = fullMask | (sz << 32);
    while ((oldval & fullMask) != fullMask) {
      if (__sync_bool_compare_and_swap(change, oldval, newval)) {
        delete[] ptr;
        break;
      }
      oldval = *change;
    }
    //asBits = fullMask;
    //asBits |= (sz << 32);
    //delete[] ptr;
  }
}


size_t CompressedCoverage::size() const {
  if ((asBits & tagMask) == tagMask) {
    return ((asBits & sizeMask) >> 2);
  } else {
    if ((asBits & fullMask) == fullMask) {
      return (uint32_t) (asBits >> 32);
    } else {
      return *((uint32_t*) getPointer());
    }
  }
}



uint8_t* CompressedCoverage::getPointer() const {
  assert ((asBits & tagMask) == 0);
  return reinterpret_cast<uint8_t*>(asBits & pointerMask);  
}


string CompressedCoverage::toString() const {
  bool isPtr = ((asBits & tagMask) == 0);
  size_t sz = size();
  bool full = isFull();

  string bits(64, '0');
  
  for (int i = 0; i < 64; i++) {
    if (asBits & (intptr_t(1) << (63-i))) {
      bits[i] = '1';
    }
  }
  
  if (isPtr) {
    ostringstream info;
    info << "Pointer: ";
    if (full) {
      info << "Full, size = ";
      info << sz;
      info << endl;
    } else {
      info << "Non-full, size = " << sz << ", not-filled = ";
      uint32_t filled = *((const uint32_t*)(getPointer()+4));
      info << filled << endl;
              
      size_t nbytes = round_to_bytes(sz);
      uint8_t *ptr = getPointer() + 8;
      string ptrbits(nbytes*8, '0');
      for (size_t i = 0; i < nbytes; i++) {
        for (size_t j = 0; j < 8; j++) {
          if (((ptr[i] & (1<<j)) >> j) == 1) {
            ptrbits[nbytes*8-1-(8*i+j)] = '1';
          }
        }
      }
      //info << ptrbits;
      bits += ptrbits;
    }
    return bits + "\n" + info.str();
  } else {
    ostringstream info;
    info << "Local array:";
    if (full) {
      info << ", Full,";
    }
    info <<" size = " << sz;
    return bits + "\n" + info.str();
  }
}


void CompressedCoverage::cover(size_t start, size_t end) {
  if (end >= size()) {
    printf("start=%zu, end=%zu\n", start, end);
    printf("size=%zu\n", size());
  }
  assert(end < size());

  if (isFull()) {
    return;
  } else {
    if ((asBits & tagMask) == tagMask) { // local array
      intptr_t s = intptr_t(3); // 0b11
      s <<= (8 + 2*start); 
      size_t val;
      for (; start <= end; start++,s<<=2) {
        val = (asBits & s) >> (8 + 2*start);
        if (val < 2) {
          val++;
          intptr_t *change = &asBits;
          //asBits &= ~s; // clear bits
          //asBits |= (val << (8 +  2*start));
          while (1) {
            intptr_t oldval = *change;
            intptr_t newval = oldval & ~s;
            newval |= ( val << (8 + 2*start));
            if (__sync_bool_compare_and_swap(change, oldval, newval)) {
              break;
            }
          }
        }
      }
      if (isFull()) {
        asBits |= fullMask;
      }
    } else {
      size_t fillednow = 0;
      uint8_t s = intptr_t(3); // 0b11
      uint8_t* ptr = getPointer() + 8;
      size_t val;
      size_t index, pos;
      for (; start <= end; start++) {
        index = start >> 2;
        pos = 2*(start & 0x03);
        val = (ptr[index] & (s << pos)) >> pos;
        if (val < 2) {
          val++;
          uint8_t *change = &ptr[index];
          while (1) {
            intptr_t oldval = *change;
            if ((oldval & (s << pos)) >> pos == 2) {
              break;
            }
            intptr_t newval = oldval & ~(s << pos);
            newval |= ( val << (pos));
            if (__sync_bool_compare_and_swap(change, oldval, newval)) {
              if (((newval & (s << pos)) >> pos)  == 2) {
                fillednow++;
              }
              break;
            }
          }
          //ptr[index] &= ~(s << pos);
          //ptr[index] |= (val << pos);
        }
      }
      uint32_t *change = (uint32_t *) (getPointer() + 4);
      while (1) {
        uint32_t oldval = *change;
        uint32_t newval = oldval - fillednow;
        assert(newval >= 0);
        if (__sync_bool_compare_and_swap(change, oldval, newval)) {
          break;
        }
      }
      
      if (isFull()) {
        releasePointer();
        assert((asBits & fullMask) == 2);
      }
    }
  }
}


uint8_t CompressedCoverage::covAt(size_t index) const {
  assert((asBits & fullMask) != fullMask);
  if ((asBits & tagMask) == tagMask) { 
    return ((static_cast<uintptr_t>(asBits) >> (8 + 2*index)) & 0x03);
  } else {
    uint8_t* ptr = getPointer() + 8;
    size_t pos = 2*(index & 0x03);
    return (ptr[index >> 2] & (0x03 << pos)) >> pos;
  }
}


size_t CompressedCoverage::lowCoverageCount() const {
  if (isFull()) {
    return 0;
  }
  size_t sz = size();
  size_t low = 0;
  for(size_t i=0; i<sz; ++i) {
    if (covAt(i) < 2) {
      ++low;
    }
  }
  return low;
}


vector<pair<int, int> > CompressedCoverage::getSplittingVector() const {
  size_t a = 0, b = 0, sz = size();
  vector<pair<int, int> > v;
    
  while (b != sz) {
    // [a,...,b-1] is a fully covered subinterval and (a,b) has been added to v
    while (a < sz && covAt(a) <= 1) {
      a++;
    }
    if (a == sz) {
      break;
    }
    b = a;
    while (b < sz && covAt(b) > 1) {
      b++; 
    }
    v.push_back(make_pair(a,b));
    a = b;
  }
  return v;
}

bool CompressedCoverage::isFull() const {
  if ((asBits & fullMask) == fullMask) {
    return true;
  }
  
  if ((asBits & tagMask) == tagMask) {
    size_t sz = size();
    return (static_cast<uintptr_t>(asBits) >> 8) == (0xAAAAAAAAAAAAAA >> 2*(28 - sz));
  } else {
    size_t uncovered = *((const uint32_t*) (getPointer()+4));
    return (uncovered == 0);
  }
}
