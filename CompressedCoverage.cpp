#include <string>
#include <cstdio>
#include <stdint.h>
#include <assert.h>
#include <sstream>
#include "CompressedCoverage.hpp"

using namespace std;

size_t round_to_bytes(const size_t len)  { return (len+3)/4; }


CompressedCoverage::CompressedCoverage(size_t _size) {
  if (_size <= size_limit) {
    asBits = intptr_t(0); // zero out
    asBits |= tagMask;  // set 0-bit to 1
    asBits |= (sizeMask & (_size << 2)); // set bits 2-6 to size;
  } else {
    uint8_t* ptr = new uint8_t[8+round_to_bytes(_size)];
    *(reinterpret_cast<uint32_t*>(ptr)) = (uint32_t) _size; // first 4 bytes store size
    *(reinterpret_cast<uint32_t*>(ptr+4)) = (uint32_t) _size; // next  4 bytes store number of uncovered bases
    memset(ptr+8, 0, round_to_bytes(_size)); // 0 out array allocated
    asPointer = ptr; // last bit is 0
    assert(getPointer() == ptr);
    assert(size() == _size);
  }
}


CompressedCoverage::~CompressedCoverage() {
  releasePointer();
}


void CompressedCoverage::releasePointer() {
  if ((asBits & tagMask) == 0 && (asBits & fullMask) != 2) {
    // release pointer
    uint8_t* ptr = getPointer();
    size_t sz = size();
    asBits = fullMask;
    asBits |= (sz << 32);
    delete[] ptr;
  }
}


size_t CompressedCoverage::size() const {
  if ((asBits & tagMask) == 1) {
    return ((asBits & sizeMask) >> 2);
  } else {
    if ((asBits & fullMask) == 2) {
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
  assert(start <= end);
  assert(end < size());
  if (isFull()) {
    return;
  } else {
    if ((asBits & tagMask) == 1) { // local array
      intptr_t s = intptr_t(3); // 0b11
      s <<= (8 + 2*start); 
      size_t val;
      for (; start <= end; start++,s<<=2) {
        val = (asBits & s) >> (8 + 2*start);
        if (val < 2) {
          val++;
          asBits &= ~s; // clear bits
          asBits |= (val << (8 +  2*start));
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
          ptr[index] &= ~(s << pos);
          ptr[index] |= (val << pos);
          if (val == 2) {
            fillednow++;
          }
        }
      }
      *((uint32_t*) (getPointer()+4)) -= fillednow;
      

      if (isFull()) {
        releasePointer();
        assert((asBits & fullMask) == 2);
      }
      
    }
  }
}


bool CompressedCoverage::isFull() const {
  if ((asBits & fullMask) == 2) {
    return true;
  }
  
  if ((asBits & tagMask) == 1) {
    size_t sz = size();
    intptr_t b = (((1 << sz)-1) << 8);
    return ((asBits & b) == (0xAAAAAAAAAAAAAA00 & b));
  } else {
    size_t uncovered = *((const uint32_t*) (getPointer()+4));
    return (uncovered == 0);
  }
}
