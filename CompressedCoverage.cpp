#include <string>
#include <cstdio>
#include <stdint.h>
#include <assert.h>
#include <sstream>
#include "CompressedCoverage.hpp"

using namespace std;

size_t round_to_bytes(const size_t len)  { return (len+3)/4; }
 

// use:  cc = CompressedCoverage(sz, full);
// post: if (sz > 0) then initialize the instance else skip initializing
CompressedCoverage::CompressedCoverage(size_t sz, bool full) {
  if (sz > 0) {
    initialize(sz, full);
  } else {
    asBits = tagMask;
  }
}


// use:  cc.initialize(sz, full);
// post: the data structure has been initialized either as a small local array on the stack
//       or a bigger array on the heap if sz > size_limit
void CompressedCoverage::initialize(size_t sz, bool full) {
  if (sz <= size_limit) {
    asBits = 0; // zero out
    asBits |= tagMask;  // set 0-bit to 1
    asBits |= (sizeMask & (sz << 2)); // set bits 2-6 to size;
    if (full) {
      asBits |= fullMask;
    }
  } else {
    if (!full) {
      uint8_t* ptr = new uint8_t[8+round_to_bytes(sz)];
      asPointer = ptr; // last bit is 0
      *(get32Pointer()) = sz;
      *(get32Pointer() + 1) = sz;
      memset(ptr+8, 0, round_to_bytes(sz)); // 0 out array allocated
    } else {
      asBits = fullMask;
      asBits |= (sz << 32);
    }
  }
  assert(sz == size());
}


// use:  delete cc; 
// post: 
CompressedCoverage::~CompressedCoverage() {
}


// use:  cc.releasePointer();  
// post: if there was data on the heap then it has been freed
void CompressedCoverage::releasePointer() {
  if ((asBits & tagMask) == 0 && (asBits & fullMask) != fullMask) {
    // release pointer
    uint8_t* ptr = get8Pointer();
    uintptr_t sz = size();
    uintptr_t *change = &asBits;
    uintptr_t oldval = *change;
    uintptr_t newval = fullMask | (sz << 32);
    while ((oldval & fullMask) != fullMask) {
      if (__sync_bool_compare_and_swap(change, oldval, newval)) {
        delete[] ptr;
        break;
      }
      oldval = *change;
    }
  }
}


// use:  i = cc.size();
// post: i is the number of kmers that cc can hold coverage for
size_t CompressedCoverage::size() const {
  if ((asBits & tagMask) == tagMask) {
    return ((asBits & sizeMask) >> 2);
  } else {
    if ((asBits & fullMask) == fullMask) {
      return asBits >> 32;
    } else {
      return *(get32Pointer());
    }
  }
}


// use:  s = cc.toString();
// post: s contains all important information about cc
string CompressedCoverage::toString() const {
  bool isPtr = ((asBits & tagMask) == 0);
  size_t sz = size();
  bool full = isFull();
  uintptr_t one(1);

  string bits(64, '0');
  
  for (int i = 0; i < 64; i++) {
    if (asBits & (one << (63-i))) {
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
      const uint32_t filled = *(getConst32Pointer() + 1);
      info << filled << endl;
              
      size_t nbytes = round_to_bytes(sz);
      uint8_t *ptr = get8Pointer() + 8;
      string ptrbits(nbytes*8, '0');
      for (size_t i = 0; i < nbytes; i++) {
        for (size_t j = 0; j < 8; j++) {
          if (((ptr[i] & (1<<j)) >> j) == 1) {
            ptrbits[nbytes*8-1-(8*i+j)] = '1';
          }
        }
      }
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


// use:  cc.cover(start, end);
// pre:  start <= end, end < cc.size()
// post: the coverage of kmers: start,...,end has been increased by one
void CompressedCoverage::cover(size_t start, size_t end) {
  assert(end < size());

  if (isFull()) {
    return;
  } else {
    if ((asBits & tagMask) == tagMask) { // local array
      uintptr_t s(3); // 0b11
      s <<= (8 + 2*start); 
      size_t val;
      for (; start <= end; start++,s<<=2) {
        val = (asBits & s) >> (8 + 2*start);
        if (val < 2) {
          val++;
          uintptr_t *change = &asBits;
          while (1) {
            uintptr_t oldval = *change;
            uintptr_t newval = oldval & ~s;
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
      uint8_t s(3); // 0b11
      uint8_t* ptr = get8Pointer() + 8;
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
            uintptr_t oldval = *change;
            if ((oldval & (s << pos)) >> pos == 2) {
              break;
            }
            uintptr_t newval = oldval & ~(s << pos);
            newval |= ( val << pos);
            if (__sync_bool_compare_and_swap(change, oldval, newval)) {
              if (((newval & (s << pos)) >> pos)  == 2) {
                fillednow++;
              }
              break;
            }
          }
        }
      }
      if (fillednow > 0) {
        // Decrease filledcounter
        __sync_add_and_fetch(get32Pointer() + 1, -fillednow);
      }
      /*
      uint32_t *change = get32Pointer() + 1;
      while (1) {
        uint32_t oldval = *change;
        uint32_t newval = oldval - fillednow;
        if (__sync_bool_compare_and_swap(change, oldval, newval)) {
          break;
        }
      }
      */
      
      if (isFull()) {
        releasePointer();
        assert((asBits & fullMask) == fullMask);
      }
    }
  }
}


// use:  k = cc.covat(index);
// pre:  0 <= index < size(), cc is not full
// post: k is the coverage at index
uint8_t CompressedCoverage::covAt(size_t index) const {
  assert((asBits & fullMask) != fullMask);
  if ((asBits & tagMask) == tagMask) { 
    return ((asBits >> (8 + 2*index)) & 0x03);
  } else {
    uint8_t* ptr = get8Pointer() + 8;
    size_t pos = 2*(index & 0x03);
    return (ptr[index >> 2] & (0x03 << pos)) >> pos;
  }
}


// use:  (low, sum) = ccov.lowCoverage();
// pre:  
// post: low is the number of kmers under coverage limtis
//       sum is the sum of these low coverages
pair<size_t, size_t> CompressedCoverage::lowCoverageInfo() const {

  if (isFull()) {
    return make_pair(0,0);
  }
  size_t sz = size();
  size_t low = 0;
  size_t sum = 0;
  for(size_t i=0; i<sz; ++i) {
    size_t cov = covAt(i);
    if (cov < 2) {
      ++low;
      sum += cov;
    }
  }
  return make_pair(low, sum);
}

 
// use:  v = ccov.splittingVector();
// pre:  ccov.isFull() == false
// post: v is a vector of pairs (a1,b1), (a2,b2),...,(an,bn)
//       where ai < aj if i < j 
//         and bi < bj if i < j
//       these pairs are all the fully covered subintervals of the corresponding contig
//       i.e. [ai,...,bi-1] is fully covered
vector<pair<int, int> > CompressedCoverage::splittingVector() const {
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


// use:  b = cc.isFull();
// post: (b == true) <==> cc is full
bool CompressedCoverage::isFull() const {
  if ((asBits & fullMask) == fullMask) {
    return true;
  }
  
  if ((asBits & tagMask) == tagMask) {
    return (asBits >> 8) == (localCoverageMask >> 2*(28 - size()));
  } else {
    return *(getConst32Pointer() + 1) == 0;
  }
}

// use: cc.setFull()
// pre:
// post: cc is full and any memory is released
void CompressedCoverage::setFull() {
  if (!isFull()) {
    releasePointer();   
  }
}
