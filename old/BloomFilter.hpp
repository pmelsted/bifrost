#ifndef BFG_BLOOMFILTER_HPP
#define BFG_BLOOMFILTER_HPP

#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>

#include "hash.hpp"
#include "libdivide.h"
#include "stdint.h"

//using namespace libdivide;
using namespace std;

static const unsigned char mask[8] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};

/*static const unsigned char BitsSetTable256[256] = 
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
};
*/


/* Short description: 
 *  - Use the very fast MurmurHash to hash keys into a Bloom Filter 
 *  - Gain speed by using libdivide when not dividing by powers of 2
 *  - Easily write a Bloom Filter to a file
 *  - Easily read a Bloom Filter from a file
 *  - Pretty print the false positive rate
 * */
class BloomFilter {
protected:
  unsigned char* table_; // bit array
  uint32_t seed_;        // for hash functions
  uint64_t size_;        // number of bits, table_ has (size_/8) number of bytes
  size_t k_;             // number of hash functions
  libdivide::divider<uint64_t> fast_div_; // fast division
public:
  BloomFilter() : seed_(0), size_(0), table_(NULL), k_(0), fast_div_() {}

  BloomFilter(size_t num, size_t bits, uint32_t seed) : seed_(seed), size_(0), table_(NULL), fast_div_() {
    //cout << "num="<<num << ", bits="<<bits;
    size_ = rndup(bits*num);
    //cout <<", size=" << size_ << endl;

    init_table();
    init_k(bits);
  } 

  ~BloomFilter() {
    clear();
  }

  size_t memory() const { 
    size_t m = sizeof(BloomFilter) + (size_ >> 3);
    fprintf(stderr, "BloomFilter:\t\t%zuMB\n",  m >> 20);
    return m;
  }

  template<typename T>
  bool contains(T x) const {
    return (search(x) == 0);
  }

  // use:  r = bf.search(x)
  // pre:
  // post: r is the number of bits that need to be set to 1 so that
  //       x is a member of bf
  template<typename T>
  size_t search(T x) const 
  {
    size_t r = k_;
    uint64_t id;
    uint64_t hash;
    uint64_t hash0; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_  , &hash0);
    hash0 |= 1; // odd number
    uint64_t hash1; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+1, &hash1);
    for (uint64_t i = 0; i < k_; i++) {
      //MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+i, &hash);
      hash = hash0 * i + hash1;
      id = hash - (hash / fast_div_) * size_; // equal to hash % size_;
      //assert(id == (hash % size_));
      if ((table_[id >> 3] & mask[id & 0x07]) != 0) {
        r--;      
      }
    }
    return r;
  }

  // use:  r = bf.insert(x)
  // pre:
  // post: x is a member of bf, r is the number of bits modified
  template<typename T>
  size_t insert(T x) {
    size_t r = 0;
    uint64_t id;
    uint64_t hash;
    unsigned char val;
    uint64_t hash0; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_  , &hash0);
    hash0 |= 1; // odd number
    uint64_t hash1; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+1, &hash1);
    for(uint64_t i = 0; i < k_; i++) {
      //MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+i,&hash);
      hash = hash0 * i + hash1;
      id = hash - (hash / fast_div_) * size_;
      //assert(id == (hash % size_));

      
      if ((table_[id>>3] & mask[id & 0x07]) == 0) {
        val = __sync_fetch_and_or(table_ + (id>>3), mask[id & 0x07]); // val is the value prior to or-ing
        // another thread could have changed it in the meantime
        if ((val & mask[id & 0x07]) == 0) {
          r++; // we changed the value 
        } 
      }
      

      /*
        while ( 1 ) {
        val = table_[id>>3];
        if ((val & mask[id & 0x07]) == 0) {
          unsigned char nval = val | mask[id & 0x07];
          if (__sync_bool_compare_and_swap(&(table_[id>>3]), val, nval)) {
            r++;
            break;
          }
        } else {
          break;
        }
      }
      */
    }
    return r;
  }

  bool WriteBloomFilter(FILE *fp) {
    // write metadata in this order (size_, seed_, k_);
    if (fwrite(&size_, sizeof(size_), 1, fp) != 1) { cout << "size_" << endl; return false;}
    if (fwrite(&seed_, sizeof(seed_), 1, fp) != 1) { cout << "seed_" << endl; return false;}
    if (fwrite(&k_,    sizeof(k_)   , 1, fp) != 1) { cout << "k_" << endl; return false;}

    // now write actual data
    if (fwrite(table_, sizeof(unsigned char), size_>>3, fp) != (size_>>3))  { cout << "table_" << endl; return false;}

    return true;
  }

  bool ReadBloomFilter(FILE *fp) {
    clear(); // free current table
    // read metadata
    if (fread(&size_, sizeof(size_), 1, fp) != 1) return false;
    if (fread(&seed_, sizeof(seed_), 1, fp) != 1) return false;
    if (fread(&k_,    sizeof(k_),    1, fp) != 1) return false;

    // allocate memory
    init_table();
    // read table
    if (fread(table_, sizeof(unsigned char), size_>>3, fp) != (size_>>3)) return false;
    // done

    return true;
  }

  size_t count() const {
    size_t c = 0;
    for (size_t i = 0; i < (size_ >> 3); i++) {
      unsigned char u = table_[i]; 
      for (size_t j = 128; j != 0; j = j>>1) {
        if ((u & j) != 0) {
          c++;
        }
      }
    }
    cout << c << " bits set out of " << size_ << " with k = " << k_ << endl;
    if (c != 0) {
      double n = size_*(-log(1.0-((double)c)/size_))/k_;
      //cout << "estimate =" << (size_t)n  << endl;
      return (size_t) n;
    } else {
      return 0;
    }
  }


  void clear() {
    if (table_ != NULL) {
      delete[] table_;
    }
    table_ = NULL;
    size_ = 0;
  }

protected:

  void init_table() {
    fast_div_ = libdivide::divider<uint64_t>(size_);
    table_ = new unsigned char[size_>>3];
    memset(table_, 0, size_>>3);
  }

  void init_k(size_t bits) {
    size_t k = (size_t) (bits*log(2));
    if (fpp(bits,k) < fpp(bits,k+1)) {
      k_ = k;
    } else {
      k_ = k+1;
    }
    cerr << "k="<<k_<<", fpp="<<fpp(bits,k_) << endl;
  }

  double fpp(size_t bits, size_t k) const {
    //    cout << bits<<","<<k<<","<<(-((double)k)/((double)bits)) << endl;
    return pow(1-exp(-((double)k)/((double)bits)),(double)k);
  }
  
  uint64_t rndup(uint64_t x) const {
    return ((x+63) >> 6)<<6;
  }
   
};

#endif // BFG_BLOOMFILTER_HPP
