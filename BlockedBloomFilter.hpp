#ifndef BFG_BLOCKEDBLOOMFILTER_HPP
#define BFG_BLOCKEDBLOOMFILTER_HPP

#include <cmath>
#include <iostream>

#include "hash.hpp"
#include "libdivide.h"



static const uint64_t mask[8] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};


/* Short description:
 *  - Extended BloomFilter which hashes into 64-bit blocks
 *    that can be accessed very fast from the CPU cache
 * */
class BlockedBloomFilter {
 private:
  uint64_t *table_;
  uint64_t blocks_;
  uint32_t seed_;
  uint64_t size_;
  size_t k_;
  libdivide::divider<uint64_t> fast_div_; // fast division

 public:
  BlockedBloomFilter() : seed_(0), size_(0), table_(NULL), k_(0), blocks_(0), fast_div_() {}
  BlockedBloomFilter(size_t num, size_t bits, uint32_t seed) : seed_(seed), size_(0), table_(NULL), fast_div_() {
    //cout << "num="<<num << ", bits="<<bits;
    size_ = rndup512(bits*num);
    blocks_ = size_/512;
    //cout <<", size=" << size_ << endl;

    init_table();
    init_k(bits);
  }

  ~BlockedBloomFilter() {
    clear();
  }


  size_t memory() const {
    size_t m = sizeof(BlockedBloomFilter) + (blocks_ / 64 );
    fprintf(stderr, "BlockedBloomFilter:\t\t%zuMB\n",  m >> 20);
    return m;
  }



  template<typename T>
  bool contains(T x) const {
    return (search(x) == 0);
  }


  template<typename T>
  size_t search(T x) const {
    size_t r = k_;
    uint64_t block; MurmurHash3_x64_64((const void *) &x, sizeof(T), seed_+2, &block);
    // block is the index of the 512 bit memory block where x would be stored
    // 0 <= block < blocks
    block = block - (block / fast_div_) * (blocks_); // block % blocks

    uint64_t hash0; MurmurHash3_x64_64((const void *) &x, sizeof(T), seed_ , &hash0);
    hash0 |= 1; // make hash0 an odd number
    uint64_t hash1; MurmurHash3_x64_64((const void *) &x, sizeof(T), seed_+1, &hash1);
    for (uint64_t i = 0; i < k_; i++) {
      // 0 <= bit < 512, which bit to set
      uint64_t bit = (hash0 * i + hash1) & 0x1ffULL; // equal to hash % 512;
      // we set bit number (id % 8) in byte (table_[block + id/8]) to 1
      uint64_t maskcheck = 1ULL << (bit & 0x3fULL);
      uint64_t loc = 8*block + (bit>>6);

      if ((table_[loc] &  maskcheck) != 0) {
        r--;
      }
    }
    return r;
  }

  template<typename T>
  size_t insert(T x) {
    size_t r = 0;
    uint64_t block; MurmurHash3_x64_64((const void *) &x, sizeof(T), seed_+2, &block);
    // block is the index of the 512 bit memory block where x would be stored
    // 0 <= block < blocks
    block = block - (block / fast_div_) * (blocks_); // block % blocks
    // Multiply block by 64 to get the first byte of the block

    uint64_t hash0; MurmurHash3_x64_64((const void *) &x, sizeof(T), seed_  , &hash0);
    hash0 |= 1; // make hash0 an odd number
    uint64_t hash1; MurmurHash3_x64_64((const void *) &x, sizeof(T), seed_+1, &hash1);
    for(uint64_t i = 0; i < k_; i++) {
      // 0 <= bit < 512, which bit to set
      uint64_t bit = (hash0 * i + hash1) & 0x1ffULL; // equal to hash % 512;
      // we set bit number (id % 8) in byte (table_[block + id/8]) to 1
      uint64_t maskcheck = 1ULL << (bit & 0x3fULL);
      uint64_t loc = 8*block + (bit>>6);

      if ((table_[loc] &  maskcheck) == 0) {
        uint64_t val = __sync_fetch_and_or(table_ + loc, maskcheck);
        if ((val & maskcheck) == 0) {
          r++;
        }
      }
    }
    return r;
  }


  bool WriteBloomFilter(FILE *fp) {
    if (fwrite(&size_,   sizeof(size_),   1, fp) != 1) { return false;}
    if (fwrite(&blocks_, sizeof(blocks_), 1, fp) != 1) {return false;}
    if (fwrite(&seed_,   sizeof(seed_),   1, fp) != 1) {return false;}
    if (fwrite(&k_,      sizeof(k_),      1, fp) != 1) {return false;}

    if (fwrite(table_, sizeof(uint64_t), 8*blocks_, fp) != (8*blocks_)) {return false;}
    return true;
  }

  bool ReadBloomFilter(FILE *fp) {
    clear();
    if (fread(&size_, sizeof(size_), 1, fp) != 1) { return false;}
    if (fread(&blocks_, sizeof(blocks_), 1, fp) != 1) { return false;}
    if (fread(&seed_, sizeof(seed_), 1, fp) != 1) { return false;}
    if (fread(&k_,    sizeof(k_),    1, fp) != 1) {return false;}

    init_table();
    if (fread(table_, sizeof(uint64_t), 8*blocks_, fp) != (8*blocks_)) {return false;}

    return true;
  }

  size_t count() const {
    unsigned char *t = (unsigned char *) table_;
    size_t c = 0;
    for (size_t i = 0; i < 64*blocks_; i++) {
      unsigned char u = t[i];
      for (size_t j = 128; j != 0; j = j>>1) {
        if ((u & j) != 0) {
          c++;
        }
      }
    }
    std::cout << c << " bits set out of " << size_ << " with k = " << k_ << std::endl;
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
    blocks_ = 0;
  }

 private:



  void init_table() {
    fast_div_ = libdivide::divider<uint64_t>(blocks_);
    table_ = new uint64_t[8*blocks_];
    memset(table_, 0, 8*blocks_);
  }

  void init_k(size_t bits) {
    size_t k = (size_t) (bits*log(2));
    if (fpp(bits,k) < fpp(bits,k+1)) {
      k_ = k;
    } else {
      k_ = k+1;
    }
    std::cerr << "k="<<k_<<", fpp="<<fpp(bits,k_) <<  std::endl;
  }

  double fpp(size_t bits, size_t k) const {
    //    cout << bits<<","<<k<<","<<(-((double)k)/((double)bits)) << endl;
    return pow(1-exp(-((double)k)/((double)bits)),(double)k);
  }

  uint64_t rndup512(uint64_t x) const {
    return ((x+511)/512)*512;
  }
};

#endif // BFG_BLOCKEDBLOOMFILTER_HPP
