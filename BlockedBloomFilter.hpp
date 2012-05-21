#include <cmath>
#include <cstring>
#include "hash.hpp"
#include "libdivide.h"
#include "stdint.h"
#include <iostream>
#include <cstdio>

using namespace libdivide;
using namespace std;

static const unsigned char mask[8] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};
 
class BlockedBloomFilter {
private:
  unsigned char* table_;
  uint32_t seed_;
  uint64_t size_;
  uint64_t blocks;
  size_t k_;
  divider<uint64_t> fast_div_;
  divider<uint64_t> block_div_;
public:
  BlockedBloomFilter() : seed_(0), size_(0), table_(NULL), k_(0), fast_div_() {}

  BlockedBloomFilter(size_t num, size_t bits, uint32_t seed) : seed_(seed), size_(0), table_(NULL), fast_div_(), block_div_(){
    cout << "num="<<num << ", bits="<<bits;
    size_ = rndup(bits*num);
	blocks = 1 + (size_ >> 9);
    cout <<", size=" << size_ << " , blocks=" << blocks << endl;

    init_table();
    init_k(bits);
  } 

  ~BlockedBloomFilter() {
    clear();
  }

  template<typename T>
  bool contains(T x)  {
    uint64_t id;
    uint64_t hash;
    uint64_t block; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+2, &block);
	block = block - (block / block_div_) * (blocks); // block % blocks 
	block *= (1 << 6);
    uint64_t hash0; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_, &hash0);
    hash0 = (hash0<<1)>>1 + 1; // odd number
    uint64_t hash1; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+1, &hash1);
    for (uint64_t i = 0; i < k_; i++) {
      //MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+i, &hash);
      hash = hash0 * i + hash1;
      id = hash - (hash / fast_div_) * (1 << 9); // equal to hash % 512;
      assert(id == (hash % 512));
      if ((table_[block + (id >> 3)] & mask[id & 0x07]) == 0) {
	return false;      
      }
    }
    return true;
  }

  template<typename T>
  void insert(T x) {
    uint64_t id;
    uint64_t hash;
    uint64_t block; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+2, &block);
	block = block - (block / block_div_) * (blocks); // block % blocks 
	block *= (1 << 6);
    uint64_t hash0; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_  , &hash0);
    hash0 = (hash0<<1)>>1 + 1; // odd number, why not hash0 |= 1?
    uint64_t hash1; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+1, &hash1);
    for(uint64_t i = 0; i < k_; i++) {
      //MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+i,&hash);
      hash = hash0 * i + hash1;
      id = hash - (hash / fast_div_) * (1 << 9); // equal to hash % 512;
      assert(id == (hash % 512));
      table_[block + (id >> 3)] |= mask[id & 0x07];
    }
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


private:

  void clear() {
    if (table_ != NULL) {
		delete[] table_;
    }
    table_ = NULL;
	blocks = 0;
    size_ = 0;
  }

  void init_table() {
    fast_div_ = libdivide::divider<uint64_t>(1 << 9);
    block_div_ = libdivide::divider<uint64_t>(blocks);
	table_ = new unsigned char[size_>>3];
	memset(table_, 0, size_ >> 3);
  }

  void init_k(size_t bits) {
    size_t k = (size_t) (bits*log(2));
    if (fpp(bits,k) < fpp(bits,k+1)) {
      k_ = k;
    } else {
      k_ = k+1;
    }
    cout << "k="<<k_<<", fpp="<<fpp(bits,k) << endl; // why not fpp(bits,k_) ?
  }

  double fpp(size_t bits, size_t k) {
    //    cout << bits<<","<<k<<","<<(-((double)k)/((double)bits)) << endl;
    return pow(1-exp(-((double)k)/((double)bits)),(double)k);
  }
  
  uint64_t rndup(uint64_t x) {
    return ((x+63) >> 6)<<6;
  }
   
};
