#ifndef BFG_BLOCKEDBLOOMFILTER_HPP
#define BFG_BLOCKEDBLOOMFILTER_HPP

#include "BloomFilter.hpp"


/* Short description: 
 *  - Extended BloomFilter which hashes into 64-bit blocks
 *    that can be accessed very fast from the CPU cache 
 * */
class BlockedBloomFilter : public BloomFilter {
private:
  uint64_t blocks;
public:
  
  BlockedBloomFilter() : blocks(0) {
    BloomFilter();
  }
  
  BlockedBloomFilter(size_t num, size_t bits, uint32_t seed) : blocks(0) {
    seed_ = seed; size_ = 0; table_ = NULL;
    size_ = rndup(bits*num);
    init_table();
    init_k(bits);
  }

  template<typename T>
  bool contains(T x) const {
    return (search(x) == 0);
  }


  template<typename T>
  size_t search(T x) const {
    size_t r = k_;
    uint64_t id;
    uint64_t hash;
    uint64_t block; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+2, &block);
    // block is the index of the 512 bit memory block where x would be stored
    // 0 <= block < blocks
    block = block - (block / fast_div_) * (blocks); // block % blocks 
    // Multiply block by 64 to get the first byte of the block 
    block <<= 6; 
    uint64_t hash0; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_, &hash0);
    hash0 |= 1; // make hash0 an odd number
    uint64_t hash1; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+1, &hash1);
    for (uint64_t i = 0; i < k_; i++) {
      //MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+i, &hash);
      hash = hash0 * i + hash1;
      // 0 <= id < 512, id represents one bit 
      id = hash & 0x1ff; // equal to hash % 512;
      // we check if bit number 1+(id % 8) in byte (table_[block + id/8]) is set
      if ((table_[block + (id >> 3)] & mask[id & 0x07]) == 0) {
	r--;
      }
    }
    return r;
  }

  template<typename T>
  size_t insert(T x) {
    size_t r = 0;
    unsigned char val;
    uint64_t id;
    uint64_t hash;
    uint64_t block; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+2, &block);
    // block is the index of the 512 bit memory block where x would be stored
    // 0 <= block < blocks
    block = block - (block / fast_div_) * (blocks); // block % blocks 
    // Multiply block by 64 to get the first byte of the block 
    block <<= 6;
    uint64_t hash0; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_  , &hash0);
    hash0 |= 1; // make hash0 an odd number
    uint64_t hash1; MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+1, &hash1);
    for(uint64_t i = 0; i < k_; i++) {
      //MurmurHash3_x64_64((const void*) &x, sizeof(T), seed_+i,&hash);
      hash = hash0 * i + hash1;
      // 0 <= id < 512, id represents one bit 
      id = hash & 0x1ff; // equal to hash % 512;
      // we set bit number 1+(id % 8) in byte (table_[block + id/8]) to 1
      if ((table_[block + (id>>3)] & mask[id & 0x07]) == 0) {
	val = __sync_fetch_and_or(table_ + block + (id>>3), mask[id & 0x07]);
	if ((val & mask[id & 0x07]) == 0) 
	  r++;
      }
	//table_[block + (id >> 3)] |= mask[id & 0x07];
    }
    return r;
  }


private:

  void clear() {
    BloomFilter::clear();
    blocks = 0;
  }

  void init_table() {
    blocks = 1 + (size_ >> 9);
    fast_div_ = libdivide::divider<uint64_t>(blocks);
    table_ = new unsigned char[size_>>3];
    memset(table_, 0, size_ >> 3);
  }
};

#endif // BFG_BLOCKEDBLOOMFILTER_HPP
