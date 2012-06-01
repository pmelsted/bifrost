#ifndef BFG_BLOCKEDBLOOMFILTER_HPP
#define BFG_BLOCKEDBLOOMFILTER_HPP

#include "BloomFilter.hpp"

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
  bool contains(T x)  {
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
      table_[block + (id >> 3)] |= mask[id & 0x07];
    }
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


#endif
