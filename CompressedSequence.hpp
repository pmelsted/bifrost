#ifndef BFG_COMPRESSED_SEQUENCE_HPP
#define BFG_COMPRESSED_SEQUENCE_HPP

#include <cstring>
#include <string>
#include <stdint.h>
#include "Kmer.hpp"

class CompressedSequence {
public: 
  CompressedSequence() : _length(0),_capacity(0),_data(0) {}
  ~CompressedSequence();
  CompressedSequence(const CompressedSequence& o);
  CompressedSequence& operator=(const CompressedSequence& o);
  explicit CompressedSequence(const char *s);
  explicit CompressedSequence(const string &s);
  explicit CompressedSequence(const Kmer &km);


  const char operator[](size_t index) const;


  size_t size() const { return _length; }
  Kmer getKmer(size_t offset) const;
  string toString() const;
  void toString(char *s) const;
  void toString(char *s, size_t offset, size_t length) const;
  string toString(size_t offset, size_t length) const;

  //  void setSequence(const CompressedSequence &o, size_t length, size_t offset = 0, bool reversed=false);
  void setSequence(const CompressedSequence &o, size_t start, size_t length, size_t offset = 0, bool reversed = false);
  void setSequence(const char *s, size_t length, size_t offset = 0, bool reversed=false);
  void setSequence(const string &s, size_t length, size_t offset = 0, bool reversed=false);
  void setSequence(const Kmer &km, size_t length, size_t offset = 0, bool reversed=false);
  
  void reserveLength(size_t new_length);

  CompressedSequence rev() const;
  size_t jump(const char *s, size_t i, size_t pos, bool reversed) const;
  //size_t endJump(char *s, size_t i, size_t dist, int pos, bool reversed) const;
  //size_t straightJump(char *s, size_t i, int pos) const;

private:
  size_t round_to_bytes(const size_t len) const { return (len+3)/4; }
  void _resize_and_copy(size_t new_cap, size_t copy_limit);
  uint32_t _length; // size of sequence
  uint32_t _capacity; // capacity of array allocated in bytes
  char *_data; // 0-based 2bit compressed dna string
};

#endif // BFG_COMPRESSED_SEQUENCE_H
