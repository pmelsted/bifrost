#ifndef BFG_KMER_HPP
#define BFG_KMER_HPP

#ifndef MAX_KMER_SIZE
 #define MAX_KMER_SIZE 32
#endif

#include <stdio.h>
#include <stdint.h>
#include <cassert>
#include <cstring>

#include "hash.hpp"


//for debug

char *int2bin(uint32_t a, char *buffer, int buf_size);




class Kmer {
 public:

  Kmer();
  Kmer(const Kmer& o);
  explicit Kmer(const char *s); 
  

  Kmer& operator=(const Kmer& o);
  
  void set_deleted();

  bool operator<(const Kmer& o) const;

  bool operator==(const Kmer& o) const;

  bool operator!=(const Kmer& o) const {
    return !(*this == o);
  }

  void set_kmer(const char *s);

  uint64_t hash() const;

  

  Kmer twin() const;

  Kmer getLink(const size_t index) const;

  Kmer forwardBase(const char b) const;

  Kmer backwardBase(const char b) const;
  
  void printBinary() const;
  
  void toString(char * s) const;

  // static functions
  static void set_k(unsigned int _k);


  static const unsigned int MAX_K = MAX_KMER_SIZE;
  static unsigned int k;

 private:
  static unsigned int k_bytes;
  //  static unsigned int k_longs;
  static unsigned int k_modmask;

  // data fields
  union {
    uint8_t bytes[MAX_K/4];
    //uint32_t longs[MAX_K/16];
  };


  // private functions
  void shiftRight(int shift);

  void shiftLeft(int shift);
};


struct KmerHash {
  size_t operator()(const Kmer &km) const {
    return km.hash();
  }
};


#endif // BFC_KMER_HPP
