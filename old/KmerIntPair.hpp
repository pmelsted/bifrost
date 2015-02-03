#ifndef BFG_KMER_INT_PAIR_HPP
#define BFG_KMER_INT_PAIR_HPP

#include "Kmer.hpp"

// TODO: switch out for a templated version
struct KmerIntPair {
  KmerIntPair() {};
  KmerIntPair(const Kmer &km, unsigned int k);

  char v[sizeof(Kmer)+sizeof(char)];
  unsigned int GetVal() const;
  void SetVal(const unsigned int k);
  const Kmer& GetKey() const;
  void SetKey(const Kmer& km);

  static const size_t KmerOffset = 0;
  static const size_t IntOffset = sizeof(Kmer);
};



struct SelectKmerKey {
  const Kmer& operator()(const KmerIntPair &p) const {
    return p.GetKey();
  }
};

struct SetKmerKey {
  void operator()(KmerIntPair *value, const Kmer& km);
};

#endif // BFG_KMER_INT_PAIR_HPP
