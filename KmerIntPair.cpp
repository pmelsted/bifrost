#include <cstdlib>

#include "Kmer.hpp"
#include "KmerIntPair.hpp"


KmerIntPair::KmerIntPair(const Kmer &km, unsigned int val) {
  SetKey(km);
  SetVal(val);
}

void KmerIntPair::SetVal(unsigned int val) {
  char val8 = (val > 0xFF) ?  0xFF : (char)val;
  //memcpy(&this->v + KmerIntPair::IntOffset, &val8, sizeof(uint8_t));
  this->v[KmerIntPair::IntOffset] = val8;
}

unsigned int KmerIntPair::GetVal() const {
  //uint8_t tmp = *reinterpret_cast<const uint8_t*>(this+KmerIntPair::IntOffset);
  return (uint8_t)this->v[KmerIntPair::IntOffset];
}

const Kmer& KmerIntPair::GetKey() const {
  return *reinterpret_cast<const Kmer*>(this + KmerIntPair::KmerOffset);
}

void KmerIntPair::SetKey(const Kmer& km) {
  memcpy(this, &km, sizeof(Kmer));
}

void SetKmerKey::operator()(KmerIntPair *value, const Kmer& km) {
  memcpy(value + KmerIntPair::KmerOffset, &km, sizeof(Kmer));
}


