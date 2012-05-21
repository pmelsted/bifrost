#include <iterator>
#include <utility>
#include "Kmer.hpp"
#include "KmerIterator.hpp"

KmerIterator& KmerIterator::operator++() 
{
  int pos_ = p_.second;
  if (!invalid_) {
    if (s_[pos_+Kmer::k] == 0) {
      invalid_ = true;
      return *this;
    } else {
      find_next(pos_,pos_+Kmer::k-1,true);
      return *this;
    }
  };
  return *this;
}

KmerIterator KmerIterator::operator++(int) {
  KmerIterator tmp(*this); 
  operator++(); 
  return tmp;
}


bool KmerIterator::operator==(const KmerIterator& o) {
  if (invalid_  || o.invalid_) {
    return invalid_ && o.invalid_;
  } else {
    return (s_ == o.s_) && (p_.second == o.p_.second);
  }
}

std::pair< Kmer, int>& KmerIterator::operator*() {
  return p_;
}

std::pair< Kmer, int>* KmerIterator::operator->() {
  return &(operator*());
}


// start
void KmerIterator::find_next(int i, int j, bool last_valid) {
  ++i;
  ++j;
  
  while (s_[j] != 0) {
    char c = s_[j];
    if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
      if (last_valid) {
	p_.first = p_.first.forwardBase(c);
	break; // default case, 
      } else {
	if (i + Kmer::k - 1 == j) {
	  p_.first= Kmer(s_+i);
	  last_valid = true;
	  break; // create k-mer from scratch
	} else {
	  ++j;
	}
      }
    } else {
      ++j;
      i = j;
      last_valid = false;
    }
  }
  if (i+Kmer::k-1 == j && s_[j] != 0) {
      p_.second = i;
  } else {
    invalid_ = true;
  }
  
}
