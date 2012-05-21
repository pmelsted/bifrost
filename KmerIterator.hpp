#ifndef BFG_KMER_ITERATOR_HPP
#define BFG_KMER_ITERATOR_HPP

#include <iterator>
#include "Kmer.hpp"


class KmerIterator : public std::iterator< std::input_iterator_tag, std::pair<Kmer, int>, int> {
public:
  KmerIterator() : s_(NULL),p_(), invalid_(true) {}
  KmerIterator(const char* s) : s_(s),p_(), invalid_(false) { find_next(-1,-1,false);}
  KmerIterator(const KmerIterator& o) : s_(o.s_), p_(o.p_), invalid_(o.invalid_) {}

  KmerIterator& operator++();
  KmerIterator operator++(int);

  bool operator==(const KmerIterator& o);
  bool operator!=(const KmerIterator& o) { return !this->operator==(o);}

  std::pair< Kmer, int>& operator*();
  std::pair< Kmer, int>* operator->();

private:
  void find_next(int i, int j, bool last_valid);
  
  const char *s_;
  std::pair<Kmer, int> p_;
  bool invalid_;
};

#endif
