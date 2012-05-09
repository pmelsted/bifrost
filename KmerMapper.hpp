#ifndef BFG_KMERMAPPER_HPP
#define BFG_KMERMAPPER_HPP

#include <vector>
#include <string>

#include "Common.hpp"
#include "Kmer.hpp"
#include "CompressedSequence.hpp"
#include "Contig.hpp"

#include "google/sparse_hash_map"
using google::sparse_hash_map;



class ContigRef {
public:
  ContigRef() : isContig(true) { ref.contig = NULL;}
  ContigRef(uint32_t id, int32_t pos) : isContig(false) {ref.idpos.id=id; ref.idpos.pos = pos;}
  bool isEmpty() { return (isContig && ref.contig == NULL);}
  
  union ContigRefUnion_t{
    struct ContigRefProper_t {
      uint32_t id; // maps to ContigArray, managed by Mapper class
      int32_t pos; // negative is reverse
    } idpos;
    Contig *contig; // not managed by ContigRef, but by KmerMapper
  } ref;
  bool isContig;
};

class KmerMapper {
public:
typedef google::sparse_hash_map<Kmer, ContigRef, KmerHash> hmap_contig_t;
  typedef hmap_contig_t::iterator iterator;
  typedef hmap_contig_t::const_iterator const_iterator;

  KmerMapper() {stride = Kmer::k;}
  KmerMapper(size_t init) : contigs(init), map(init) {stride = Kmer::k;}
  ~KmerMapper();

  ContigRef addContig(const string &s);
  ContigRef addContig(const char *s);
  
  ContigRef joinContigs(ContigRef a, ContigRef b);
  //  ContigRef extendContig(ContigRef a, const string &s);
  //  ContigRef extendContig(ContigRef a, const char *s);

  ContigRef find(const Kmer km); // maybe change 
  //const_iterator end() { return map.end(); }

  size_t stride; // store every stride-th kmer  
private:
  ContigRef find_rep(ContigRef a);

  hmap_contig_t map;
  vector<ContigRef> contigs;
};




#endif // BFG_KMERMAPPER_H
