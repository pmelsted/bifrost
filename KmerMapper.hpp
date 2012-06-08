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



/* Short description: 
 *  - A ContigRef can be:
 *    1) An empty reference
 *    2) A reference to a contig
 *    3) A reference to position inside another ContigRef
 *  - A chain of references always ends at a contig
 *  */
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


/* Short description: 
 *  - Map kmers in a contig to ContigRefs
 *  - Store which ContigRefs refer to Contigs
 *  - Get contig at the end of a ContigRef chain
 *  - Get a ContigRef that a kmer maps to
 *  - Join contigs while maintaining right references
 *  - Pretty print a reference map of a Contig
 *  */
class KmerMapper {
public:
typedef google::sparse_hash_map<Kmer, ContigRef, KmerHash> hmap_contig_t;
  typedef hmap_contig_t::iterator iterator;
  typedef hmap_contig_t::const_iterator const_iterator;

  KmerMapper(size_t init = 1000) :  map(init) {stride = Kmer::k; contigs.reserve(init);}
  ~KmerMapper();

  ContigRef addContig(const string &s);
  ContigRef addContig(const char *s);
  
  ContigRef joinContigs(ContigRef a, ContigRef b);
  //  ContigRef extendContig(ContigRef a, const string &s);
  //  ContigRef extendContig(ContigRef a, const char *s);

  ContigRef find(const Kmer km); // maybe change 
  //const_iterator end() { return map.end(); }

  size_t stride; // store every stride-th kmer  
  const size_t size() const {
    return map.size();
  }

  ContigRef getContig(const size_t id) const;
  ContigRef getContig(const ContigRef ref) const;
  void printContig(const size_t id);
private:
  ContigRef find_rep(ContigRef a) const;

  hmap_contig_t map;
  vector<ContigRef> contigs;
};

#endif // BFG_KMERMAPPER_H
