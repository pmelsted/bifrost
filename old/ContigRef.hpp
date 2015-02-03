#ifndef BFG_CONTIGREF_HPP
#define BFG_CONTIGREF_HPP

#include "Common.hpp"
#include "Kmer.hpp"
#include "KmerMapper.hpp"
#include "CompressedSequence.hpp"
#include "CompressedCoverage.hpp"

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
  
  union ContigRefUnion_t {
    struct ContigRefProper_t {
      uint32_t id; // maps to ContigArray, managed by Mapper class
      int32_t pos; // negative is reverse
    } idpos;
    Contig *contig; // not managed by ContigRef, but by KmerMapper
  } ref;
  bool isContig;
};

#endif // BFG_CONTIGREF_HPP
