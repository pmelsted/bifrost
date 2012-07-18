%module graph
%{
#define SWIG_FILE_WITH_INIT
#include "Kmer.hpp"
#include "Contig.hpp"
#include "CompressedSequence.hpp"
#include "CompressedCoverage.hpp"
#include "BloomFilter.hpp"

%}
%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"
%include "Kmer.hpp"
%include "Contig.hpp"
%include "CompressedSequence.hpp"
%include "CompressedCoverage.hpp"

%include "BloomFilter.hpp"


using namespace std;




%ignore int2bin;
%extend Kmer {
  char *__str__() {
    static char tmp[Kmer::MAX_K+1];
    $self->toString(tmp);
    return &tmp[0];
  }
}

%extend Contig {
  const char *__str__() {
    static string s = $self->seq.toString();
    return s.c_str();
  }
}


%extend BloomFilter {
  bool __contains__(const Kmer km) {
    return $self->contains(km);
  }

  void open(const char *fn) {
    FILE *f = fopen(fn,"r");
    $self->ReadBloomFilter(f);
    fclose(f);
  }
} 
