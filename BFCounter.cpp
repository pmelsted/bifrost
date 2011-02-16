#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <stdint.h>

#include <sys/stat.h>


#include <functional>


#include "google/sparse_hash_map"
#include "google/sparsehash/sparsehashtable.h"
#include "fastq.hpp"
#include "kmer.hpp"
#include "bloom_filter.hpp"
#include "hash.hpp"

using namespace std;
using google::sparse_hash_map;

struct ProgramOptions {
  size_t k;
  size_t nkmers;
  vector<string> files;
};

struct KmerIntPair {
  KmerIntPair() {};
  KmerIntPair(const Kmer &km, unsigned int k);

  char v[sizeof(Kmer)+sizeof(char)];
  unsigned int GetVal() const;
  void SetVal(const unsigned int k);
  const Kmer& GetKey() const;
  void SetKey(const Kmer& km);

  static const size_t KmerOffset;
  static const size_t IntOffset;
};


const size_t KmerIntPair::KmerOffset = 0;
const size_t KmerIntPair::IntOffset = sizeof(Kmer);

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

struct SelectKmerKey {
  const Kmer& operator()(const KmerIntPair &p) const {
    return p.GetKey();
  }
};

struct SetKmerKey {
  void operator()(KmerIntPair *value, const Kmer& km) {
    memcpy(value + KmerIntPair::KmerOffset, &km, sizeof(Kmer));
  }
};

void PrintUsage() {
  cerr << "Usage: BFCounter k nkmers fasta [fasta ...]" << endl;
}

void ParseOptions(int argc, char **argv, ProgramOptions &opt) {
  if (argc > 3) {
    stringstream ss(argv[1]);
    if ((ss >> opt.k).fail()) {
      cerr << "Error: k not an integer" << endl;
      PrintUsage();
      exit(1);
    }

    stringstream ss2(argv[2]);
    if ((ss2 >> opt.nkmers).fail()) {
      cerr << "Error: nkmers not an integer" << endl;
      PrintUsage();
      exit(1);
    }

    for (int i = 3; i < argc; i++) {
      opt.files.push_back(argv[i]);
    }
    
    struct stat stFileInfo;
    vector<string>::iterator it;
    int intStat;
    for(it = opt.files.begin(); it != opt.files.end(); ++it) {
      intStat = stat(it->c_str(), &stFileInfo);
      if (intStat != 0) {
	cerr << "Error: file not found, " << *it << endl;
	PrintUsage();
	exit(1);
      }
    }

  } else {
    PrintUsage();
    exit(1);
  }
}

typedef google::sparse_hashtable<KmerIntPair, Kmer, KmerHash, SelectKmerKey, SetKmerKey, std::equal_to<Kmer>, std::allocator<KmerIntPair> > hmap_t;

void debugKmer(const ProgramOptions &opt) {
  const char *s = "TCACAGTGTTGAACCTTTGTTTGGATGGAGCAGTTAGTGTTGAACCTTTGTTTGGATGGAGCAGTTAGTGTTGAACCTTTGTTTGGATGGAGCAGTTAGTGTTGAACCTTTGTTTGGATGGAGCAGTT";
  char tmp[1024];
  int len = string(s).size();

  Kmer::set_k(opt.k);
  size_t k = opt.k;

  hmap_t kmap;
  Kmer km(s);
  hmap_t::iterator it;
  cerr << s << endl;
  for (int i = 0; i <= len-k; i++) {
    if (i > 0) {
      km = km.forwardBase(s[i+k-1]);
    }
    km.toString(tmp); cerr << tmp;

    it = kmap.find(km);
    if (it!= kmap.end()) {
      cerr << " found " << it->GetVal();
      it->SetVal(it->GetVal()+1);
      cerr << " -> " << it->GetVal() << endl;
    } else {
      kmap.insert(KmerIntPair(km,1));
      cerr << " inserted " << endl; // << kmap.find(km)->GetVal() << endl;
    }
  }
  cerr << endl << endl;

  for (it = kmap.begin(); it != kmap.end(); ++it) {
    km = it->GetKey(); km.toString(tmp);
    cerr << tmp << " " << it->GetVal() << endl;;
  }

  
}

void CountBF(const ProgramOptions &opt) {
  Kmer::set_k(opt.k);
  size_t k = Kmer::k;

  //typedef sparse_hash_map<Kmer, uint32_t, KmerHash> hmap_t;


  hmap_t kmap;
  
  bloom_filter BF(opt.nkmers, (size_t) 4, (unsigned long) time(NULL));


  char name[8196],s[8196];
  size_t name_len,len;

  uint64_t n_read = 0;
  uint64_t num_kmers = 0;  
  uint64_t filtered_kmers = 0;
  uint64_t total_cov = 0;


  FastqFile FQ(opt.files);

  while (FQ.read_next(name, &name_len, s, &len, NULL) >= 0) {
    // add code to handle N's
    Kmer km(s);
    for (size_t i = 0; i <= len-k; ++i) {
      num_kmers++;
      if (i > 0) {
	km = km.forwardBase(s[i+k-1]);
      }
      Kmer tw = km.twin();
      Kmer rep = (km < tw) ? km : tw;
      if (BF.contains(rep)) {
	// has no effect if already in map
	pair<hmap_t::iterator,bool> ref = kmap.insert(KmerIntPair(rep,0));
      } else {
	BF.insert(rep);
      }
    }
    ++n_read;

    if (n_read % 1000000 == 0) {
      cerr << "processed " << n_read << " reads" << endl;
    }
  }

  cerr << "re-open all files" << endl;
  // close all files, reopen and get accurate counts;
  FQ.reopen();
  hmap_t::iterator it;

  while (FQ.read_next(name, &name_len, s, &len, NULL) >= 0) {
    //cerr << "read " << len << " characters" << endl;
    Kmer km(s);
    for (size_t i = 0; i <= len-k; ++i) {
      
      if (i > 0) {
	km = km.forwardBase(s[i+k-1]);
      }
      Kmer tw = km.twin();
      Kmer rep = (km < tw) ? km : tw;

      it = kmap.find(rep);
      if (it != kmap.end()) {
	//cerr << "Val before: " <<  it->GetVal();
	it->SetVal(it->GetVal()+1); // add 1 to count
	total_cov += 1;
	//cerr << ", after: " << it->GetVal() << endl;
      }
    }
  }
  
  FQ.close();

  cerr << "closed all files" << endl;

  // the hash map needs an invalid key to mark as deleted
  Kmer km_del;
  km_del.set_deleted();
  kmap.set_deleted_key(km_del);
  size_t n_del =0 ;

  for(it = kmap.begin(); it != kmap.end(); ) {
    if (it->GetVal() <= 1) {
      hmap_t::iterator del(it);
      ++it;
      // remove k-mer that got through the bloom filter
      kmap.erase(del);
      ++n_del;
    } else {
      ++it;
    }
  }

  total_cov -= n_del;

  cerr << "processed " << num_kmers << " kmers in " << n_read  << " reads"<< endl;
  cerr << "found " << kmap.size() << " non-filtered kmers, removed " << n_del << endl;
  filtered_kmers = num_kmers - total_cov;
    
  cerr << "total coverage " << total_cov << ", estimated number of kmers " << filtered_kmers << endl;
  cerr << "average coverage " << (total_cov / ((double) kmap.size())) << endl;
  cerr << num_kmers << endl  << filtered_kmers << endl << kmap.size() << endl;
}




int main(int argc, char **argv) {

  //cerr << "sizeof(kmer) " << sizeof(Kmer) << ", sizeof(pair) " << sizeof(KmerIntPair) << endl;

  //parse command line options
  ProgramOptions opt;
  ParseOptions(argc,argv,opt);

  //call main function
  CountBF(opt);
  //debugKmer(opt);
}
