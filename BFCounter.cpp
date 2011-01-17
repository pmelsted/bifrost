#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <stdint.h>

#include <sys/stat.h>


#include "google/sparse_hash_map"

#include "fastq.hpp"
#include "kmer.hpp"
#include "bloom_filter.hpp"
#include "hash.hpp"

using namespace std;
using google::sparse_hash_map;

struct ProgramOptions {
  int k;
  int nkmers;
  vector<string> files;
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


void CountBF(const ProgramOptions &opt) {
  Kmer::set_k(opt.k);
  size_t k = Kmer::k;

  typedef sparse_hash_map<Kmer, uint32_t, KmerHash> hmap_t;

  hmap_t kmap;

  bloom_filter BF(opt.nkmers, (size_t) 4, (unsigned long) time(NULL));

  char name[8196],s[8196];
  size_t name_len,len;

  size_t n_read = 0;
  
  FastqFile FQ(opt.files);
  size_t num_kmers = 0;
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
	pair<hmap_t::iterator,bool> ref = kmap.insert(make_pair(rep,0));
      } else {
	BF.insert(rep);
      }
    }
    ++n_read;

    if (n_read % 1000000 == 0) {
      cerr << "processed " << n_read << " reads" << endl;
    }
  }

  // close all files, reopen and get accurate counts;
  FQ.reopen();
  hmap_t::iterator it;

  while (FQ.read_next(name, &name_len, s, &len, NULL) >= 0) {
    Kmer km(s);
    for (size_t i = 0; i <= len-k; ++i) {
      if (i > 0) {
	km = km.forwardBase(s[i+k-1]);
      }
      Kmer tw = km.twin();
      Kmer rep = (km < tw) ? km : tw;

      it = kmap.find(rep);
      if (it != kmap.end()) {
	it->second += 1; // add 1 count
      }
    }
  }
  
  FQ.close();

  // the hash map needs an invalid key to mark as deleted
  Kmer km_del;
  km_del.set_deleted();
  kmap.set_deleted_key(km_del);
  size_t n_del =0 ;

  for(it = kmap.begin(); it != kmap.end(); ) {
    if (it->second <= 1) {
      hmap_t::iterator del(it);
      ++it;
      // remove k-mer that got through the bloom filter
      kmap.erase(del);
      ++n_del;
    } else {
      ++it;
    }
  }

  cerr << "processed " << num_kmers << " kmers in " << n_read << endl;
  cerr << "found " << kmap.size() << " unique kmers, removed " << n_del << endl;
  
}


int main(int argc, char **argv) {
  //parse command line options
  ProgramOptions opt;
  ParseOptions(argc,argv,opt);

  //call main function
  CountBF(opt);
}
