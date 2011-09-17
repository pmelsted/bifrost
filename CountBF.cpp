#include <iostream>
#include <cstdlib>
#include <ctime>

#include <sstream>
#include <vector>
#include <string>

#include <stdint.h>
#include <sys/stat.h>
#include <functional>

#include <getopt.h>

#include "Common.hpp"
#include "CountBF.hpp"

#include "HashTables.hpp"
#include "fastq.hpp"
#include "Kmer.hpp"
#include "bloom_filter.hpp"


// structs for getopt



struct CountBF_ProgramOptions {
  size_t k;
  size_t nkmers;
  string output;
  bool verbose;
  vector<string> files;

  CountBF_ProgramOptions() : k(0), nkmers(0), verbose(false) {}
};

void CountBF_PrintUsage() {
  cerr << "Usage ...";
}




void CountBF_ParseOptions(int argc, char **argv, CountBF_ProgramOptions &opt) {
  int verbose_flag = 0;
  const char* opt_string = "n:k:o:";
  static struct option long_options[] =
  {
    {"verbose", no_argument,  &verbose_flag, 1},
    {"kmer-size", required_argument, 0, 'k'},
    {"num-kmers", required_argument, 0, 'n'},
    {"output", required_argument, 0, 'o'},
    {0,0,0,0}
  };

  int option_index = 0;
  int c;

  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0: break;
    case 'k': 
      opt.k = atoi(optarg); 
      break;
    case 'o': 
      opt.output = optarg;
      break;
    case 'n': 
      opt.nkmers = atoi(optarg);
      break;
    default: break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
  
  if (verbose_flag) {
    opt.verbose = true;
  }
}


bool CountBF_CheckOptions(const CountBF_ProgramOptions &opt) {
  bool ret = true;

  if (opt.k <= 0 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error, invalid value for kmer-size: " << opt.k << endl;
    cerr << "Values must be between 1 and " << (MAX_KMER_SIZE-1) << endl;
    ret = false;
  }

  if (opt.nkmers <= 0) {
    cerr << "Error, invalid value for num-kmers: " << opt.nkmers << endl;
    cerr << "Values must be positive integers" << endl;
    ret = false;
  }

  if (opt.files.size() == 0) {
    cerr << "Need to specify files for input" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    vector<string>::const_iterator it;
    int intStat;
    for(it = opt.files.begin(); it != opt.files.end(); ++it) {
      intStat = stat(it->c_str(), &stFileInfo);
      if (intStat != 0) {
	cerr << "Error: file not found, " << *it << endl;
	ret = false;
      }
    }
  }
  
  //TODO: check if we have permission to write to outputfile
  
  return ret;

}


void CountBF(int argc, char **argv) {
  
  CountBF_ProgramOptions opt;
  CountBF_ParseOptions(argc,argv,opt);
  
  if (!CountBF_CheckOptions(opt)) {
    CountBF_PrintUsage();
    exit(1);
  }
  
  Kmer::set_k(opt.k);
  size_t k = Kmer::k;


  // create hash table and bloom filter
  hmap_t kmap;
  bloom_filter BF(opt.nkmers, (size_t) 4, (unsigned long) time(NULL));


  
  char name[8196],s[8196];
  size_t name_len,len;

  uint64_t n_read = 0;
  uint64_t num_kmers = 0;  
  uint64_t filtered_kmers = 0;
  uint64_t total_cov = 0;

  // loops over all files
  FastqFile FQ(opt.files);

  // for each read
  while (FQ.read_next(name, &name_len, s, &len, NULL) >= 0) {
    // TODO: add code to handle N's, currently all N's are mapped to A
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

    if (opt.verbose && n_read % 1000000 == 0) {
      cerr << "processed " << n_read << " reads" << endl;
    }
  }

  cerr << "re-open all files" << endl;
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

  if (opt.verbose) {
    cerr << "processed " << num_kmers << " kmers in " << n_read  << " reads"<< endl;
    cerr << "found " << kmap.size() << " non-filtered kmers, removed " << n_del << endl;
    filtered_kmers = num_kmers - total_cov;
    
    cerr << "total coverage " << total_cov << ", estimated number of kmers " << filtered_kmers << endl;
    cerr << "average coverage " << (total_cov / ((double) kmap.size())) << endl;
    cerr << num_kmers << endl  << filtered_kmers << endl << kmap.size() << endl;
  }

  if (opt.verbose) {
    cerr << "Writing hash table to file " << opt.output << " .. "; cerr.flush();
  }
  FILE* f = fopen(opt.output.c_str(), "wb");
  if (f == NULL) {
    cerr << "Error could not write to file!" << endl;
  } else {
    // first metadata for hash table
    kmap.write_metadata(f);
    // then the actual hashtable
    kmap.write_nopointer_data(f);
    fclose(f);
    f = NULL;
  }
  if (opt.verbose) {
    cerr << " .. done" << endl;
  }
}
