#include <iostream>
#include <cstdio>
#include <string>
#include <getopt.h>

#include <sys/stat.h>

#include "Common.hpp"
#include "DumpBF.hpp"
#include "HashTables.hpp"
#include "Kmer.hpp"
#include "KmerIntPair.hpp"

struct DumpBF_ProgramOptions {
  size_t k;
  bool verbose;
  string input;
  string output;
  DumpBF_ProgramOptions() : k(0), verbose(false) {}
};

void DumpBF_PrintUsage() {
  // empty...
}

void DumpBF_ParseOptions(int argc, char **argv, DumpBF_ProgramOptions &opt) {
  int verbose_flag = 0;
  const char* opt_string = "k:i:o:";
  static struct option long_options[] =
    {
      {"verbose", no_argument,  &verbose_flag, 1},
      {"kmer-size", required_argument, 0, 'k'},
      {"input", required_argument, 0, 'i'},
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
    case 'i': 
      opt.input = optarg;
      break;
    default: break;
    }
  }
  
  if (verbose_flag) {
    opt.verbose = true;
  }
}

bool DumpBF_CheckOptions(const DumpBF_ProgramOptions &opt) {
  bool ret = true;
  
  if (opt.k <= 0 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error, invalid value for kmer-size: " << opt.k << endl;
    cerr << "Values must be between 1 and " << (MAX_KMER_SIZE-1) << endl;
    ret = false;
  }

  struct stat stFileInfo;
  int intStat;
  // check inputfile
  intStat = stat(opt.input.c_str(), &stFileInfo);
  if (intStat != 0) {
    cerr << "Error: input file not found, " << opt.input << endl;
    ret = false;
  }

  // TODO: check if we can write to output.

  return ret;
}


void DumpBF(int argc, char **argv) {
  DumpBF_ProgramOptions opt;
  DumpBF_ParseOptions(argc, argv, opt);

  if (!DumpBF_CheckOptions(opt)) {
    DumpBF_PrintUsage();
    exit(1);
  }

  Kmer::set_k(opt.k);
  size_t k = Kmer::k;

  // load hashtable from file
  hmap_t kmap;
  FILE* f = fopen(opt.input.c_str(), "rb");
  if (f == NULL) {
    cerr << "Could not open file " << opt.input << endl;
    exit(1);
  } else {
    kmap.read_metadata(f);
    kmap.read_nopointer_data(f);
    fclose(f);
    f = NULL;
  }

  // set deleted key
  Kmer km_del;
  km_del.set_deleted();
  kmap.set_deleted_key(km_del);


  // open output file
  FILE* of = fopen(opt.output.c_str(), "w");
  if (of == NULL) {
    cerr << "Could not open file for writing, " << opt.output << endl;
    exit(1);
  } else {
    // buffer
    char buf[1024];
    size_t cov;
    Kmer km;
    hmap_t::iterator it,it_end;
    it_end = kmap.end();
    for (it = kmap.begin(); it != kmap.end(); ++it) {
      km = it->GetKey();
      cov = it->GetVal();
      km.toString(&buf[0]);
      fprintf(of, "%s\t%zu\n",buf, cov);
    }
  }
}
