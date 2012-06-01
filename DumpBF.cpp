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
  bool quake;
  string input;
  string output;
  DumpBF_ProgramOptions() : k(0), verbose(false), quake(false) {}
};

void DumpBF_PrintUsage() {
  cerr << "BFGraph " << BFG_VERSION << endl << endl;
  cerr << "Writes k-mer occurrences into a tab-separated text file." << endl << endl;
  cerr <<     "Usage: BFGraph dump [options] ...";
  cerr << endl << endl <<
    "-k, --kmer-size=INT             Size of k-mers, at most " << (int) (Kmer::MAX_K-1)<< endl << 
    "-i, --input=STRING              Filename for input file from count" << endl << 
    "-o, --output=STRING             Filename for output" << endl <<
    "    --quake                     Write q-mers and quality scores to file for use with Quake (default=FALSE)" << endl <<
    "    --verbose                   Print lots of messages during run" << endl << endl
    ;

}

void DumpBF_ParseOptions(int argc, char **argv, DumpBF_ProgramOptions &opt) {
  int verbose_flag = 0;
  int quake_flag = 0;
  const char* opt_string = "k:i:o:";
  static struct option long_options[] =
    {
      {"verbose", no_argument,  &verbose_flag, 1},
      {"quake" , no_argument, &quake_flag, 1},
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
  if (quake_flag) {
    opt.quake = true;
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

void DumpBF_Quake(const DumpBF_ProgramOptions &opt) {
  hmapq_t kmap;

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
  
  Kmer km_del;
  km_del.set_deleted();
  kmap.set_deleted_key(km_del);



  FILE *of = fopen(opt.output.c_str(), "w");
  if (of == NULL) {
    cerr << "Could not open file for writing, " << opt.output << endl;
    exit(1);
  } else {
    char buf[1024];
    float qval;
    Kmer km;
    hmapq_t::iterator it, it_end;
    it_end = kmap.end();
    for (it = kmap.begin(); it != kmap.end(); ++it) {
      km = it->first;
      qval = it->second;
      km.toString(&buf[0]);
      fprintf(of, "%s\t%.6f\n",buf,qval);
    }

    fclose(of);
  }
  
}

void DumpBF_Normal(const DumpBF_ProgramOptions &opt) {
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
    
    fclose(of);
  }
}

void DumpBF(int argc, char **argv) {
  DumpBF_ProgramOptions opt;
  DumpBF_ParseOptions(argc, argv, opt);

  if (argc < 2) {
    DumpBF_PrintUsage();
    exit(1);
  }

  if (!DumpBF_CheckOptions(opt)) {
    cerr << endl << endl;
    DumpBF_PrintUsage();
    exit(1);
  }

  Kmer::set_k(opt.k);
  //size_t k = Kmer::k;

  if (opt.quake) {
    DumpBF_Quake(opt);
  } else {
    DumpBF_Normal(opt);
  }
}
