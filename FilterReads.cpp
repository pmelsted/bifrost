#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>

#include <sstream>
#include <vector>
#include <string>

#include <stdint.h>
#include <sys/stat.h>
#include <functional>

#include <getopt.h>

#include "Common.hpp"
#include "FilterReads.hpp"

#include "HashTables.hpp"
#include "fastq.hpp"
#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "BloomFilter.hpp"


// structs for getopt

struct FilterReads_ProgramOptions {
  size_t k;
  size_t nkmers;
  size_t nkmers2;
  string output;
  bool verbose;
  size_t bf;
  size_t bf2;
  vector<string> files;
  FilterReads_ProgramOptions() : k(0), nkmers(0), nkmers2(0), verbose(false), bf(4), bf2(8) {}
};

// use:  FilterReads_PrintUsage();
// pre:   
// post: Information about how the `./BFGraph filter` program is used 
//       has been printed to cerr
void FilterReads_PrintUsage() {
  cerr << "BFGraph " << BFG_VERSION << endl << endl;
  cerr << "Filters errors in fastq or fasta files and saves results" << endl << endl;
  cerr << "Usage: BFGraph filter [options] ... FASTQ files";
  cerr << endl << endl <<
    "-k, --kmer-size=INT             Size of k-mers, at most " << (int) (Kmer::MAX_K-1)<< endl << 
    "-n, --num-kmers=LONG            Estimated number of k-mers (upper bound)" << endl <<
    "-N, --num-kmer2=LONG            Estimated number of k-mers in genome (upper bound)" << endl <<
    "-o, --output=STRING             Filename for output" << endl <<
    "-b, --bloom-bits=INT            Number of bits to use in Bloom filter (default=4)" << endl <<
    "-B, --bloom-bits2=INT          Number of bits to use in second Bloom filter (default=8)" << endl <<
    "    --verbose                   Print lots of messages during run" << endl << endl
    ;
}


// use:  FilterReads_ParseOptions(argc, argv, opt);
// pre:  argc is the parameter count, argv is a list of valid parameters and
//       opt can contain the parsed parameters
// post: All the parameters from argv have been parsed into opt
void FilterReads_ParseOptions(int argc, char **argv, FilterReads_ProgramOptions &opt) {
  int verbose_flag = 0;
  const char* opt_string = "n:N:k:o:b:B:";
  static struct option long_options[] =
  {
    {"verbose", no_argument,  &verbose_flag, 1},
    {"kmer-size", required_argument, 0, 'k'},
    {"num-kmers", required_argument, 0, 'n'},
    {"num-kmers2", required_argument, 0, 'N'},
    {"output", required_argument, 0, 'o'},
    {"bloom-bits", required_argument, 0, 'b'},
    {"bloom-bits2", required_argument, 0, 'B'},
    {0,0,0,0}
  };

  int option_index = 0;
  int c;
  stringstream ss;
  stringstream ss2;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);
    //cout << "debug: c="<<c<<", optarg="<<optarg << endl;
    if (c == -1) {
      break;
    }

    switch (c) {
    case 0: 
      break;
    case 'k': 
      opt.k = atoi(optarg); 
      break;
    case 'o': 
      opt.output = optarg;
      break;
    case 'n': 
      ss << optarg;
      ss >> opt.nkmers;
      break;
    case 'N':
      ss2 << optarg;
      ss2>> opt.nkmers2;
      break;
    case 'b':
      opt.bf = atoi(optarg);
      break;
    case 'B':
      opt.bf2 = atoi(optarg);
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


// use:  b = FilterReads_CheckOptions(opt);
// pre:  opt contains parameters for the program BFGraph
// post: (b == true)  <==>  the parameters are valid
bool FilterReads_CheckOptions(FilterReads_ProgramOptions &opt) {
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

  if (opt.nkmers2 <= 0) {
    cerr << "Error, invalid value for num-kmers2: " << opt.nkmers2 << endl;
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
  
  if (opt.bf <= 0) {
    cerr << "Invalid value for bloom filter size" << endl;
    ret = false;
  }

  if (opt.bf2 <= 0) {
    cerr << "Invalid value for bloom filter size for second set" << endl;
    ret = false;
  }

  //TODO: check if we have permission to write to outputfile
  
  return ret;

}


// use:  FilterReads_PrintSummary(opt);
// pre:  opt has information about Kmer size, Bloom Filter sizes
// post: Information about the two Bloom Filters has been printed to cerr 
void FilterReads_PrintSummary(const FilterReads_ProgramOptions &opt) {
  double fp;
  cerr << "Using bloom filter size: " << opt.bf << " bits" << endl;
  cerr << "Estimated false positive rate: ";
  fp = pow(pow(.5,log(2.0)),(double) opt.bf);
  cerr << fp << endl;  
  
  cerr << "Using bloom filter size for second set: " << opt.bf2 << " bits" << endl;
  cerr << "Estimated false positive rate for second set: ";
  fp = pow(pow(.5,log(2.0)),(double) opt.bf2);
  cerr << fp << endl;
}


// use:  FilterReads_Normal(opt);
// pre:  opt has information about Kmer size, Bloom Filter sizes,
//       lower bound of Kmer count, upper bound of Kmer count
//       input file name strings and output file name
// post: Input files have been opened and and all the reads 
//       have been broken into kmers of given Kmer size. 
//       The kmers have been filtered through two Bloom Filters
//       and those that survived through the second Bloom Filter have
//       been written into the output file
void FilterReads_Normal(const FilterReads_ProgramOptions &opt) {
  // create hash table and bloom filter
  size_t k = Kmer::k;
  BloomFilter BF(opt.nkmers, (size_t) opt.bf, (uint32_t) time(NULL));
  BloomFilter BF2(opt.nkmers2, (size_t) opt.bf2, (uint32_t) time(NULL));
  
  char name[8196],s[8196];
  size_t name_len,len;

  uint64_t n_read = 0;
  uint64_t num_kmers = 0;  
  uint64_t num_ins = 0;

  // loops over all files
  FastqFile FQ(opt.files);

  // for each read
  while (FQ.read_next(name, &name_len, s, &len, NULL, NULL) >= 0) {
    // TODO: add code to handle N's, currently all N's are mapped to A

    KmerIterator it(s), it_end;
    Kmer km;
    for (;it != it_end; ++it) {
      ++num_kmers;
      km = it->first;
      Kmer tw = km.twin();
      Kmer rep = km.rep();
      if (BF.contains(rep)) {
        if (!BF2.contains(rep)) {
          BF2.insert(rep);
          ++num_ins;
        }
      } else {
          BF.insert(rep);
      }
    }
    ++n_read;

    if (opt.verbose && n_read % 1000000 == 0) {
      cerr << "processed " << n_read << " reads" << endl;
    }
  }
  
 
  FQ.close();

  cerr << "closed all files" << endl;

  if (opt.verbose) {
    cerr << "processed " << num_kmers << " kmers in " << n_read  << " reads"<< endl;
    cerr << "found " << num_ins << " non-filtered kmers" << endl;
  }

  if (opt.verbose) {
    cerr << "Writing bloom filter to file " << opt.output << " .. "; cerr.flush();
    cerr << "Bloom filter size is" << num_ins << endl;
  }


  FILE* f = fopen(opt.output.c_str(), "wb");
  if (f == NULL) {
    cerr << "Error could not write to file!" << endl;
  } else {
  // first metadata for bloom filter
  // then the actual filter
  if (!BF2.WriteBloomFilter(f)) {
    cerr << "Error writing data to file!" << endl;
  }
  fclose(f);
  f = NULL;
  }
  if (opt.verbose) {
    cerr << " done" << endl;
  }
}

// use:  FilterReads(argc, argv);
// pre:  argc is the number of arguments in argv and argv includes 
//       arguments for the Filtering
// post: If the number of arguments is correct and the arguments are valid
//       the Filtering has been completed
void FilterReads(int argc, char **argv) {
  
  FilterReads_ProgramOptions opt;
  FilterReads_ParseOptions(argc,argv,opt);

  if (argc < 2) {
    FilterReads_PrintUsage();
    exit(1);
  }
  
  if (!FilterReads_CheckOptions(opt)) {
    FilterReads_PrintUsage();
    exit(1);
  }
  
  // set static global k-value
  Kmer::set_k(opt.k);

  if (opt.verbose) {
    FilterReads_PrintSummary(opt);
  }

  FilterReads_Normal(opt);

}
