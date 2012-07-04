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
  uint32_t seed;
  bool verbose;
  size_t bf;
  size_t bf2;
  vector<string> files;
  FilterReads_ProgramOptions() : k(0), nkmers(0), nkmers2(0), verbose(false), bf(4), bf2(8), seed(0) {}
};

// use:  FilterReads_PrintUsage();
// pre:   
// post: Information about how to "filter reads" has been printed to cerr
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
      "-B, --bloom-bits2=INT           Number of bits to use in second Bloom filter (default=8)" << endl <<
      "-s, --seed=INT                  Seed used for randomisazaion (default time based)" << endl <<
      "    --verbose                   Print lots of messages during run" << endl << endl
      ;
}


// use:  FilterReads_ParseOptions(argc, argv, opt);
// pre:  argc is the parameter count, argv is a list of valid parameters for 
//       "filtering reads" and opt is ready to contain the parsed parameters
// post: All the parameters from argv have been parsed into opt
void FilterReads_ParseOptions(int argc, char **argv, FilterReads_ProgramOptions &opt) {
  int verbose_flag = 0;
  const char* opt_string = "n:N:k:o:b:B:s:";
  static struct option long_options[] =
  {
    {"verbose", no_argument,  &verbose_flag, 1},
    {"kmer-size", required_argument, 0, 'k'},
    {"num-kmers", required_argument, 0, 'n'},
    {"num-kmers2", required_argument, 0, 'N'},
    {"output", required_argument, 0, 'o'},
    {"bloom-bits", required_argument, 0, 'b'},
    {"bloom-bits2", required_argument, 0, 'B'},
    {"seed", required_argument, 0, 's'},
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
    case 's':
      opt.seed = atoi(optarg);
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
// pre:  opt contains parameters for "filtering reads"
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
  /**
   *  outline of algorithm
   *   - create two bloom filters, BF and BF2
   *   - for each read in all files 
   *     - for all kmers in read
   *       - if kmer in BF 
   *         - insert kmer into BF2
   *       - else 
   *         - insert kmer into BF
   *  now BF2 contains at least all the kmers that appear once
   */


  uint32_t seed = opt.seed;
  if (seed == 0) {
    seed = (uint32_t) time(NULL);
  }
  BloomFilter BF(opt.nkmers, (size_t) opt.bf, seed);
  BloomFilter BF2(opt.nkmers2, (size_t) opt.bf2, seed+1); // use different seeds
  
  char name[8196],s[8196];
  size_t name_len,len;

  uint64_t n_read = 0;
  uint64_t num_kmers = 0;  
  uint64_t num_ins = 0;

  FastqFile FQ(opt.files);


  size_t chunk_size = 20000, reads_now = 0;
  vector<string> readv;
  bool done = false;
  #ifdef _OPENMP
  //omp_set_num_threads(1);
  #endif

  while (!done) {
    readv.clear();
    reads_now = 0;
    while (reads_now < chunk_size) {
      if (FQ.read_next(name, &name_len, s, &len, NULL, NULL) >= 0) {
        readv.push_back(string(s));
        ++n_read;
        ++reads_now;
      } else {
        done = true;
        break;
      }
    }

    KmerIterator iter, iterend;
    #pragma omp parallel for private(iter) shared(iterend, readv, BF, reads_now) reduction(+:num_ins)
    for (size_t index = 0; index < reads_now; ++index) {
      iter = KmerIterator(readv[index].c_str());
      for (;iter != iterend; ++iter) {
        //++num_kmers;
        Kmer km = iter->first;
        Kmer rep = km.rep();
        size_t r = BF.search(rep);        
        if (r == 0) {
          if (!BF2.contains(rep)) {
            BF2.insert(rep);
            ++num_ins;
          }
        } else {
          if (BF.insert(rep) == r) {
            ++num_ins;
          } else {
            if (opt.verbose) {
              cerr << "clash!" << endl;
            }
            BF2.insert(rep); // better safe than sorry
          }
        }
      }
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
//       arguments for "filtering the reads", including filenames
// post: If the number of arguments is correct and the arguments are valid
//       the "reads have been filtered" and written to a file 
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
