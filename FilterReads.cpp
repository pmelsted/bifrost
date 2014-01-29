#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "Common.hpp"
#include "FilterReads.hpp"
#include "HashTables.hpp"
#include "fastq.hpp"
#include "Kmer.hpp"
#include "KmerIterator.hpp"
//#include "BloomFilter.hpp"
#include "BlockedBloomFilter.hpp"


struct FilterReads_ProgramOptions {
  bool verbose;
  size_t threads, read_chunksize, k, nkmers, nkmers2;
  FILE *outputfile;
  string output;
  size_t bf, bf2;
  uint32_t seed;
  vector<string> files;
  FilterReads_ProgramOptions() : verbose(false), threads(1), k(0), nkmers(0), nkmers2(0), \
                                 outputfile(NULL), bf(4), bf2(8), seed(0), read_chunksize(20000) {}
};

// use:  FilterReads_PrintUsage();
// pre:   
// post: Information about how to filter reads has been printed to cerr
void FilterReads_PrintUsage() {
  cerr << "BFGraph " << BFG_VERSION << endl;
  cerr << "Filters errors in fastq or fasta files and saves results to a file specified by -o or --output" << endl << endl;
  cerr << "Usage: BFGraph filter [options] ... FASTQ files";
  cerr << endl << endl << "Options:" << endl <<
      "  -v, --verbose               Print lots of messages during run" << endl <<
      "  -t, --threads=INT           Number of threads to use (default 1)" << endl << 
      "  -c, --chunk-size=INT        Read chunksize to split betweeen threads (default 20000 for multithreaded else 1)" << endl <<
      "  -k, --kmer-size=INT         Size of k-mers, the same value as used for filtering reads" << endl << 
      "  -n, --num-kmers=LONG        Estimated number of k-mers (upper bound)" << endl <<
      "  -N, --num-kmer2=LONG        Estimated number of k-mers in genome (upper bound)" << endl <<
      "  -o, --output=STRING         Filename for output" << endl <<
      "  -b, --bloom-bits=INT        Number of bits to use in Bloom filter (default=4)" << endl <<
      "  -B, --bloom-bits2=INT       Number of bits to use in second Bloom filter (default=8)" << endl <<
      "  -s, --seed=INT              Seed used for randomization (default time based)" 
  << endl << endl;
}


// use:  FilterReads_ParseOptions(argc, argv, opt);
// pre:  argc is the parameter count, argv is a list of valid parameters for 
//       "filtering reads" and opt is ready to contain the parsed parameters
// post: All the parameters from argv have been parsed into opt
void FilterReads_ParseOptions(int argc, char **argv, FilterReads_ProgramOptions &opt) {
  const char* opt_string = "vt:k:n:N:o:b:B:s:c:";
  static struct option long_options[] = {
      {"verbose",     no_argument,       0, 'v'},
      {"threads",     required_argument, 0, 't'},
      {"chunk-size",  required_argument, 0, 'c'},
      {"kmer-size",   required_argument, 0, 'k'},
      {"num-kmers",   required_argument, 0, 'n'},
      {"num-kmers2",  required_argument, 0, 'N'},
      {"output",      required_argument, 0, 'o'},
      {"bloom-bits",  required_argument, 0, 'b'},
      {"bloom-bits2", required_argument, 0, 'B'},
      {"seed",        required_argument, 0, 's'},
      {0,             0,                 0,  0 }
  };

  int option_index = 0, c;
  stringstream ss, ss2;
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
    switch (c) {
      case 0: 
        break;
      case 'v': 
        opt.verbose = true; 
        break;
      case 't':
        opt.threads = atoi(optarg);
        break;
      case 'k': 
        opt.k = atoi(optarg); 
        break;
      case 'n': 
        ss << optarg;
        ss >> opt.nkmers;
        break;
      case 'N':
        ss2 << optarg;
        ss2 >> opt.nkmers2;
        break;
      case 'o': 
        opt.output = optarg;
        break;
      case 'c':
        opt.read_chunksize = atoi(optarg);
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
  while (optind < argc) {
    opt.files.push_back(argv[optind++]);
  }
}


// use:  b = FilterReads_CheckOptions(opt);
// pre:  opt contains parameters for "filtering reads"
// post: (b == true)  <==>  the parameters are valid
bool FilterReads_CheckOptions(FilterReads_ProgramOptions &opt) {
  bool ret = true;
  
  size_t max_threads = 1;
  #ifdef _OPENMP
    max_threads = omp_get_max_threads();
  #endif
  if (opt.threads == 0 || opt.threads > max_threads) {
    cerr << "Error: Invalid number of threads " << opt.threads;
    if (max_threads == 1) {
      cerr << ", can only use 1 thread on this system" << endl;
    } else {
      cerr << ", need a number between 1 and " << max_threads << endl;
    }
    ret = false;
  }
  
  if (opt.read_chunksize == 0) {
    cerr << "Error: Invalid chunk-size: " << opt.read_chunksize
         << ", need a number greater than 0" << endl;
    ret = false;
  } else if (opt.threads == 1) {
    cerr << "Setting chunksize to 1 because of only 1 thread" << endl;
    opt.read_chunksize = 1;
  }

  if (opt.k <= 0 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error, invalid value for kmer-size: " << opt.k << endl;
    cerr << "Values must be between 1 and " << (MAX_KMER_SIZE-1) << endl;
    ret = false;
  }

  if (opt.nkmers <= 0) {
    cerr << "Error, invalid value for num-kmers (parameter -n): " << opt.nkmers << endl;
    cerr << "Values must be positive integers" << endl;
    ret = false;
  }

  if (opt.nkmers2 <= 0) {
    cerr << "Error, invalid value for num-kmers2 (parameter -N):" << opt.nkmers2 << endl;
    cerr << "Values must be positive integers" << endl;
    ret = false;
  }
  
  opt.outputfile = fopen(opt.output.c_str(), "wb");

  if (opt.outputfile == NULL) {
    cerr << "Error: Could not open file for writing, " << opt.outputfile << endl;
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

  return ret;
}


// use:  FilterReads_PrintSummary(opt);
// pre:  opt has information about Kmer size, Bloom Filter sizes
// post: Information about the two Bloom Filters has been printed to cerr 
void FilterReads_PrintSummary(const FilterReads_ProgramOptions &opt) {
  double fp;
  cerr << "Kmer size: " << opt.k << endl
       << "Chunksize: " << opt.read_chunksize << endl
       << "Using bloom filter size: " << opt.bf << " bits per element" << endl
       << "Estimated false positive rate: ";
  fp = pow(pow(.5,log(2.0)),(double) opt.bf);
  cerr << fp << endl;  
  
  cerr << "Using bloom filter size for second set: " << opt.bf2 << " bits per element" << endl;
  cerr << "Estimated false positive rate for second set: ";
  fp = pow(pow(.5,log(2.0)),(double) opt.bf2);
  cerr << fp << endl;
}


// use:  FilterReads_Normal(opt);
// pre:  opt has information about Kmer size, Bloom Filter sizes,
//       lower bound of Kmer count, upper bound of Kmer count
//       input file name strings and outputfile
// post: Input files have been opened and and all the reads 
//       have been broken into kmers of given Kmer size. 
//       The kmers have been filtered through two Bloom Filters
//       and those that survived through the second Bloom Filter have
//       been written into the outputfile
void FilterReads_Normal(const FilterReads_ProgramOptions &opt) {
  /**
   *  outline of algorithm
   *    create two bloom filters, BF and BF2
   *    for each read in all files 
   *      for all kmers in read
   *        if kmer in BF 
   *          insert kmer into BF2
   *        else 
   *          insert kmer into BF
   *  now BF2 contains at least all the kmers that appear once
   */

  #ifdef _OPENMP
    omp_set_num_threads(opt.threads);
  #endif

  uint32_t seed = opt.seed;
  if (seed == 0) {
    seed = (uint32_t) time(NULL);
  }
  
  BloomFilter BF(opt.nkmers, (size_t) opt.bf, seed);
  BloomFilter BF2(opt.nkmers2, (size_t) opt.bf2, seed + 1); // use different seeds
  
  bool done = false;
  char name[8192], s[8192];
  size_t name_len, len, read_chunksize = opt.read_chunksize;
  uint64_t n_read = 0, num_kmers = 0, num_ins = 0;

  FastqFile FQ(opt.files);
  vector<string> readv;

  while (!done) {
    readv.clear();
    size_t reads_now = 0;
    while (reads_now < read_chunksize) {
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
    #pragma omp parallel for private(iter) shared(iterend, readv, BF, reads_now) reduction(+:num_ins,num_kmers)
    for (size_t index = 0; index < reads_now; ++index) {
      iter = KmerIterator(readv[index].c_str());
      for (;iter != iterend; ++iter) {
        ++num_kmers;
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
  cerr << "Closed all fasta/fastq files" << endl;

  if (opt.verbose) {
    cerr << "processed " << num_kmers << " kmers in " << n_read  << " reads"<< endl;
    cerr << "found " << num_ins << " non-filtered kmers" << endl;
  }

  if (opt.verbose) {
    cerr << "Writing bloom filter to " << opt.output << endl
         << "Bloom filter size is " << num_ins << endl;
  }


  // First write metadata for bloom filter to opt.outputfile,
  // then the actual filter
  if (!BF2.WriteBloomFilter(opt.outputfile)) {
    cerr << "Error writing data to file: " << opt.output << endl;
  }

  fclose(opt.outputfile);

  if (opt.verbose) {
    cerr << " done" << endl;
  }

  cout << "Bloomfilter 1 count: " << BF.count() << endl;
  cout << "Bloomfilter 2 count: " << BF2.count() << endl;
    
}

// use:  FilterReads(argc, argv);
// pre:  argc is the number of arguments in argv and argv includes 
//       arguments for filtering reads, including filenames
// post: If the number of arguments is correct and the arguments are valid
//       the reads have been filtered and written to a file 
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
