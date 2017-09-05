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

#include <thread>
#include <atomic>

#include "Common.hpp"
#include "FilterReads.hpp"
#include "fastq.hpp"
#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "BlockedBloomFilter.hpp"

#include "minHashIterator.hpp"
#include "RepHash.hpp"


struct FilterReads_ProgramOptions {
  bool verbose;
  size_t threads, read_chunksize, k, g, nkmers, nkmers2;
  FILE *outputfile;
  string output;
  size_t bf, bf2;
  uint32_t seed;
  bool ref;
  vector<string> files;
  FilterReads_ProgramOptions() : verbose(false), threads(1), k(0), g(21), nkmers(0), nkmers2(0), \
    outputfile(NULL), bf(4), bf2(8), seed(0), read_chunksize(10000), ref(false) {}
};

// use:  FilterReads_PrintUsage();
// pre:
// post: Information about how to filter reads has been printed to cerr
void FilterReads_PrintUsage() {
  cout << "BFGraph " << BFG_VERSION << endl;
  cout << "Filters errors in fastq or fasta files and saves results to a file specified by -o or --output" << endl << endl;
  cout << "Usage: BFGraph filter [options] ... FASTQ files";
  cout << endl << endl << "Options:" << endl <<
       "  -v, --verbose               Print lots of messages during run" << endl <<
       "  -t, --threads=INT           Number of threads to use (default 1)" << endl <<
       "  -c, --chunk-size=INT        Read chunksize to split betweeen threads (default 10000 for multithreaded else 1)" << endl <<
       "  -k, --kmer-size=INT         Size of k-mers" << endl <<
       "      --ref                   Reference mode, no filtering use only num_kmers and bloom-bits" << endl <<
       "  -g, --min-size=INT          Size of minimizers (default=21)" << endl <<
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
void FilterReads_ParseOptions(int argc, char **argv, FilterReads_ProgramOptions& opt) {
  const char *opt_string = "vt:k:g:n:N:o:b:B:s:c:";
  static struct option long_options[] = {
    {"verbose",     no_argument,       0, 'v'},
    {"threads",     required_argument, 0, 't'},
    {"chunk-size",  required_argument, 0, 'c'},
    {"kmer-size",   required_argument, 0, 'k'},
    {"min-size",     no_argument,      0, 'g'},
    {"num-kmers",   required_argument, 0, 'n'},
    {"num-kmers2",  required_argument, 0, 'N'},
    {"output",      required_argument, 0, 'o'},
    {"bloom-bits",  required_argument, 0, 'b'},
    {"bloom-bits2", required_argument, 0, 'B'},
    {"seed",        required_argument, 0, 's'},
    {"ref",         no_argument,       0,  0 },
    {0,             0,                 0,  0 }
  };

  int option_index = 0, c;
  stringstream ss, ss2;
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
    switch (c) {
    case 0:
      if (strcmp(long_options[option_index].name, "ref") == 0) {
        opt.ref = true;
      }
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
    case 'g':
      opt.g = atoi(optarg);
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
bool FilterReads_CheckOptions(FilterReads_ProgramOptions& opt) {
  bool ret = true;

  size_t max_threads = std::thread::hardware_concurrency();

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
    /*cerr << "Setting chunksize to 1 because of only 1 thread" << endl;
      opt.read_chunksize = 1;*/
  }

  if (opt.k <= 0 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error, invalid value for kmer-size: " << opt.k << endl;
    cerr << "Values must be between 1 and " << (MAX_KMER_SIZE-1) << endl;
    ret = false;
  }

  if (opt.g <= 0 || opt.g >= MAX_KMER_SIZE) {
    cerr << "Error, invalid value for min-size: " << opt.g << endl;
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

  if (opt.ref) {
    opt.bf2 = 0;
    opt.nkmers2 = 0;
  }

  return ret;
}


// use:  FilterReads_PrintSummary(opt);
// pre:  opt has information about Kmer size, Bloom Filter sizes
// post: Information about the two Bloom Filters has been printed to cerr
void FilterReads_PrintSummary(const FilterReads_ProgramOptions& opt) {
  double fp;
  cerr << "Kmer size: " << opt.k << endl
       << "Chunksize: " << opt.read_chunksize << endl
       << "Using bloom filter size: " << opt.bf << " bits per element" << endl
       << "Estimated false positive rate: ";
  fp = pow(pow(.5,log(2.0)),(double) opt.bf);
  cerr << fp << endl;

  if (!opt.ref) {
    cerr << "Using bloom filter size for second set: " << opt.bf2 << " bits per element" << endl;
    cerr << "Estimated false positive rate for second set: ";
    fp = pow(pow(.5,log(2.0)),(double) opt.bf2);
    cerr << fp << endl;
  }
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
/*void FilterReads_Normal(const FilterReads_ProgramOptions& opt) {

    BlockedBloomFilter BF(opt.nkmers, (size_t) opt.bf);
    BlockedBloomFilter BF2(opt.nkmers2, (size_t) opt.bf2); // use different seeds

    bool done = false;
    char name[8192];
    string s;
    size_t name_len, len, read_chunksize = opt.read_chunksize;
    uint64_t n_read = 0;
    atomic<uint64_t> num_kmers(0), num_ins(0);

    FastqFile FQ(opt.files);
    vector<string> readv;

    const bool neighbor_hash = true;

    // Main worker thread
    auto worker_function = [&](vector<string>::const_iterator a, vector<string>::const_iterator b) {

        uint64_t l_num_kmers = 0, l_num_ins = 0;
        char* str;

        // for each input
        for (auto x = a; x != b; ++x) {

            str = const_cast<char*>(x->c_str());

            KmerHashIterator<RepHash> it_kmer_h(str, x->length(), opt.k), it_kmer_h_end;
            minHashIterator<RepHash> it_min(str, x->length(), opt.k, Minimizer::g, RepHash(), neighbor_hash);

            for (int last_pos = -1; it_kmer_h != it_kmer_h_end; ++it_kmer_h, ++it_min, ++l_num_kmers) {

                std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                // If one or more k-mer were jumped because contained non-ACGT char.
                if (p_.second != last_pos + 1){
                    str = &str[p_.second];
                    it_min = minHashIterator<RepHash>(str, x->length() - p_.second, opt.k, Minimizer::g, RepHash(), neighbor_hash);
                }

                last_pos = p_.second;
                uint64_t min_hr = it_min.getHash();

                if (!opt.ref) {

                    if (BF.search_and_insert(p_.first, min_hr)){
                        ++l_num_ins;
                        //BF.inc_count_block(min_hr, Minimizer(&str[it_min.getPosition()]), Kmer(&x->c_str()[p_.second]));
                    }
                    else if (BF2.search_and_insert(p_.first, min_hr)){
                        ++l_num_ins;
                        //BF2.inc_count_block(min_hr, Minimizer(&str[it_min.getPosition()]), Kmer(&x->c_str()[p_.second]));
                    }
                }
                else {

                    //BF.insert(p_.first, min_hr);
                    ++l_num_ins;
                }
            }
        }

        // atomic adds
        num_kmers += l_num_kmers;
        num_ins += l_num_ins;
    };

    while (!done) {

        readv.clear();
        size_t reads_now = 0;

        while (reads_now < read_chunksize) {

            if (FQ.read_next(name, &name_len, s, &len, NULL, NULL) >= 0) {

                readv.emplace_back(s);
                ++n_read;
                ++reads_now;
            }
            else {

                done = true;
                break;
            }
        }

        vector<thread> workers;
        // create worker threads
        auto rit = readv.begin();
        size_t batch_size = readv.size()/opt.threads;
        size_t leftover   = readv.size()%opt.threads;

        for (size_t i = 0; i < opt.threads; i++) {

            size_t jump = batch_size + ((i < leftover ) ? 1 : 0);

            auto rit_end(rit);
            advance(rit_end, jump);
            workers.push_back(thread(worker_function, rit, rit_end));

            rit = rit_end;
        }

        assert(rit==readv.end());

        for (auto &t : workers) t.join();
    }

    //BF.print_count_blocks(opt.output, "BF1");
    //BF2.print_count_blocks(opt.output, "BF2");

    FQ.close();

    //if (opt.verbose) {

        cerr << "Closed all fasta/fastq files" << endl;
        cerr << "processed " << num_kmers << " kmers in " << n_read  << " reads"<< endl;
        cerr << "found " << num_ins << " non-filtered kmers" << endl;
        cerr << "Writing bloom filter to " << opt.output << endl << "Bloom filter size is " << num_ins << endl;
    //}

    // First write metadata for bloom filter to opt.outputfile then the actual filter
    if (!opt.ref) {

        if (!BF2.WriteBloomFilter(opt.outputfile)) cerr << "Error writing data to file: " << opt.output << endl;
        else if (!BF.WriteBloomFilter(opt.outputfile)) cerr << "Error writing data to file: " << opt.output << endl;
    }

    fclose(opt.outputfile);
}*/

void FilterReads_Normal(const FilterReads_ProgramOptions& opt) {
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

    BlockedBloomFilter BF(opt.nkmers, (size_t) opt.bf);
    BlockedBloomFilter BF2(opt.nkmers2, (size_t) opt.bf2); // use different seeds

    bool done = false;
    const bool multi_threaded = (opt.threads != 1);
    char name[8192];
    string s;
    size_t name_len, len, read_chunksize = opt.read_chunksize;
    uint64_t n_read = 0;
    atomic<uint64_t> num_kmers(0), num_ins(0);

    FastqFile FQ(opt.files);
    vector<string> readv;

    // Main worker thread
    auto worker_function = [&](vector<string>::iterator a, vector<string>::iterator b) {

        uint64_t l_num_kmers = 0, l_num_ins = 0;

        // for each input
        for (auto x = a; x != b; ++x) {

            const char* str = x->c_str();
            int len = x->length();

            KmerHashIterator<RepHash> it_kmer_h(str, len, opt.k), it_kmer_h_end;
            minHashIterator<RepHash> it_min(str, len, opt.k, Minimizer::g, RepHash(), true);

            for (int last_pos = -1; it_kmer_h != it_kmer_h_end; ++it_kmer_h, ++it_min, ++l_num_kmers) {

                std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                // If one or more k-mer were jumped because contained non-ACGT char.
                if (p_.second != last_pos + 1)
                    it_min = minHashIterator<RepHash>(&str[p_.second], len - p_.second, opt.k, Minimizer::g, RepHash(), true);

                last_pos = p_.second;
                uint64_t min_hr = it_min.getHash();

                if (!opt.ref) {

                    if (BF.search_and_insert(p_.first, min_hr, multi_threaded)) ++l_num_ins;
                    else if (BF2.search_and_insert(p_.first, min_hr, multi_threaded)) ++l_num_ins;
                }
                else {

                    //BF.insert(p_.first, min_hr);
                    ++l_num_ins;
                }
            }
        }

        // atomic adds
        num_kmers += l_num_kmers;
        num_ins += l_num_ins;
    };

    while (!done) {

        readv.clear();
        size_t reads_now = 0;

        while (reads_now < read_chunksize) {

            if (FQ.read_next(name, &name_len, s, &len, NULL, NULL) >= 0) {

                readv.emplace_back(s);
                ++n_read;
                ++reads_now;
            }
            else {

                done = true;
                break;
            }
        }

        vector<thread> workers;
        // create worker threads
        auto rit = readv.begin();
        size_t batch_size = readv.size()/opt.threads;
        size_t leftover   = readv.size()%opt.threads;

        for (size_t i = 0; i < opt.threads; i++) {

            size_t jump = batch_size + ((i < leftover ) ? 1 : 0);

            auto rit_end(rit);
            advance(rit_end, jump);
            workers.push_back(thread(worker_function, rit, rit_end));

            rit = rit_end;
        }

        assert(rit==readv.end());

        for (auto &t : workers) t.join();
    }

    FQ.close();

    //if (opt.verbose) {

        cerr << "Closed all fasta/fastq files" << endl;
        cerr << "processed " << num_kmers << " kmers in " << n_read  << " reads"<< endl;
        cerr << "found " << num_ins << " non-filtered kmers" << endl;
        cerr << "Writing bloom filter to " << opt.output << endl;
        cerr << "Number of blocks in Bloom filter is " << BF2.getNbBlocks() << endl;
    //}

    // First write metadata for bloom filter to opt.outputfile then the actual filter
    if (!opt.ref) {

        if (!BF2.WriteBloomFilter(opt.outputfile)) cerr << "Error writing data to file: " << opt.output << endl;
        else if (!BF.WriteBloomFilter(opt.outputfile)) cerr << "Error writing data to file: " << opt.output << endl;
    }

    fclose(opt.outputfile);
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
    Minimizer::set_g(opt.g);

    if (opt.verbose) FilterReads_PrintSummary(opt);

    FilterReads_Normal(opt);
}
