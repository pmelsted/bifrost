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
#include <utility>
#include <vector>

#include "Common.hpp"
#include "CompressedSequence.hpp"
#include "Contig.hpp"
#include "KmerMapper.hpp"
#include "BloomFilter.hpp"
#include "ContigMethods.hpp"
#include "HashTables.hpp"
#include "KmerIterator.hpp"
#include "fastq.hpp"


struct BuildContigs_ProgramOptions {
  size_t k;
  size_t read_chunksize;
  size_t contig_size; // not configurable
  string input;
  string output;
  bool verbose;
  size_t threads;
  size_t stride;
  bool stride_set;
  vector<string> files;
  BuildContigs_ProgramOptions() : k(0), verbose(false) , contig_size(1000000), read_chunksize(1000), \
                                  threads(1), stride(0), stride_set(false) {}
};


// use:  BuildContigs_PrintUsage();
// pre:   
// post: Information about how to "build contigs" has been inted to cerr
void BuildContigs_PrintUsage() {
  cerr << "BFGraph " << BFG_VERSION << endl << endl;
  cerr << "Filters errors in fastq or fasta files and saves results" << endl << endl;
  cerr << "Usage: BFGraph contigs [options] ... FASTQ files";
  cerr << endl << endl <<
      "-k, --kmer-size=INT             Size of k-mers, at most " << (int) (Kmer::MAX_K-1)<< endl << 
      "-s, --stride=INT                Distance between saved kmers when mapping (default is the kmer size)" << endl << 
      "-c, --chunk-size=INT            Read chunksize to split betweeen threads (default 1000)" << endl <<
      "-t, --threads=INT               Number of threads to use (default 1)" << endl << 
      "-i, --input=STRING              Filtered reads" << endl <<
      "-o, --output=STRING             Filename for output" << endl <<
      "    --verbose                   Print lots of messages during run" << endl << endl
      ;
}


// use:  BuildContigs_ParseOptions(argc, argv, opt);
// pre:  argc is the parameter count, argv is a list of valid parameters for 
//       "building contigs" and opt is ready to contain the parsed parameters
// post: All the parameters from argv have been parsed into opt
void BuildContigs_ParseOptions(int argc, char **argv, BuildContigs_ProgramOptions &opt) {
  int verbose_flag = 0;
  const char* opt_string = "k:o:i:c:t:s:";
  static struct option long_options[] =
      {
        {"verbose", no_argument,  &verbose_flag, 1},
        {"kmer-size", required_argument, 0, 'k'},
        {"stride", optional_argument, 0, 's'},
        {"chunk-size", optional_argument, 0, 'c'},
        {"threads", optional_argument, 0, 't'},
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
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
      case 'c':
        opt.read_chunksize = atoi(optarg);
        break;
      case 'k': 
        opt.k = atoi(optarg); 
        break;
      case 'o': 
        opt.output = optarg;
        break;
      case 'i':
        opt.input = optarg;
        break;
      case 't':
        opt.threads = atoi(optarg);
        break;
      case 's':
        opt.stride = atoi(optarg);
        opt.stride_set = true;
        break;
      default: break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; ++i) {
    opt.files.push_back(argv[i]);
  }
  
  if (verbose_flag) {
    opt.verbose = true;
  }
}


// use:  b = BuildContigs_CheckOptions(opt);
// pre:  opt contains parameters for "building contigs"
// post: (b == true)  <==>  the parameters are valid
bool BuildContigs_CheckOptions(BuildContigs_ProgramOptions &opt) {
  bool ret = true;

  if (opt.k <= 0 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error, invalid value for kmer-size: " << opt.k << endl;
    cerr << "Values must be between 1 and " << (MAX_KMER_SIZE-1) << endl;
    ret = false;
  }

  if (opt.read_chunksize <= 0) {
    cerr << "Error, invalid value for chunk-size: " << opt.read_chunksize << endl;
    cerr << "Value must be greater than 0" << endl;
    ret = false;
  }
  
  if (opt.input.empty()) {
    cerr << "Input file missing" << endl;
  } else {
    struct stat inFileInfo;
    if (stat(opt.input.c_str(), &inFileInfo) != 0) {
      cerr << "Error input file not found " << opt.input << endl;
      ret = false;
    }
  }

  if (opt.threads < 1) {
    cerr << "Invalid number of threads " << opt.threads << ", need a number >= 1" << endl;
    ret = false;
  }
  
  if (opt.stride_set) {
    if (opt.stride < 1) {
    cerr << "Invalid value for stride " << opt.threads << ", need a number >= 1" << endl;
    ret = false;
    }
  } else {
    opt.stride = opt.k;
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


// use:  BuildContigs_PrintSummary(opt);
// pre:  opt has information about Kmer size, input file and output file
// post: Information about the Kmer size and the input and output files 
//       has been printed to cerr 
void BuildContigs_PrintSummary(const BuildContigs_ProgramOptions &opt) {
  cerr << "Kmer size " << opt.k << endl
       << "Chunksize " << opt.read_chunksize << endl
       << "Reading input file " << opt.input << endl
       << "input files: " << endl;
  vector<string>::const_iterator it;
  for (it = opt.files.begin(); it != opt.files.end(); ++it) {
    cerr << "  " << *it << endl;
  }
}

void printMemoryUsage(BloomFilter &bf, KmerMapper &mapper) {
  size_t total = 0;
  cerr << "   -----  Memory usage  -----   " << endl;
  total += bf.memory();
  total += mapper.memory();
  total >>= 20;
  cerr << "Total:\t\t\t" << total << "MB" << endl;
}

// use:  BuildContigs_Normal(opt);
// pre:  opt has information about Kmer size, input file and output file
// post: The contigs have been written to the output file 
void BuildContigs_Normal(const BuildContigs_ProgramOptions &opt) {
  /**
   *  outline of algorithm
   *   - open bloom filter file
   *   - create contig datastructures
   *   - for each read
   *     - for all kmers in read
   *        - if kmer is in bf
   *          - if it maps to contig
   *            - try to jump over as many kmers as possible
   *          - else
   *            - create new contig from kmers in both directions
   *            - from this kmer while there is only one possible next kmer
   *            - with respect to the bloom filter
   *            - when the contig is ready, 
   *            - try to jump over as many kmers as possible
   */

  size_t num_threads = opt.threads, max_threads = 1;
  
  #pragma omp parallel
  {
    #pragma omp master 
    {
      #ifdef _OPENMP
        max_threads = omp_get_num_threads(); 
      #endif 
    }
  }
  
  if (num_threads > max_threads) {
    cerr << "Using " << max_threads << " thread(s) instead of " << num_threads << " due to number of cores" << endl;
    num_threads = max_threads;
  } else {
    cerr << "Using " << num_threads << " thread(s)" << endl;
  }


  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
  #endif


  BloomFilter bf;
  FILE* f = fopen(opt.input.c_str(), "rb");
  if (f == NULL) {
    cerr << "Error, could not open file " << opt.input << endl;
    exit(1);
  } 

  if (!bf.ReadBloomFilter(f)) {
    cerr << "Error reading bloom filter from file " << opt.input << endl;
    fclose(f);
    f = NULL;
    exit(1);
  } else {
    fclose(f);
    f = NULL;
  }

  KmerMapper mapper(opt.contig_size, opt.stride);

  KmerIterator iter, iterend;
  FastqFile FQ(opt.files);

  char name[8192], s[8192];
  size_t kmernum, name_len, len, k = Kmer::k;
  int32_t cmppos;
  uint64_t n_read = 0; 

  vector<string> readv;
  vector<NewContig> *smallv, *parray = new vector<NewContig>[num_threads];
  bool done = false;
  size_t reads_now, read_chunksize = opt.read_chunksize;
  cerr << "starting real work" << endl;
  int round = 0;
  while (!done) {
    readv.clear();
    reads_now = 0;

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
    ++round;
    cerr << "starting round " << round << endl;

    #pragma omp parallel default(shared) private(kmernum,cmppos,smallv) shared(mapper,parray,readv,bf,reads_now,k)
    {
      KmerIterator iter, iterend;
      size_t threadnum = 0;
      #ifdef _OPENMP
        threadnum= omp_get_thread_num();
      #endif
      smallv = &parray[threadnum];


      #pragma omp for nowait
      for(size_t index=0; index < reads_now; ++index) {
        const char *cstr = readv[index].c_str();
        iter = KmerIterator(cstr);
        Kmer km, rep;

        if (iter != iterend) {
          km = iter->first;
          rep = km.rep();
        }

        while (iter != iterend) {
          if (!bf.contains(rep)) { // km is not in the graph
            // jump over it
            iter.raise(km, rep);
          } else {
            CheckContig cc = check_contig(bf, mapper, km);
            if (cc.cr.isEmpty()) {
              // Map contig (or increase coverage) from this sequence after this thread finishes
              MakeContig mc = make_contig(bf, mapper, km);
              string seq = mc.seq;
              
              size_t index = iter->second;
              size_t cmpindex = index + k;
              size_t seqindex = mc.pos;
              size_t seqcmpindex = seqindex + k;
              iter.raise(km, rep);
              
              while (cstr[cmpindex] == seq[seqcmpindex] && cstr[cmpindex] != '\0') {
                assert(cstr[cmpindex] != 'N');
                Kmer k1(&cstr[cmpindex-k]);
                Kmer k2(&seq[seqcmpindex-k]);
                assert(k1 == k2);
                ++cmpindex;
                ++seqcmpindex;
                iter.raise(km,rep);
              }

              // cstr[index,...,cmpindex-1] == seq[seqindex,...,seqcmpindex-1]
              // Coverage of seq[seqindex,...,seqcmpindex-k] will by increase by one later in master thread
              smallv->push_back(NewContig(seq,seqindex,seqcmpindex-k));
            } else {
              Contig *contig = mapper.getContig(cc.cr).ref.contig;
              int32_t pos = cc.cr.ref.idpos.pos;

              getMappingInfo(cc.repequal, pos, cc.dist, k, kmernum, cmppos);
              bool reversed = (pos >= 0) != cc.repequal;
              int jumpi = 1 + iter->second + contig->seq.jump(cstr, iter->second + k, cmppos, reversed);

              if (reversed) {
                assert(contig->seq.getKmer(kmernum) == km.twin());
              } else {
                assert(contig->seq.getKmer(kmernum) == km);
              }
              contig->cover(kmernum,kmernum);

              int32_t direction = reversed ? -1 : 1;
              kmernum += direction;
              iter.raise(km, rep);

              while (iter != iterend && iter->second < jumpi) {
                assert(cstr[iter->second+k-1] != 'N');
                assert(kmernum >= 0);
                assert(kmernum < contig->numKmers());

                // Update coverage
                contig->cover(kmernum,kmernum);

                kmernum += direction;
                iter.raise(km, rep);
              }
            }
          }     
        }
      }
    }

    for (size_t i=0; i < num_threads; i++) {
      for (vector<NewContig>::iterator it=parray[i].begin(); it != parray[i].end(); ++it) {
        // The kmer did not map when it was added to this vector
        // so we make the contig if it has not been made yet
        const char *seq = it->seq.c_str();
        Contig *contig;

        Kmer km(seq); 
        ContigRef mapcr = mapper.find(km);
        
        if(mapcr.isEmpty()) {
          // The contig has not been mapped so we map it and increase coverage
          // of the kmers that came from the read
          
          size_t id = mapper.addContig(seq); 

          contig = mapper.getContig(id).ref.contig;
          size_t limit = it->end;
          contig->cover(it->start,limit);
        } else {
          // The contig has been mapped so we only increase the coverage of the
          // kmers that came from the read
          contig = mapper.getContig(mapcr).ref.contig;

          int32_t pos = mapcr.ref.idpos.pos;
          bool repequal = (km == km.rep());

          getMappingInfo(repequal, pos, 0, k, kmernum, cmppos); // 0 because km is the first or the last kmer in the contig
          bool reversed = ((pos >= 0) != repequal);
          size_t start = it->start, end = it->end;
          if (reversed) {
            assert(contig->seq.getKmer(kmernum) == km.twin());
            kmernum -= it->start;
            contig->cover(kmernum - (end - start), kmernum);
          }
          else {
            assert(contig->seq.getKmer(kmernum) == km);
            kmernum += it->start;
            contig->cover(kmernum, kmernum + end - start);
          }
        }
      }
      parray[i].clear();

    }
    
    cerr << " end of round" << endl;
    cerr << " processed " << mapper.size() << " contigs" << endl;
  }

  size_t contigsBefore = mapper.contigCount();
  pair<pair<size_t, size_t>, size_t> contigDiff = mapper.splitAndJoinContigs();
  int contigsAfter = contigsBefore + contigDiff.first.first - contigDiff.first.second - contigDiff.second;
  if (opt.verbose) {
    cerr << "Before split and join: " << contigsBefore << " contigs" << endl;
    cerr << "After split and join: " << contigsAfter << " contigs" <<  endl;
    cerr << "Contigs splitted: " << contigDiff.first.first << endl;
    cerr << "Contigs deleted: " << contigDiff.first.second << endl;
    cerr << "Contigs joined: " << contigDiff.second << endl;
    cerr << "Number of reads " << n_read  << ", kmers stored " << mapper.size() << endl;
    printMemoryUsage(bf, mapper);
  }
  mapper.writeContigs(opt.output);
  delete [] parray;
}



// use:  BuildContigs(argc, argv);
// pre:  argc is the number of arguments in argv and argv includes 
//       arguments for "building the contigs", including filenames
// post: If the number of arguments is correct and the arguments are valid
//       the "contigs have been built" and written to a file
void BuildContigs(int argc, char** argv) {
  
  BuildContigs_ProgramOptions opt;

  BuildContigs_ParseOptions(argc,argv,opt);

  if (argc < 2) {
    BuildContigs_PrintUsage();
    exit(1);
  }
  
  if (!BuildContigs_CheckOptions(opt)) {
    BuildContigs_PrintUsage();
    exit(1);
  }
  
  // set static global k-value
  Kmer::set_k(opt.k);

  if (opt.verbose) {
    BuildContigs_PrintSummary(opt);
  }

  BuildContigs_Normal(opt);

}
