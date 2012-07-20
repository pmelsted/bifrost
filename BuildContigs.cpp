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
  bool verbose;
  size_t threads, k;
  string freads, output, contigfilename, graphfilename;
  FILE *contigfile, *graphfile;
  size_t stride;
  bool stride_set;
  size_t read_chunksize;
  size_t contig_size; // not configurable
  vector<string> files;
  BuildContigs_ProgramOptions() : verbose(false), threads(1), k(0), stride(0), stride_set(false), \
                                  contigfile(NULL), graphfile(NULL), read_chunksize(1000), contig_size(1000000) {}
};

// use:  BuildContigs_PrintUsage();
// pre:   
// post: Information about how to "build contigs" has been inted to cerr
void BuildContigs_PrintUsage() {
  cerr << endl << "BFGraph " << BFG_VERSION << endl;
  cerr << "Creates contigs from filtered fasta/fastq files and saves results" << endl << endl;
  cerr << "Usage: BFGraph contigs [options] ... FASTQ files";
  cerr << endl << endl << "Options:" << endl <<
      "  -v, --verbose               Print lots of messages during run" << endl <<
      "  -t, --threads=INT           Number of threads to use (default 1)" << endl << 
      "  -c, --chunk-size=INT        Read chunksize to split betweeen threads (default 1000 for multithreaded else 1)" << endl <<
      "  -k, --kmer-size=INT         Size of k-mers, at most " << (int) (Kmer::MAX_K-1)<< endl << 
      "  -f, --filtered=STRING       File with filtered reads" << endl <<
      "  -o, --output=STRING         Prefix for output files" << endl <<
      "  -s, --stride=INT            Distance between saved kmers when mapping (default is kmer-size)"
  << endl << endl;
}


// use:  BuildContigs_ParseOptions(argc, argv, opt);
// pre:  argc is the parameter count, argv is a list of valid parameters for 
//       "building contigs" and opt is ready to contain the parsed parameters
// post: All the parameters from argv have been parsed into opt
void BuildContigs_ParseOptions(int argc, char **argv, BuildContigs_ProgramOptions &opt) {
  const char* opt_string = "vt:k:f:o:c:s:";
  static struct option long_options[] = {
      {"verbose",    no_argument,       0, 'v'},
      {"threads",    required_argument, 0, 't'},
      {"kmer-size",  required_argument, 0, 'k'},
      {"filtered",   required_argument, 0, 'f'},
      {"output",     required_argument, 0, 'o'},
      {"chunk-size", required_argument, 0, 'c'},
      {"stride",     required_argument, 0, 's'},
      {0,            0,                 0,  0 }
  };

  int option_index = 0, c;
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
      case 'f':
        opt.freads = optarg;
        break;
      case 'o': 
        opt.output = optarg;
        break;
      case 'c':
        opt.read_chunksize = atoi(optarg);
        break;
      case 's':
        opt.stride = atoi(optarg);
        opt.stride_set = true;
        break;
      default: break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  while (optind < argc) {
    opt.files.push_back(argv[optind++]);
  }
}


// use:  b = BuildContigs_CheckOptions(opt);
// pre:  opt contains parameters for "building contigs"
// post: (b == true)  <==>  the parameters are valid
bool BuildContigs_CheckOptions(BuildContigs_ProgramOptions &opt) {
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

  if (opt.k == 0 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error: Invalid kmer-size: " << opt.k 
         << ", need a number between 1 and " << (MAX_KMER_SIZE-1) << endl;
    ret = false;
  }
  
  if (opt.freads.empty()) {
    cerr << "Error: File with filtered reads missing" << endl;
  } else {
    struct stat freadsFileInfo;
    if (stat(opt.freads.c_str(), &freadsFileInfo) != 0) {
      cerr << "Error: File not found " << opt.freads << endl;
      ret = false;
    }
  }
  
  opt.contigfilename = opt.output + ".contigs";
  opt.graphfilename = opt.output + ".graph";
  opt.contigfile = fopen(opt.contigfilename.c_str(), "w");
  opt.graphfile = fopen(opt.graphfilename.c_str(), "w"); 
  
  if (opt.contigfile == NULL) {                                                                                                                  
    cerr << "Error: Could not open file for writing, " << opt.contigfile << endl;
    ret = false;
  } 
  if (opt.graphfile == NULL) {
    cerr << "Error: Could not open file for writing, " << opt.graphfile << endl;
    ret = false;
  } 

  
  if (opt.stride_set) {
    if (opt.stride == 0) {
      cerr << "Error: Invalid stride " << opt.stride << ", need a number >= 1" << endl;
      ret = false;
    }
  } else {
    opt.stride = opt.k;
  } 

  if (opt.files.size() == 0) {
    cerr << "Error: Missing fasta/fastq input files" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    vector<string>::const_iterator it;
    int intStat;
    for(it = opt.files.begin(); it != opt.files.end(); ++it) {
      intStat = stat(it->c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: File not found, " << *it << endl;
        ret = false;
      }
    }
  }
  
  return ret;
}


// use:  BuildContigs_PrintSummary(opt);
// pre:  opt has information about Kmer size, input file and output file
// post: Information about the Kmer size and the input and output files 
//       has been printed to cerr 
void BuildContigs_PrintSummary(const BuildContigs_ProgramOptions &opt) {
  cerr << "Kmer size: " << opt.k << endl
       << "Chunksize: " << opt.read_chunksize << endl
       << "Stride: " << opt.stride << endl
       << "Reading file with filtered reads: " << opt.freads << endl
       << "fasta/fastq files: " << endl;
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
  cerr << "Total:\t\t\t" << total << "MB" << endl << endl;
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

  size_t num_threads = opt.threads;
  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
  #endif

  BloomFilter bf;
  FILE *f = fopen(opt.freads.c_str(), "rb");
  if (f == NULL) {
    cerr << "Error, could not open file " << opt.freads << endl;
    exit(1);
  } 

  if (!bf.ReadBloomFilter(f)) {
    cerr << "Error reading bloom filter from file " << opt.freads << endl;
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

  size_t read_chunksize = opt.read_chunksize;
  int round = 0;
  cerr << "Starting real work ....." << endl << endl;

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
    ++round;
    if (read_chunksize > 1) {
      //cerr << "starting round " << round << endl;
    }

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

              getMappingInfo(cc.repequal, pos, cc.dist, kmernum, cmppos);
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
        /* We must use this when breaking on self-loops in find_contig_forward, 
           because then km is not neccessarily the first or the last km then in the contig */
        CheckContig cc = check_contig(bf, mapper, km);
        
        if(cc.cr.isEmpty()) {
          // The contig has not been mapped so we map it and increase coverage
          // of the kmers that came from the read
          
          size_t id = mapper.addContig(seq); 

          contig = mapper.getContig(id).ref.contig;
          size_t limit = it->end;
          contig->cover(it->start,limit);
        } else {
          // The contig has been mapped so we only increase the coverage of the
          // kmers that came from the read
          contig = mapper.getContig(cc.cr).ref.contig;
          int covlength = contig->numKmers();

          int32_t pos = cc.cr.ref.idpos.pos;
          bool repequal = (km == km.rep());

          getMappingInfo(repequal, pos, cc.dist, kmernum, cmppos); // 
          bool reversed = ((pos >= 0) != repequal);
          size_t start = it->start, end = it->end;
          if (reversed) {
            assert(contig->seq.getKmer(kmernum) == km.twin());
            kmernum -= it->start;
            int left = kmernum - (end - start);
            size_t right = kmernum;
            /* We could use % covlength in the cover function instead of this if-else catastrophe */
            if (left < 0) {
              // Maps to a self-looping contig
              // We have to update coverage "out of bounds" (out from the right continues at the beginning)
              assert(0 <= left + covlength);
              assert(cc.dist > 0); 
              contig->cover(0, right);
              contig->cover(left + covlength, covlength - 1);
            } else {
              contig->cover(left, right);
            }
          }
          else {
            assert(contig->seq.getKmer(kmernum) == km);
            kmernum += it->start;
            size_t left = kmernum;
            size_t right = kmernum + end - start;
            /* We could use % covlength in the cover function instead of this if-else catastrophe */
            if (right >= covlength) {
              // Maps to a self-looping contig
              // We have to update coverage "out of bounds" (out from the right continues at the beginning)
              assert(cc.dist > 0); 
              assert(right < 2*covlength);
              contig->cover(0, right - covlength);
              contig->cover(left, covlength - 1);
            } else {
              contig->cover(left, right);
            }
          }
        }
      }
      parray[i].clear();

    }
    
    if (read_chunksize > 1) {
      cerr << " end of round" << endl;
      cerr << " processed " << mapper.contigCount() << " contigs" << endl;
    }
  }
  FQ.close();
  cerr << "Closed all fasta/fastq files" << endl;

  size_t contigsBefore = mapper.contigCount();
  cerr << "Splitting and joining the contigs" << endl;
  pair<pair<size_t, size_t>, size_t> contigDiff = mapper.splitAndJoinContigs();
  int contigsAfter = contigsBefore + contigDiff.first.first - contigDiff.first.second - contigDiff.second;
  if (opt.verbose) {
    cerr << "Before split and join: " << contigsBefore << " contigs" << endl;
    cerr << "After split and join: " << contigsAfter << " contigs" <<  endl;
    cerr << "Contigs splitted: " << contigDiff.first.first << endl;
    cerr << "Contigs deleted: " << contigDiff.first.second << endl;
    cerr << "Contigs joined: " << contigDiff.second << endl;
    cerr << "Number of reads " << n_read  << ", kmers stored " << mapper.size() << endl << endl;
    printMemoryUsage(bf, mapper);
  }

  mapper.writeContigs(opt.contigfile, opt.graphfile);
  cerr << "Writing contigs to file: " << opt.contigfilename << endl
       << "Writing the graph to file: " << opt.graphfilename << endl;

  fclose(opt.contigfile);
  fclose(opt.graphfile);
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
