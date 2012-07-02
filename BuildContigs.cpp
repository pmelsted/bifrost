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

#include "boost/tuple/tuple.hpp"

#include "BloomFilter.hpp"
#include "Common.hpp"
#include "CompressedSequence.hpp"
#include "Contig.hpp"
#include "HashTables.hpp"
#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "KmerMapper.hpp"
#include "fastq.hpp"

using boost::tuples::tie;

struct NewContig {
  string seq;
  size_t start;
  size_t end;
  NewContig(string s, size_t a, size_t b) : seq(s), start(a), end(b) {}
};
  

pair<Kmer, size_t> find_contig_forward(BloomFilter &bf, Kmer km, string* s);
pair<ContigRef, pair<size_t, bool> > check_contig(BloomFilter &bf, KmerMapper &mapper, Kmer km);
pair<string, size_t> make_contig(BloomFilter &bf, KmerMapper &mapper, Kmer km);
void getMappingInfo(const bool repequal, const int32_t pos, const size_t dist, const size_t k, size_t &kmernum, int32_t &cmppos);


static const char alpha[4] = {'A','C','G','T'};
static const char beta[4] = {'T','G','A','C'}; // c -> beta[(c & 7) >> 1] maps: 'A' <-> 'T', 'C' <-> 'G'

struct BuildContigs_ProgramOptions {
  size_t k;
  size_t read_chunksize;
  size_t contig_size; // not configurable
  string input;
  string output;
  bool verbose;
  size_t threads;
  vector<string> files;
  BuildContigs_ProgramOptions() : k(0), verbose(false) , contig_size(1000000), read_chunksize(0), threads(1) {}
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
      "-c, --chunk-size=INT            Read chunksize to split betweeen threads" << endl <<
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
  const char* opt_string = "k:o:i:c:t:";
  static struct option long_options[] =
      {
        {"verbose", no_argument,  &verbose_flag, 1},
        {"kmer-size", required_argument, 0, 'k'},
        {"chunk-size", required_argument, 0, 'c'},
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 't'},
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
  cerr << "Kmer size" << opt.k << endl
       << "Reading input file " << opt.input << endl
       << "Writing to output " << opt.output << endl
       << "input files: " << endl;
  vector<string>::const_iterator it;
  for (it = opt.files.begin(); it != opt.files.end(); ++it) {
    cerr << "  " << *it << endl;
  }
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

  size_t num_threads = 1;
  #ifdef _OPENMP
    omp_set_num_threads(opt.threads);
  #endif
  
  #pragma omp parallel
  {
    #pragma omp master 
    {
      #ifdef _OPENMP
        num_threads = omp_get_num_threads(); 
      #endif 
    }
  }


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

  KmerMapper mapper(opt.contig_size);

  KmerIterator iter, iterend;
  FastqFile FQ(opt.files);
  ContigRef cr, mapcr;

  bool repequal, reversed;
  char name[8192], s[8192];
  size_t kmernum, dist, name_len, len, k = Kmer::k;
  int32_t pos, cmppos;
  uint64_t n_read = 0; 
  pair<size_t, bool> disteq;

  vector<string> readv;
  vector<NewContig> *smallv, *parray = new vector<NewContig>[num_threads];
  bool done = false;
  size_t reads_now, read_chunksize = opt.read_chunksize;
  cerr << "using chunksize " << read_chunksize << endl;
  cerr << "starting real work" << endl;

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

    #pragma omp parallel default(shared) private(dist,mapcr,disteq,kmernum,cmppos,smallv,repequal) shared(mapper,parray,readv,bf,reads_now,k)
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
            tie(mapcr, disteq) = check_contig(bf, mapper, km);
            if (mapcr.isEmpty()) {
              // Map contig (or increase coverage) from this sequence after this thread finishes
              pair<string, size_t> cpair = make_contig(bf, mapper, km);
              string seq = cpair.first;
              
              size_t index = iter->second;
              size_t cmpindex = index + k;
              size_t seqindex = cpair.second;
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
              // increase coverage of cov[seqindex,...,seqcmpindex-k] by one, later in master thread
              smallv->push_back(NewContig(seq,seqindex,seqcmpindex-k));
            } else {
              Contig *contig = mapper.getContig(mapcr).ref.contig;
              tie(dist, repequal) = disteq;
              int32_t pos = mapcr.ref.idpos.pos;

              getMappingInfo(repequal, pos, dist, k, kmernum, cmppos);
              bool reversed = (pos >= 0) != repequal;
              int jumpi = 1 + iter->second + contig->seq.jump(cstr, iter->second + k, cmppos, reversed);

              if (reversed) {
                assert(contig->seq.getKmer(kmernum) == km.twin());
              } else {
                assert(contig->seq.getKmer(kmernum) == km);
              }
              /*
              uint8_t *change = &contig->cov[kmernum];
              uint8_t oldval = *change; 
              while (1) {
                if (oldval < 0xff) {
                  if(__sync_bool_compare_and_swap(change, oldval, oldval +1)) {
                    break;
                  }
                } else {
                  break;
                }
                oldval = *change; 
              }
              */
              // The new CompressedCoverage class
              contig->covp->cover(kmernum,kmernum);

              int32_t direction = reversed ? -1 : 1;
              kmernum += direction;
              iter.raise(km, rep);

              while (iter != iterend && iter->second < jumpi) {
                assert(cstr[iter->second+k-1] != 'N');
                assert(kmernum >= 0);
                assert(kmernum < contig->covlength);
                /*
                uint8_t *change = &contig->cov[kmernum];
                uint8_t oldval = *change; 
                while (1) {
                  if (oldval < 0xff) {
                    if (__sync_bool_compare_and_swap(change, oldval, oldval +1)) {
                      break;
                    }
                  } else {
                    break;
                  }
                  oldval = *change; 
                }
                */
                // The new CompressedCoverage class
                contig->covp->cover(kmernum,kmernum);

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

        Kmer km(seq); // Is this definitely the kmer we want to check? 
        tie(mapcr, disteq) = check_contig(bf, mapper, km);
        
        if(mapcr.isEmpty()) {
          // The contig has not been mapped so we map it and increase coverage
          // of the kmers that came from the read
          size_t id = mapper.addContig(seq);
          contig = mapper.getContig(id).ref.contig;
          //tie(mapcr, disteq) = check_contig(bf, mapper, km);
          size_t limit = it->end;
          /*
          for (size_t index=it->start; index <= limit ; ++index) {
            Kmer covkm = Kmer(seq+index);
            assert(contig->seq.getKmer(index) == covkm);
            if (contig->cov[index] < 0xff) {
              contig->cov[index] += 1;
            }
          }
          */
          // The new CompressedCoverage class
          contig->covp->cover(it->start,limit);
        } else {
          // The contig has been mapped so we only increase the coverage of the
          // kmers that came from the read
          contig = mapper.getContig(mapcr).ref.contig;

          tie(dist, repequal) = disteq;
          assert(dist == 0); // km is the first or the last kmer in the contig
          pos = mapcr.ref.idpos.pos;

          getMappingInfo(repequal, pos, dist, k, kmernum, cmppos);
          reversed = ((pos >= 0) != repequal);
          int32_t direction = reversed ? -1 : 1;
          size_t start = it->start, end = it->end;
          if (reversed) {
            assert(contig->seq.getKmer(kmernum) == km.twin());
            kmernum -= it->start;
            // The new CompressedCoverage class
            contig->covp->cover(kmernum-(end-start),kmernum);
          }
          else {
            assert(contig->seq.getKmer(kmernum) == km);
            kmernum += it->start;
            // The new CompressedCoverage class
            contig->covp->cover(kmernum,kmernum+end-start);
          }
          /*
          while (start <= end) {
            if (contig->cov[kmernum] < 0xff) {
              contig->cov[kmernum] += 1;
            }
            kmernum += direction;
            ++start;
          }
          */
        }
      }
      parray[i].clear();
    }
  }
  // Print the good contigs
  size_t contigsBefore = mapper.contigCount();
  int contigDiff = mapper.splitAndJoinContigs();
  int contigsAfter = contigsBefore + contigDiff;
  mapper.writeContigs(opt.output);
  cerr << "Before split and join: " << contigsBefore << " contigs" << endl;
  cerr << "After split and join: " << contigsAfter << " contigs" <<  endl;
  cerr << "Number of reads " << n_read  << ", kmers stored " << mapper.size() << endl;
  cerr << "Used " << num_threads << " threads and chunksize " << read_chunksize << endl;
  delete [] parray;
}


// use:  getMappingInfo(repequal, pos, dist, k, kmernum, cmppos)
// pre:  
// post: cmppos is the first character after the kmer-match at position pos
void getMappingInfo(const bool repequal, const int32_t pos, const size_t dist, const size_t k, size_t &kmernum, int32_t &cmppos) {
  // Now we find the right location of the kmer inside the contig
  // to increase coverage 
  if (pos >= 0) {
    if (repequal) {
      cmppos = pos - dist + k;
      kmernum = cmppos - k;
    } else {
      cmppos = pos - 1 + dist;
      kmernum = cmppos +1;
    }
  } else {
    if (repequal) {
      cmppos = -pos + dist -k;
      kmernum = cmppos +1;
    } else {
      cmppos = -pos + 1 - dist; // Original: (-pos +1 -k) - dist + k
      kmernum = cmppos - k;
    }
  }
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



// use:  (cr, (dist, eq)) = check_contig_(bf,km,mapper);
// pre:  
// post: if km does not map to a contig: cr.isEmpty() == true and dist == 0
//       else: km is in a contig which cr maps to and dist is the distance 
//             from km to the mapping location 
//             (eq == true):  km has the same direction as the contig
//             else:  km has the opposite direction to the contig
pair<ContigRef, pair<size_t, bool> > check_contig(BloomFilter &bf, KmerMapper &mapper, Kmer km) {
  ContigRef cr = mapper.find(km);
  if (!cr.isEmpty()) {
    return make_pair(cr, make_pair(0, km == km.rep())); 
  }
  int i, j;
  size_t dist = 1;
  bool found = false;
  Kmer end = km;
  while (dist < mapper.stride) {
    size_t fw_count = 0;
    j = -1;
    for (i = 0; i < 4; ++i) {
      Kmer fw_rep = end.forwardBase(alpha[i]).rep();
      if (bf.contains(fw_rep)) {
        j = i;
        ++fw_count;
        if (fw_count > 1) {
          break;
        }
      }
    }

    if (fw_count != 1) {
      break;
    }

    Kmer fw = end.forwardBase(alpha[j]);

    size_t bw_count = 0;
    for (i = 0; i < 4; ++i) {
      Kmer bw_rep = fw.backwardBase(alpha[i]).rep();
      if (bf.contains(bw_rep)) {
        ++bw_count;
        if (bw_count > 1) {
          break;
        }
      }
    }

    assert(bw_count >= 1);
    if (bw_count != 1) {
      break;
    }
    cr = mapper.find(fw);
    end = fw;
    if (!cr.isEmpty()) {
      found = true;
      break;
    }
    ++dist;
  }
  if (found) {
    return make_pair(cr, make_pair(dist, end == end.rep())); 
  } else {
    return make_pair(ContigRef(), make_pair(0,0));
  }
}

// use:  seq, pos = make_contig(bf, mapper, km, s);
// pre:  km is not contained in a mapped contig in mapper 
// post: Finds the forward and backward limits of the contig
//       which contains km  according to the bloom filter bf and puts it into seq
//       pos is the position where km maps into this contig
pair<string, size_t> make_contig(BloomFilter &bf, KmerMapper &mapper, Kmer km) {
  size_t k  = Kmer::k;
  string seq, seq_fw(k, 0), seq_bw(k, 0);
  //pair<Kmer, size_t> p_fw = find_contig_forward(bf, km, &seq_fw);
  find_contig_forward(bf, km, &seq_fw);
  pair<Kmer, size_t> p_bw = find_contig_forward(bf, km.twin(), &seq_bw);
  ContigRef cr_tw_end = mapper.find(p_bw.first);
  assert(cr_tw_end.isEmpty());

  if (p_bw.second > 1) {
    seq.reserve(seq_bw.size() + seq_fw.size() - k);
    // copy reverse part of seq_bw not including k
    for (size_t j = seq_bw.size() - 1; j >= k; --j) {
      seq.push_back(beta[(seq_bw[j] & 7) >> 1]);
    }
    // append seq_fw
    seq += seq_fw;
  } else {
    seq = seq_fw;
  }
  return make_pair(seq, p_bw.second-1);
}

// use:  tie(end,dist) = find_contig_forward(bf,km,s);
// pre:  
// post: km is contained in a contig c with respect to the
//       bloom filter graph bf and end is the forward endpoint (wrt km direction)
//       and c contains dist kmers until the end (including km)
//       if s is not NULL the sequence of the contig is stored in s
pair<Kmer, size_t> find_contig_forward(BloomFilter &bf, Kmer km, string* s) {
  int j;
  size_t i,dist = 1;
  vector<char> v;

  Kmer end = km;

  assert(bf.contains(km.rep()));
  if (s != NULL) {
    char t[Kmer::MAX_K+1];
    km.toString(t);
    for (i = 0; i < Kmer::k; ++i) {
      v.push_back(t[i]);
    }
  }
  
  while (true) {
    assert(bf.contains(end.rep()));
    size_t fw_count = 0;
    j = -1;
    for (i = 0; i < 4; ++i) {
      Kmer fw_rep = end.forwardBase(alpha[i]).rep();
      if (bf.contains(fw_rep)) {
        j = i;
        ++fw_count;
        if (fw_count > 1) {
          break;
        }
      }
    }

    if (fw_count != 1) {
      break;
    }
    
    Kmer fw = end.forwardBase(alpha[j]);
    assert(0 <= j && j < 4);
    assert(bf.contains(fw.rep()));

    size_t bw_count = 0;
    for (i = 0; i < 4; ++i) {
      Kmer bw_rep = fw.backwardBase(alpha[i]).rep();
      if (bf.contains(bw_rep)) {
        ++bw_count;
        if (bw_count > 1) {
          break;
        }
      }
    }

    assert(bw_count >= 1);
    if (bw_count != 1) {
      break;
    }
    
    end = fw;
    ++dist;
    if (s != NULL) {
      v.push_back(alpha[j]);
    }
  }

  if (s != NULL) {
    s->clear();
    s->reserve(v.size());
    s->insert(s->begin(), v.begin(), v.end());
  }
  return make_pair(end, dist);
}
