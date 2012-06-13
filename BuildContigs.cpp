#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>
#include <utility>

#include <sstream>
#include <vector>
#include <string>

#include <stdint.h>
#include <sys/stat.h>
#include <functional>

#include <getopt.h>

#include "Common.hpp"

#include "HashTables.hpp"
#include "fastq.hpp"
#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "BloomFilter.hpp"
#include "CompressedSequence.hpp"
#include "KmerMapper.hpp"
#include "Contig.hpp"


pair<Kmer, size_t> find_contig_forward(BloomFilter &bf, Kmer km, string* s);
ContigRef kmer_maps_to_contig(BloomFilter &bf, Kmer km, KmerMapper &mapper);

struct BuildContigs_ProgramOptions {
  size_t k;
  size_t contig_size; // not configurable
  string input;
  string output;
  bool verbose;
  vector<string> files;
  BuildContigs_ProgramOptions() : k(0), verbose(false) , contig_size(1000000) {}
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
      "-i, --input=STRING            Filtered reads" << endl <<
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
  const char* opt_string = "k:o:i:";
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
      case 'i':
        opt.input = optarg;
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

  
  if (opt.input.empty()) {
    cerr << "Input file missing" << endl;
  } else {
    struct stat inFileInfo;
    if (stat(opt.input.c_str(), &inFileInfo) != 0) {
      cerr << "Error input file not found " << opt.input << endl;
      ret = false;
    }
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

  size_t k  = Kmer::k;
  char name[8196],s[8196];
  size_t name_len, len;

  uint64_t n_read = 0;
  uint64_t num_kmers = 0;  
  uint64_t num_ins = 0;
 
  FastqFile FQ(opt.files);
  
  Kmer km, rep;
  KmerIterator iter, iterend;

  cerr << "starting real work" << endl;

  Contig *contig;
  size_t i, jumpi;
  int32_t contiglen, pos, cmppos;
  bool repequal, reversed;

  while (FQ.read_next(name, &name_len, s, &len, NULL, NULL) >= 0) {
    iter = KmerIterator(s);
    n_read++;
    i = 0;

    if (iter != iterend) {
      km = iter->first;
      rep = km.rep();
    }

    while (iter != iterend) {
      if (!bf.contains(rep)) { // kmer i is not in the graph
        i++; // jump over it
        iter.raise(km, rep);
      } else {
        ContigRef cr = mapper.find(rep);

        if (cr.isEmpty()) {
          // The kmer does not map but the contig could although exist 

          pair<Kmer, size_t> p_fw,p_bw;
          p_fw = find_contig_forward(bf, km, NULL); 
          
          ContigRef cr_end = mapper.find(p_fw.first);

          // Check whether we have made this contig or not
          if (cr_end.isEmpty()) {
            assert(kmer_maps_to_contig(bf, km, mapper).isEmpty());
            // We have not made this contig, lets make it
            
            string seq, seq_fw(k,0), seq_bw(k,0);
            // Find the forward and backward limits of this contig
            // according to the bloom filter
            p_fw = find_contig_forward(bf, km, &seq_fw);
            p_bw = find_contig_forward(bf, km.twin(), &seq_bw);
            ContigRef cr_tw_end = mapper.find(p_bw.first);
            assert(cr_tw_end.isEmpty());

            if (p_bw.second > 1) {
              seq.reserve(seq_bw.size() + seq_fw.size() - k);
              // copy reverse part of seq_bw not including k
              // TODO: refactor this to util package
              for (int j = seq_bw.size()-1; j>=k; j--) {
                char c = seq_bw[j];
                char cc = 'N';
                switch (c) {
                  case 'A': cc = 'T'; break;
                  case 'C': cc = 'G'; break;
                  case 'G': cc = 'C'; break;
                  case 'T': cc = 'A'; break;
                }
                seq.push_back(cc);
              }
              // append seq_fw
              seq += seq_fw;
            } else {
              seq = seq_fw;
            }

            mapper.addContig(seq);
            
            //cerr << seq.size() << endl;
            ContigRef found = mapper.find(p_bw.first);
            assert(!found.isEmpty());
            if (!found.isEmpty()) {
              mapper.printContig(found.ref.idpos.id);
            }
            cr_end = mapper.find(p_fw.first);
          }

          assert(!kmer_maps_to_contig(bf, km, mapper).isEmpty());

          // cr_end is a contigRef pointing to the last position of the contig containing the kmer
	  
          // Now we jump as far ahead as we can
          pos = cr_end.ref.idpos.pos; // position of rep of end-kmer
          cmppos = -1; // position of next nucleotide after kmer to compare.
          contig = mapper.getContig(cr_end).ref.contig;
          repequal = (p_fw.first == p_fw.first.rep()); // is end-kmer rep?
          reversed = ((pos >= 0) != repequal); // is end-kmer is reverse direction of contig?
          
          if (pos >= 0) {
            if (repequal) {
              cmppos = pos - p_fw.second + 1 + k;
            } else {
              assert(pos == 0);
              cmppos = pos -1 + p_fw.second -1;
            }
          } else {
            if (repequal) {
              cmppos = p_fw.second - 2;
            } else {
              cmppos = (-pos + 1 - k) - p_fw.second + 1 + k;
            }
          }

          //maxi = i + p_fw.second;
          jumpi = 1 + i + contig->seq.jump(s, i+k, cmppos, reversed);
          i++;
          iter.raise(km, rep);
          //while (iter != iterend && bf.contains(rep) && i < maxi) {
          while (iter != iterend && i < jumpi) {
            i++;
            iter.raise(km, rep);
          }
          //assert(i == jumpi);

        } else {
          assert(!kmer_maps_to_contig(bf, km, mapper).isEmpty());
          // The kmer maps to a contig, how much can we jump through it?
          contig = mapper.getContig(cr).ref.contig;
          pos = cr.ref.idpos.pos;
          contiglen = (int32_t) contig->seq.size();
          
          cmppos = -1;
          repequal = km == rep;
          reversed = (pos >= 0) != repequal;

          if (pos >= 0) {
            assert(contiglen-pos-k >= 0);
            if (repequal) {
              cmppos = pos + k;
              //maxi = i + contiglen-pos-k + 1;
            } else {
              cmppos = pos -1;
              //maxi = i + pos +1; 
            }
          } else {
            assert(-pos >= k-1);
            if (repequal) {
              cmppos = -pos -k;
              //maxi = i + 2 - (pos + k);
            } else {
              cmppos = -pos+1;
              //maxi = i + contiglen - cmppos + 1;
            }
          }
          
          jumpi = 1 + i + contig->seq.jump(s, i + k, cmppos, reversed);
          i++;
          iter.raise(km, rep);
          //while (iter != iterend && bf.contains(rep) && i < maxi) {
          while (iter != iterend && i < jumpi) {
            i++;
            iter.raise(km, rep);
          }
          //assert(i == jumpi);
        }
      }     
    }
  }
  cerr << "Number of reads " << n_read  << ", kmers stored " << mapper.size()<< endl;
  
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


static const char alpha[4] = {'A','C','G','T'};

ContigRef kmer_maps_to_contig(BloomFilter &bf, Kmer km, KmerMapper &mapper) {
  ContigRef cr = mapper.find(km);
  if (!cr.isEmpty()) {
    return cr;
  }
  size_t dist = 1;
  Kmer fw, bw, end = km;
  while (dist <= Kmer::k) {
    size_t fw_count = 0;
    int j = -1;
    for (int i = 0; i < 4; i++) {
      Kmer fw_rep = end.forwardBase(alpha[i]).rep();
      if (bf.contains(fw_rep)) {
        j = i;
        fw_count++;
        if (fw_count > 1) {
          break;
        }
      }
    }

    if (fw_count != 1) {
      break;
    }

    fw = end.forwardBase(alpha[j]);

    size_t bw_count = 0;
    for (int i = 0; i < 4; i++) {
      Kmer bw_rep = fw.backwardBase(alpha[i]).rep();
      if (bf.contains(bw_rep)) {
        bw_count++;
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
    dist++;
    cr = mapper.find(end);
    if (!cr.isEmpty()) {
      return cr; 
    }
  }
  return cr;
}

// use:  (end,dist) = find_contig_forward(bf,km,s);
// pre:  
// post: km is contained in a contig c with respect to the
//       bloom filter graph bf and end is the forward endpoint (wrt km direction)
//       and c contains dist kmers until the end (including km)
//       if s is not NULL the sequence of the contig is stored in s
pair<Kmer, size_t> find_contig_forward(BloomFilter &bf, Kmer km, string* s) {
  //TODO: add string return option
  assert(bf.contains(km.rep()));
  vector<char> v;
  if (s != NULL) {
    char t[Kmer::MAX_K+1];
    km.toString(t);
    for (int i = 0; i < Kmer::k; i++) {
      v.push_back(t[i]);
    }
  }
  Kmer end = km;
  
  size_t dist = 1;

  Kmer fw,bw;

  while (true) {
    assert(bf.contains(end.rep()));
    size_t fw_count = 0;
    int j = -1;
    for (int i = 0; i < 4; i++) {
      Kmer fw_rep = end.forwardBase(alpha[i]).rep();
      if (bf.contains(fw_rep)) {
        j = i;
        fw_count++;
        if (fw_count > 1) {
          break;
        }
      }
    }

    if (fw_count != 1) {
      break;
    }
    
    fw = end.forwardBase(alpha[j]);
    assert(0 <= j && j < 4);
    assert(bf.contains(fw.rep()));

    size_t bw_count = 0;
    for (int i = 0; i < 4; i++) {
      Kmer bw_rep = fw.backwardBase(alpha[i]).rep();
      if (bf.contains(bw_rep)) {
        bw_count++;
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
    dist++;
    if (s != NULL) {
      v.push_back(alpha[j]);
    }
  }

  if (s != NULL) {
    s->clear();
    s->reserve(v.size());
    s->insert(s->begin(),v.begin(),v.end());
  }
  return make_pair(end,dist);
}
