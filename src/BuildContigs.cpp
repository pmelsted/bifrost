#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <vector>

#include <thread>
#include <atomic>

#include "Common.hpp"
#include "CompressedSequence.hpp"
#include "Contig.hpp"
#include "BlockedBloomFilter.hpp"
#include "ContigMethods.hpp"
#include "KmerIterator.hpp"
#include "fastq.hpp"
#include "ContigMapper.hpp"
#include "KmerHashTable.h"

#include "minHashIterator.hpp"
#include "RepHash.hpp"

struct BuildUnitigs_ProgramOptions {
  bool verbose;
  size_t threads, k, g;
  string freads, output, graphfilename;
  size_t read_chunksize;
  size_t unitig_size; // not configurable
  vector<string> files;
  bool clipTips;
  bool deleteIsolated;
  BuildUnitigs_ProgramOptions() : verbose(false), threads(1), k(0), g(23), \
    read_chunksize(10000), unitig_size(1000000), clipTips(false), deleteIsolated(false) {}
};

// use:  BuildUnitigs_PrintUsage();
// pre:
// post: Information about the correct parameters to build unitigs has been printed to cerr
void BuildUnitigs_PrintUsage() {
  cout << endl << "BFGraph " << BFG_VERSION << endl;
  cout << "Create a compacted de Bruijn graph from filtered FASTA/FASTQ files and save it to a GFA file" << endl << endl;
  cout << "Usage: BFGraph unitigs [arguments] ... FASTQ_files";
  cout << endl << endl << "Required arguments:" << endl <<
       "  -t, --threads=INT           Number of threads (default is 1)" << endl <<
       "  -k, --kmer-size=INT         Size of k-mers, same as for filtering reads (default is 31)" << endl <<
       "  -g, --min-size=INT          Size of minimizers, same as for filtering reads (default is 23)" << endl <<
       "  -f, --filtered=STRING       Filtered reads file" << endl <<
       "  -o, --output=STRING         Prefix for output GFA file" << endl <<
       "  -c, --chunk-size=INT        Read chunksize to split betweeen threads (default is 10000)" << endl <<
       endl << "Optional arguments:" << endl <<
       "  -n, --clip-tips             Clip tips shorter than k k-mers in length" << endl <<
       "  -d, --rm-isolated           Delete isolated unitigs shorter than k k-mers in length" << endl <<
       "  -v, --verbose               Print lots of messages during run" << endl;
}


// use:  BuildUnitigs_ParseOptions(argc, argv, opt);
// pre:  argc is the parameter count, argv is a list of valid parameters
//       like BuildUnitigs_PrintUsage describes and opt is ready to contain the parsed parameters
// post: All the parameters from argv have been parsed into opt
void BuildUnitigs_ParseOptions(int argc, char **argv, BuildUnitigs_ProgramOptions& opt) {
  const char *opt_string = "vt:k:g:f:o:c:nd";
  static struct option long_options[] = {
    {"verbose", no_argument, 0, 'v'},
    {"threads", required_argument, 0, 't'},
    {"kmer-size", required_argument, 0, 'k'},
    {"min-size", required_argument, 0, 'g'},
    {"filtered", required_argument, 0, 'f'},
    {"output", required_argument, 0, 'o'},
    {"chunk-size", required_argument, 0, 'c'},
    {"clip-tips", no_argument, 0, 'n'},
    {"rm-isolated", no_argument, 0, 'd'},
    {0, 0, 0, 0}
  };

  int option_index = 0, c;
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
    switch (c) {
    case 0: break;
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
    case 'f':
      opt.freads = optarg;
      break;
    case 'o':
      opt.output = optarg;
      break;
    case 'c':
      opt.read_chunksize = atoi(optarg);
      break;
    case 'n':
      opt.clipTips = true;
      break;
    case 'd':
      opt.deleteIsolated = true;
      break;
    default: break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  while (optind < argc) {
    opt.files.push_back(argv[optind++]);
  }
}


// use:  b = BuildUnitigs_CheckOptions(opt);
// pre:  opt contains parameters for building unitigs
// post: (b == true)  <==>  the parameters are valid
bool BuildUnitigs_CheckOptions(BuildUnitigs_ProgramOptions& opt) {
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

  if (opt.k == 0 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error: Invalid kmer-size: " << opt.k
         << ", need a number between 1 and " << (MAX_KMER_SIZE-1) << endl;
    ret = false;
  }

  if (opt.g <= 0 || opt.g >= MAX_KMER_SIZE) {
    cerr << "Error, invalid value for min-size: " << opt.g << endl;
    cerr << "Values must be between 1 and " << (MAX_KMER_SIZE-1) << endl;
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

  opt.graphfilename = opt.output + ".gfa";
  FILE *fp;
  if ((fp = fopen(opt.graphfilename.c_str(), "w")) == NULL) {
    cerr << "Error: Could not open file for writing: " << opt.graphfilename << endl;
    ret = false;
  } else {
    fclose(fp);
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


// use:  BuildUnitigs_PrintSummary(opt);
// pre:  opt has information about Kmer size, input file and output file
// post: Information about the Kmer size and the input and output files has been printed to cerr
void BuildUnitigs_PrintSummary(const BuildUnitigs_ProgramOptions& opt) {
  cerr << "Kmer size: " << opt.k << endl
       << "Chunksize: " << opt.read_chunksize << endl
       << "Reading file with filtered reads: " << opt.freads << endl
       << "fasta/fastq files: " << endl;
  vector<string>::const_iterator it;
  for (it = opt.files.begin(); it != opt.files.end(); ++it) {
    cerr << "  " << *it << endl;
  }
}

size_t cstrMatch(const char* a, const char* b) {

    char* _a = const_cast<char*>(a);
    char* _b = const_cast<char*>(b);

    while ((*_a != '\0') && (*_b != '\0') && (*_a == *_b)){ _a++; _b++; }

    return _a - a;
}

// use:  BuildUnitigs_Normal(opt);
// pre:  opt has information about Kmer size, input file and output file
// post: The unitigs have been written to the output file
void BuildUnitigs_Normal(const BuildUnitigs_ProgramOptions& opt) {
    /**
    *  outline of algorithm:
    *  open bloom filter file
    *    create unitig datastructures
    *    for each read
    *      for all kmers in read
    *        if kmer is in bf
    *          if it maps to unitig
    *            try to jump over as many kmers as possible
    *          else
    *            create new unitig from kmers in both directions \
    *            from this kmer while there is only one possible next kmer \
    *            with respect to the bloom filter
    *            when the unitig is ready,
    *            try to jump over as many kmers as possible
    */

    BlockedBloomFilter bf;

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

    }
    else {

        fclose(f);
        f = NULL;
    }

    UnitigMapper cmap;

    cmap.mapBloomFilter(&bf);

    KmerIterator iter, iterend;
    FastqFile FQ(opt.files);

    char name[8192];
    string s;
    size_t name_len, len;
    uint64_t n_read = 0;

    size_t read_chunksize = opt.read_chunksize;
    vector<string> readv;

    if (opt.verbose) cerr << "Starting real work ....." << endl << endl;

    const bool multi_threaded = (opt.threads != 1);

    tiny_vector<Kmer, 3>* fp_candidate = new tiny_vector<Kmer, 3>[cmap.bf->getNbBlocks()];
    KmerHashTable<bool> ignored_km_tips;

    auto worker_function = [&opt, &cmap, &fp_candidate](vector<string>::const_iterator a, vector<string>::const_iterator b,
                                                        vector<NewUnitig>* smallv,
                                                        vector<tuple<bool, uint64_t, Kmer>>* v_fp_cand,
                                                        vector<Kmer>* l_ignored_km_tip) {

        uint64_t it_min_h, last_it_min_h;
        std::pair<uint64_t*, uint64_t*> block_bf;

        Kmer km;
        RepHash rep;

        // for each input
        for (auto x = a; x != b; ++x) {

            const char* s_x = x->c_str();

            KmerHashIterator<RepHash> it_kmer_h(s_x, x->length(), opt.k), it_kmer_h_end;
            preAllocMinHashIterator<RepHash> it_min(s_x, x->length(), opt.k, opt.g, rep, true);

            for (int last_pos_km = -2; it_kmer_h != it_kmer_h_end; it_kmer_h++, it_min++) {

                std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                if (p_.second != last_pos_km + 1){ // If one or more k-mer were jumped because contained non-ACGT char.

                    km = Kmer(&s_x[p_.second]);

                    it_min += (last_pos_km == -2 ? p_.second : (p_.second - last_pos_km) - 1);
                    it_min_h = it_min.getHash();

                    block_bf = cmap.bf->getBlock(it_min_h);
                }
                else {

                    km.selfForwardBase(s_x[p_.second + opt.k - 1]);

                    it_min_h = it_min.getHash();
                    if (it_min_h != last_it_min_h) block_bf = cmap.bf->getBlock(it_min_h);
                }

                last_pos_km = p_.second;
                last_it_min_h = it_min_h;

                size_t r = cmap.bf->contains_block(p_.first, block_bf);

                if (r != 0){

                    UnitigMap cm = cmap.findUnitig(km, *x, p_.second, it_min);

                    if (cm.isEmpty) { // kmer did not map, push into queue for next unitig generation round

                        bool selfLoop = false;
                        bool isIsolated = false;

                        string newseq;

                        size_t pos_match = cmap.findUnitigSequence(km, newseq, selfLoop, isIsolated, *l_ignored_km_tip); //Build unitig from Bloom filter

                        if (isIsolated){ // According to the BF, k-mer is isolated in the graph and is a potential false positive

                            const uint64_t block = ((r == 1 ? block_bf.first : block_bf.second) - cmap.bf->getTable_ptr()) / NB_ELEM_BLOCK;

                            Kmer km_rep = km.rep();
                            const tiny_vector<Kmer, 3>& v = fp_candidate[block];
                            tuple<bool, uint64_t, Kmer> t_fp_cand = make_tuple(true, block, km_rep);

                            for (auto km_tmp : v){ // Go over global list of existing FP candidates

                                if (km_tmp == km_rep){ // If already stored as a FP candidate, it is a TP

                                    std::get<0>(t_fp_cand) = false; // K-mer must be removed from list of FP candidate
                                    break;
                                }
                            }

                            if (!std::get<0>(t_fp_cand)){ // If k-mer was a FP candidate now turned into a TP, insert into data structure

                                smallv->emplace_back(km, *x, p_.second, selfLoop ? std::string() : newseq);
                                v_fp_cand->emplace_back(t_fp_cand); // TP will be removed from list of FP later
                            }
                            else { // Need to make sure the FP candidate was not already inserted locally as a FP candidate (else it is a TP)

                                bool found_fp = false;

                                for (auto& t_fp_cand_tmp : *v_fp_cand){ // Go over local (to the thread) list of existing FP candidates

                                    if (std::get<2>(t_fp_cand_tmp) == km_rep){ // If present as FP candidate

                                        std::get<0>(t_fp_cand_tmp) = false; // Indicate FP must be removed from list of false positive candidate
                                        found_fp = true;
                                        break;
                                    }
                                }

                                if (found_fp) smallv->emplace_back(km, *x, p_.second, selfLoop ? std::string() : newseq);
                                else v_fp_cand->emplace_back(t_fp_cand);
                            }
                        }
                        else {

                            smallv->emplace_back(km, *x, p_.second, selfLoop ? std::string() : newseq);
                            it_kmer_h += cstrMatch(&s_x[p_.second + opt.k], &newseq.c_str()[pos_match + opt.k]);
                        }
                    }
                    else {

                        cmap.mapRead(cm);
                        it_kmer_h += cm.len - 1;
                    }
                }
            }
        }
    };

    vector<vector<NewUnitig>> parray(opt.threads);
    vector<vector<tuple<bool, uint64_t, Kmer>>> v_fp_cand_threads(opt.threads);
    vector<vector<Kmer>> v_ignored_km_tip_thread(opt.threads);

    int round = 0;
    bool done = false;

    while (!done) {

        readv.clear();

        size_t reads_now = 0;

        while (reads_now < read_chunksize) {

            if (FQ.read_next(name, &name_len, s, &len, NULL, NULL) >= 0){

                readv.emplace_back(s);
                ++reads_now;
                ++n_read;
            }
            else {

                done = true;
                break;
            }
        }

        ++round;

        if (read_chunksize > 1 && opt.verbose) cerr << "starting round " << round << endl;

        // run parallel code
        vector<thread> workers;
        auto rit = readv.begin();
        size_t batch_size = readv.size() / opt.threads;
        size_t leftover   = readv.size() % opt.threads;

        for (size_t i = 0; i < opt.threads; i++) {

            size_t jump = batch_size + ((i < leftover ) ? 1 : 0);
            auto rit_end(rit);

            advance(rit_end, jump);
            workers.push_back(thread(worker_function, rit, rit_end, &parray[i],
                                     &v_fp_cand_threads[i], &v_ignored_km_tip_thread[i]));
            rit = rit_end;
        }

        assert(rit == readv.end());

        for (auto& t : workers) t.join();

        for (auto &v : parray) { // for each thread

            for (auto &x : v) cmap.addUnitigSequence(x.km, x.read, x.pos, x.seq, v_ignored_km_tip_thread[0]); // add each unitig for this thread

            v.clear(); //clear the map
        }

        for (auto &v : v_fp_cand_threads) { // for each thread

            for (auto &x : v){ //for each FP cadndiate to add or delete

                Kmer km_fp = std::get<2>(x);

                tiny_vector<Kmer, 3>& v = fp_candidate[std::get<1>(x)];

                size_t i = 0;

                for (; i < v.size(); i++){ //Search vector for k-mer

                    if (v[i] == km_fp) break;
                }

                if (std::get<0>(x) == false){ // If FP is to delete

                    if (i < v.size()) v.remove(i); // If k-mer has not already been deleted

                    UnitigMap cm = cmap.find(km_fp); // Find FP to delete (was inserted into data structure as TP)

                    if (cm.isEmpty){
                        cerr << "Deleted false positive candidate was not inserted into data structure" << endl;
                        exit(1);
                    }

                    cmap.mapRead(cm); //Map the TP a second time;
                }
                else if (i >= v.size()){ //

                    if (!multi_threaded || cmap.find(km_fp).isEmpty) v.push_back(km_fp); //Add k-mer to vector if not already present
                }
                else { // K-mer already added as a FP candidate by other thread: it is now a TP

                    v.remove(i); //// K-mer deleted from list of FP candidates

                    const string str_km = km_fp.toString();

                    cmap.addUnitigSequence(km_fp, str_km, 0, str_km, v_ignored_km_tip_thread[0]);

                    UnitigMap cm = cmap.find(km_fp); // Find FP to delete (was inserted into data structure as TP)

                    if (cm.isEmpty){
                        cerr << "Deleted false positive candidate was not inserted into data structure" << endl;
                        exit(1);
                    }

                    cmap.mapRead(cm); //Map the TP a second time;
                }
            }

            v.clear(); //clear the map
        }

        for (auto &v : v_ignored_km_tip_thread) { // for each thread

            for (auto x : v) ignored_km_tips.insert(make_pair(x, false));

            v.clear();
        }

        if (read_chunksize > 1 && opt.verbose ) {

            cerr << " end of round" << endl;
            cerr << " processed " << cmap.unitigCount() << " unitigs" << endl;
        }
    }

    FQ.close();

    bf.clear();

    parray.clear();
    v_fp_cand_threads.clear();
    v_ignored_km_tip_thread.clear();

    delete[] fp_candidate;

    if (opt.verbose) {

        cerr << "Closed all fasta/fastq files" << endl;
        cerr << "Splitting unitigs" << endl;
    }

    size_t unitigsBefore = cmap.unitigCount();

    cerr << endl << "--- Splitting unitigs (1/2) ---" << endl;
    pair<size_t, size_t> unitigSplit = cmap.splitAllUnitigs();

    int unitigsAfter1 = cmap.unitigCount();

    cerr << endl << "--- Splitting unitigs (2/2) ---" << endl;
    cmap.check_fp_tips(ignored_km_tips);
    ignored_km_tips.clear();

    int unitigsAfter2 = cmap.unitigCount();

    //if (opt.verbose) {

        cerr << "Before split: " << unitigsBefore << " unitigs" << endl;
        cerr << "After split (1/2): " << unitigsAfter1 << " unitigs" <<  endl;
        cerr << "After split (2/2): " << unitigsAfter2 << " unitigs" <<  endl;
        cerr << "Unitigs split: " << unitigSplit.first << endl;
        cerr << "Unitigs deleted: " << unitigSplit.second << endl;
    //}

    cerr << endl << "--- Joining unitigs ---" << endl;
    size_t joined = cmap.joinAllUnitigs();

    int unitigsAfter3 = cmap.unitigCount();

    //if (opt.verbose) {
        cerr << "After join: " << unitigsAfter3 << " unitigs" << endl;
        cerr << "Joined " << joined << " unitigs" << endl;
    //}

    //cmap.checkIntegrity();

    if (opt.deleteIsolated || opt.clipTips){

        cerr << endl << "--- Removing isolated unitigs and/or clipping tips ---" << endl;

        vector<Kmer> v_joins;

        size_t removed = cmap.removeUnitigs(opt.deleteIsolated, opt.clipTips, v_joins);

        if (opt.clipTips) joined = cmap.joinAllUnitigs(&v_joins);
        else joined = 0;

        int unitigsAfter4 = cmap.unitigCount();

        cerr << "After: " << unitigsAfter4 << " unitigs" << endl;
        cerr << "Removed " << removed << " unitigs" << endl;
        cerr << "Joined " << joined << " unitigs" << endl;

        v_joins.clear();
    }

    cerr << endl << "--- Creating GFA ---" << endl;
    cmap.writeGFA(opt.graphfilename);
}



// use:  BuildUnitigs(argc, argv);
// pre:  argc is the number of arguments in argv and argv includes
//       arguments for "building the unitigs", including filenames
// post: If the number of arguments is correct and the arguments are valid
//       the "unitigs have been built" and written to a file
void BuildUnitigs(int argc, char **argv) {

    BuildUnitigs_ProgramOptions opt;

    BuildUnitigs_ParseOptions(argc,argv,opt);

    if (argc < 2) {

        BuildUnitigs_PrintUsage();
        exit(1);
    }

    if (!BuildUnitigs_CheckOptions(opt)) {

        BuildUnitigs_PrintUsage();
        exit(1);
    }

    // set static global k-value
    Kmer::set_k(opt.k);
    Minimizer::set_g(opt.g);

    if (opt.verbose) BuildUnitigs_PrintSummary(opt);

    BuildUnitigs_Normal(opt);

}
