#include "CompactedDBG.hpp"
#include "ColoredCDBG.hpp"

using namespace std;

// use:  PrintVersion();
// post: The version of the program has been printed to cout
void PrintVersion() {
    cout <<  BFG_VERSION << endl;
}

// use:  PrintCite();
// post: Information of how to cite this software has been printed to cerr
void PrintCite() {
    cout << "The paper describing this software has not been published." << endl;
}

void PrintUsage() {

    cout << endl << "Bifrost " << BFG_VERSION << endl << endl;
    cout << "Highly Parallel and Memory Efficient Colored and Compacted de Bruijn Graph Construction" << endl << endl;
    cout << "Usage: Bifrost [Parameters] -o <output_prefix> -f <file_1> ..." << endl << endl;
    cout << "Mandatory parameters with required argument:" << endl << endl <<
    "  -f, --input-files        Input sequence files (FASTA or FASTQ, possibly gziped) and/or graph files (GFA)" << endl <<
    "  -o, --output-file        Prefix for output file (GFA output by default)" << endl << endl <<
    "Optional parameters with required argument:" << endl << endl <<
    "  -t, --threads            Number of threads (default is 1)" << endl <<
    "  -k, --kmer-length        Length of k-mers (default is 31)" << endl <<
    //"  -c, --load-colors        Input color sets file (default is no color)" << endl <<
    "  -g, --min-length         Length of minimizers (default is 23)" << endl <<
    "  -n, --num-kmers          Estimated number of different k-mers in input files (default: KmerStream estimation)" << endl <<
    "  -N, --num-kmers2         Estimated number of different k-mers occurring twice or more in the input files (default: KmerStream estimation)" << endl <<
    "  -b, --bloom-bits         Number of Bloom filter bits per k-mer occurring at least once in the input files (default is 14)" << endl <<
    "  -B, --bloom-bits2        Number of Bloom filter bits per k-mer occurring at least twice in the input files (default is 14)" << endl <<
    "  -l, --load-mbbf          Filename for input Blocked Bloom Filter, skips filtering step (default is no input)" << endl <<
    "  -w, --write-mbbf         Filename for output Blocked Bloom Filter (default is no output)" << endl <<
    "  -s, --chunk-size         Read chunk size per thread (default is 64)" << endl <<
    endl << "Optional parameters with no argument:" << endl << endl <<
    "  -p, --produce-colors     Produce a colored and compacted de Bruijn graph (default is without colors)" << endl <<
    "  -r, --reference          Reference mode, no filtering" << endl <<
    "  -i, --clip-tips          Clip tips shorter than k k-mers in length" << endl <<
    "  -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length" << endl <<
    "  -m, --keep-mercy         Keep low coverage k-mers connecting tips" << endl <<
    "  -a, --fasta              Output file is in FASTA format (only sequences) instead of GFA" << endl <<
    "  -v, --verbose            Print information messages during construction" << endl <<
    endl;
}

void parse_ProgramOptions(int argc, char **argv, CCDBG_Build_opt& opt) {

    const char* opt_string = "f:o:c:n:N:t:k:g:b:B:l:w:s:pridmav";

    static struct option long_options[] = {

        {"input-files",     required_argument,  0, 'f'},
        {"output-file",     required_argument,  0, 'o'},
        {"load-colors",     required_argument,  0, 'c'},
        {"num-kmers",       required_argument,  0, 'n'},
        {"num-kmers2",      required_argument,  0, 'N'},
        {"threads",         required_argument,  0, 't'},
        {"kmer-length",     required_argument,  0, 'k'},
        {"min-length",      required_argument,  0, 'g'},
        {"bloom-bits",      required_argument,  0, 'b'},
        {"bloom-bits2",     required_argument,  0, 'B'},
        {"load-mbbf",       required_argument,  0, 'l'},
        {"write-mbbf",      required_argument,  0, 'w'},
        {"chunk-size",      required_argument,  0, 's'},
        {"produce-colors",  no_argument,        0, 'p'},
        {"reference",       no_argument,        0, 'r'},
        {"clip-tips",       no_argument,        0, 'i'},
        {"del-isolated",    no_argument,        0, 'd'},
        {"keep-mercy",      no_argument,        0, 'm'},
        {"fasta",           no_argument,        0, 'a'},
        {"verbose",         no_argument,        0, 'v'},
        {0,                 0,                  0,  0 }
    };

    int option_index = 0, c;

    while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

        switch (c) {

            case 'p':
                opt.outputColors = true;
                break;
            case 'r':
                opt.reference_mode = true;
                break;
            case 'v':
                opt.verbose = true;
                break;
            case 'i':
                opt.clipTips = true;
                break;
            case 'd':
                opt.deleteIsolated = true;
                break;
            case 'm':
                opt.useMercyKmers = true;
                break;
            case 'a':
                opt.outputGFA = false;
                break;
            case 't':
                opt.nb_threads = atoi(optarg);
                break;
            case 'k':
                opt.k = atoi(optarg);
                break;
            case 'g':
                opt.g = atoi(optarg);
                break;
            case 's':
                opt.read_chunksize = atoi(optarg);
                break;
            case 'n':
                opt.nb_unique_kmers = atoi(optarg);
                break;
            case 'N':
                opt.nb_non_unique_kmers = atoi(optarg);
                break;
            case 'b':
                opt.nb_bits_unique_kmers_bf = atoi(optarg);
                break;
            case 'B':
                opt.nb_bits_non_unique_kmers_bf = atoi(optarg);
                break;
            case 'o':
                opt.prefixFilenameOut = optarg;
                break;
            case 'w':
                opt.outFilenameBBF = optarg;
                break;
            case 'l':
                opt.inFilenameBBF = optarg;
                break;
            case 'f': {

                for (optind--; (optind < argc) && (*argv[optind] != '-'); ++optind){

                      opt.filename_seq_in.push_back(argv[optind]);
                }

                break;
            }
            case 'c': {

                for (optind--; (optind < argc) && (*argv[optind] != '-'); ++optind){

                      opt.filename_colors_in.push_back(argv[optind]);
                }

                break;
            }
            default: break;
        }
    }
}

bool check_ProgramOptions(CCDBG_Build_opt& opt) {

    bool ret = true;

    size_t max_threads = std::thread::hardware_concurrency();

    if (opt.nb_threads <= 0){

        cerr << "Error: Number of threads cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "Error: Number of threads cannot be greater than or equal to " << max_threads << endl;
        ret = false;
    }

    if (opt.read_chunksize <= 0) {

        cerr << "Error: Chunk size of reads to share among threads cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (opt.k <= 0){

        cerr << "Error: Length k of k-mers cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (opt.k >= MAX_KMER_SIZE){

        cerr << "Error: Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        ret = false;
    }

    if (opt.g <= 0){

        cerr << "Error: Length g of minimizers cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (opt.g >= MAX_KMER_SIZE){

        cerr << "Error: Length g of minimizers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        ret = false;
    }

    if (!opt.reference_mode && (opt.nb_non_unique_kmers > opt.nb_unique_kmers)){

        cerr << "Error: The estimated number of non unique k-mers ";
        cerr << "cannot be greater than the estimated number of unique k-mers" << endl;
        ret = false;
    }

    if (opt.outFilenameBBF.length() != 0){

        FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

        if (fp == NULL) {

            cerr << "Error: Could not open file for writing output Blocked Bloom filter: " << opt.outFilenameBBF << endl;
            ret = false;
        }
        else fclose(fp);
    }

    if (opt.inFilenameBBF.length() != 0){

        FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

        if (fp == NULL) {

            cerr << "Error: Could not read file input Blocked Bloom filter: " << opt.inFilenameBBF << endl;
            ret = false;
        }
        else fclose(fp);
    }

    const string out = opt.prefixFilenameOut + (opt.outputGFA ? ".gfa" : ".fasta");

    FILE* fp = fopen(out.c_str(), "w");

    if (fp == NULL) {

        cerr << "Error: Could not open file for writing output graph in GFA format: " << out << endl;
        ret = false;
    }
    else {

        fclose(fp);
        if (remove(out.c_str()) != 0) cerr << "Error: Could not remove temporary file " << out << endl;
    }

    if (opt.filename_seq_in.size() == 0) {

        cerr << "Error: Missing input files" << endl;
        ret = false;
    }
    else {

        struct stat stFileInfo;
        vector<string> filename_seq_in_tmp;

        char* buffer = new char[4096]();

        for (vector<string>::const_iterator it = opt.filename_seq_in.begin(); it != opt.filename_seq_in.end(); ++it) {

            const int iStat = stat(it->c_str(), &stFileInfo);

            if (iStat != 0) {

                cerr << "Error: File not found, " << *it << endl;
                ret = false;
            }

            const string s_ext = it->substr(it->find_last_of(".") + 1);

            if ((s_ext == "txt")){

                FILE* fp = fopen(it->c_str(), "r");

                if (fp != NULL){

                    fclose(fp);

                    ifstream ifs_file_txt(*it);
                    istream i_file_txt(ifs_file_txt.rdbuf());

                    while (i_file_txt.getline(buffer, 4096)){

                        fp = fopen(buffer, "r");

                        if (fp == NULL) {

                            cerr << "Error: Could not open file " << buffer << " for reading" << endl;
                            ret = false;
                        }
                        else {

                            fclose(fp);
                            filename_seq_in_tmp.push_back(string(buffer));
                        }
                    }

                    ifs_file_txt.close();
                }
                else {

                    cerr << "Error: Could not open file " << *it << " for reading" << endl;
                    ret = false;
                }
            }
            else filename_seq_in_tmp.push_back(*it);
        }

        opt.filename_seq_in = move(filename_seq_in_tmp);

        delete[] buffer;
    }

    if (opt.filename_colors_in.size() != 0) {

        struct stat stFileInfo;
        vector<string>::const_iterator it;
        int intStat;

        for (it = opt.filename_colors_in.begin(); it != opt.filename_colors_in.end(); ++it) {

            intStat = stat(it->c_str(), &stFileInfo);

            if (intStat != 0) {
                cerr << "Error: File not found, " << *it << endl;
                ret = false;
            }
        }
    }

    return ret;
}

int main(int argc, char **argv){

    if (argc < 2) PrintUsage();
    else {

        /*ColoredCDBG<> ccdbg;

        ccdbg.read("/home/wiz/Documents/my_git/bfgraph/build/test", 4, true);*/

        CCDBG_Build_opt opt;

        opt.reference_mode = false; // We dont know yet if we want colors or not
        opt.outputColors = false; // We dont know yet if we want colors or not

        parse_ProgramOptions(argc, argv, opt); // Parse input parameters

        if (check_ProgramOptions(opt)){ // Check if input parameters are valid

            if (opt.outputColors || (opt.filename_colors_in.size() != 0)){ // If colors in or out

                opt.reference_mode = true; //If the user wants a colored cdBG, it is reference mode
                opt.outputColors = true; //If the user wants a colored cdBG, it is reference mode

                ColoredCDBG<> cdbg(opt.k, opt.g);

                cdbg.buildGraph(opt);
                cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
                cdbg.buildColors(opt);
                cdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);
            }
            else {

                CDBG_Build_opt opt_ = opt.getCDBG_Build_opt();

                CompactedDBG<> cdbg(opt_.k, opt_.g);

                cdbg.build(opt_);
                cdbg.simplify(opt_.deleteIsolated, opt_.clipTips, opt_.verbose);
                cdbg.write(opt_.prefixFilenameOut, opt_.nb_threads, opt_.outputGFA, opt_.verbose);
            }
        }
    }
}
