#include "CompactedDBG.hpp"
#include "ColoredCDBG.hpp"

using namespace std;

// use:  PrintVersion();
// post: The version of the program has been printed to cout
void PrintVersion() {
    cout <<  BFG_VERSION << endl;
}

void PrintUsage() {

    cout << endl << "Bifrost " << BFG_VERSION << endl << endl;
    cout << "Highly parallel construction and indexing of colored and compacted de Bruijn graphs" << endl << endl;
    cout << "Usage: Bifrost [COMMAND] [GENERAL_PARAMETERS] [COMMAND_PARAMETERS]" << endl << endl;
    cout << "[COMMAND]:" << endl << endl <<
    "   build                   Build a compacted de Bruijn graph, with or without colors" << endl <<
    "   update                  Update a compacted (possible colored) de Bruijn graph with new sequences" << endl << endl;
    cout << "[GENERAL_PARAMETERS]:" << endl << endl;
    cout << "   > Mandatory with required argument:" << endl << endl <<
    "   -s, --input-seq-file     Input sequence file (FASTA/FASTQ possibly gzipped)" << endl <<
    "                            Multiple files can be provided as a list in a TXT file (one file per line)" << endl <<
    "                            K-mers with exactly 1 occurrence in the input sequence files will be discarded" << endl <<
    "   -r, --input-ref-file     Input reference file (FASTA/FASTQ possibly gzipped and GFA)" << endl <<
    "                            Multiple files can be provided as a list in a TXT file (one file per line)" << endl <<
    "                            All k-mers of the input reference files are used" << endl <<
    "   -o, --output-file        Prefix for output file(s)" << endl <<
    endl << "   > Optional with required argument:" << endl << endl <<
    "   -t, --threads            Number of threads (default is 1)" << endl <<
    endl << "   > Optional with no argument:" << endl << endl <<
    "   -i, --clip-tips          Clip tips shorter than k k-mers in length" << endl <<
    "   -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length" << endl <<
    "   -v, --verbose            Print information messages during execution" << endl << endl;
    cout << "[COMMAND_PARAMETERS]: build" << endl << endl;
    cout << "   > Optional with required argument:" << endl << endl <<
    "   -k, --kmer-length        Length of k-mers (default is 31)" << endl <<
    "   -m, --min-length         Length of minimizers (default is 23)" << endl <<
    "   -b, --bloom-bits         Number of Bloom filter bits per k-mer with 1+ occurrences in the input files (default is 14)" << endl <<
    "   -B, --bloom-bits2        Number of Bloom filter bits per k-mer with 2+ occurrences in the input files (default is 14)" << endl <<
    "   -l, --load-mbbf          Input Blocked Bloom Filter file, skips filtering step (default is no input)" << endl <<
    "   -w, --write-mbbf         Output Blocked Bloom Filter file (default is no output)" << endl <<
    "   -u, --chunk-size         Read chunk size per thread (default is 64)" << endl <<
    endl << "   > Optional with no argument:" << endl << endl <<
    "   -c, --colors             Color the compacted de Bruijn graph (default is no coloring)" << endl <<
    "   -y, --keep-mercy         Keep low coverage k-mers connecting tips" << endl <<
    "   -a, --fasta              Output file is in FASTA format (only sequences) instead of GFA" << endl << endl;
    cout << "[COMMAND_PARAMETERS]: update" << endl << endl;
    cout << "   > Mandatory with required argument:" << endl << endl <<
    "   -g, --input-graph-file   Input graph file to update (GFA format)" << endl <<
    endl << "   > Optional with required argument:" << endl << endl <<
    "   -f, --input-color-file   Input color file associated with the input graph file to update" << endl <<
    "   -k, --kmer-length        Length of k-mers (default is read from input graph file if built with Bifrost or 31)" << endl <<
    "   -m, --min-length         Length of minimizers (default is read from input graph file if built with Bifrost or 23)" << endl << endl;
}

void parse_ProgramOptions(int argc, char **argv, CCDBG_Build_opt& opt) {

    int option_index = 0, c;

    //const char* opt_string = "s:r:g:o:t:k:m:n:N:b:B:l:w:u:f:idvcya";
    const char* opt_string = "s:r:g:o:t:k:m:b:B:l:w:u:f:idvcya";

    static struct option long_options[] = {

        {"input-seq-file",     required_argument,  0, 's'},
        {"input-ref-file",     required_argument,  0, 'r'},
        {"input-graph-file",    required_argument,  0, 'g'},
        {"output-file",         required_argument,  0, 'o'},
        {"threads",             required_argument,  0, 't'},
        {"kmer-length",         required_argument,  0, 'k'},
        {"min-length",          required_argument,  0, 'm'},
        //{"num-kmers",           required_argument,  0, 'n'},
        //{"num-kmers2",          required_argument,  0, 'N'},
        {"bloom-bits",          required_argument,  0, 'b'},
        {"bloom-bits2",         required_argument,  0, 'B'},
        {"load-mbbf",           required_argument,  0, 'l'},
        {"write-mbbf",          required_argument,  0, 'w'},
        {"chunk-size",          required_argument,  0, 'u'},
        {"input-color-file",    required_argument,  0, 'f'},
        {"clip-tips",           no_argument,        0, 'i'},
        {"del-isolated",        no_argument,        0, 'd'},
        {"verbose",             no_argument,        0, 'v'},
        {"colors",              no_argument,        0, 'c'},
        {"keep-mercy",          no_argument,        0, 'y'},
        {"fasta",               no_argument,        0, 'a'},
        {0,                     0,                  0,  0 }
    };

    if (strcmp(argv[1], "build") == 0) opt.build = true;
    else if (strcmp(argv[1], "update") == 0) opt.update = true;

    if (opt.build || opt.update){

        while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

            switch (c) {

                case 's':
                    opt.filename_seq_in.push_back(optarg);
                    break;
                case 'r':
                    opt.filename_ref_in.push_back(optarg);
                    break;
                case 'g':
                    opt.filename_graph_in = optarg;
                    break;
                case 'f':
                    opt.filename_colors_in = optarg;
                    break;
                case 'o':
                    opt.prefixFilenameOut = optarg;
                    break;
                case 't':
                    opt.nb_threads = atoi(optarg);
                    break;
                case 'k':
                    opt.k = atoi(optarg);
                    break;
                case 'm':
                    opt.g = atoi(optarg);
                    break;
                /*case 'n':
                    opt.nb_unique_kmers = atoi(optarg);
                    break;
                case 'N':
                    opt.nb_non_unique_kmers = atoi(optarg);
                    break;*/
                case 'b':
                    opt.nb_bits_unique_kmers_bf = atoi(optarg);
                    break;
                case 'B':
                    opt.nb_bits_non_unique_kmers_bf = atoi(optarg);
                    break;
                case 'w':
                    opt.outFilenameBBF = optarg;
                    break;
                case 'l':
                    opt.inFilenameBBF = optarg;
                    break;
                case 'u':
                    opt.read_chunksize = atoi(optarg);
                    break;
                case 'i':
                    opt.clipTips = true;
                    break;
                case 'd':
                    opt.deleteIsolated = true;
                    break;
                case 'v':
                    opt.verbose = true;
                    break;
                case 'c':
                    opt.outputColors = true;
                    break;
                case 'y':
                    opt.useMercyKmers = true;
                    break;
                case 'a':
                    opt.outputGFA = false;
                    break;
                default: break;
            }
        }
    }
}

bool check_ProgramOptions(CCDBG_Build_opt& opt) {

    bool ret = true;

    size_t max_threads = std::thread::hardware_concurrency();

    auto check_files = [&](vector<string>& v_files) {

        vector<string> files_tmp;

        char* buffer = new char[4096]();

        for (const auto& file : v_files) {

            if (!check_file_exists(file)) {

                cerr << "Error: File " << file << " not found." << endl;
                ret = false;
            }
            else {

                const string s_ext = file.substr(file.find_last_of(".") + 1);

                if ((s_ext == "txt")){

                    FILE* fp = fopen(file.c_str(), "r");

                    if (fp != NULL){

                        fclose(fp);

                        ifstream ifs_file_txt(file);
                        istream i_file_txt(ifs_file_txt.rdbuf());

                        while (i_file_txt.getline(buffer, 4096)){

                            fp = fopen(buffer, "r");

                            if (fp == NULL) {

                                cerr << "Error: Could not open file " << buffer << " for reading." << endl;
                                ret = false;
                            }
                            else {

                                fclose(fp);
                                files_tmp.push_back(string(buffer));
                            }
                        }

                        ifs_file_txt.close();
                    }
                    else {

                        cerr << "Error: Could not open file " << file << " for reading." << endl;
                        ret = false;
                    }
                }
                else files_tmp.push_back(file);
            }
        }

        v_files = move(files_tmp);

        delete[] buffer;
    };

    // Check general parameters

    if (!opt.build && !opt.update){

        cerr << "Error: No command selected (can be 'build' or 'update')." << endl;
        ret = false;
    }

    if (opt.nb_threads <= 0){

        cerr << "Error: Number of threads cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "Error: Number of threads cannot be greater than or equal to " << max_threads << "." << endl;
        ret = false;
    }

    if (opt.k <= 0){

        cerr << "Error: Length k of k-mers cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.k >= MAX_KMER_SIZE){

        cerr << "Error: Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << "." << endl;
        ret = false;
    }

    if (opt.g <= 0){

        cerr << "Error: Length m of minimizers cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.g > opt.k - 2){

        cerr << "Error: Length m of minimizers cannot exceed k - 2 (" << (opt.k - 2) << ")." << endl;
        ret = false;
    }

    const string out = opt.prefixFilenameOut + (opt.outputGFA ? ".gfa" : ".fasta");

    FILE* fp = fopen(out.c_str(), "w");

    if (fp == NULL) {

        cerr << "Error: Could not open file for writing output graph in GFA format: " << out << "." << endl;
        ret = false;
    }
    else {

        fclose(fp);
        if (remove(out.c_str()) != 0) cerr << "Error: Could not remove temporary file " << out << "." << endl;
    }

    if ((opt.filename_seq_in.size() + opt.filename_ref_in.size()) == 0) {

        cerr << "Error: Missing input files." << endl;
        ret = false;
    }
    else {

        check_files(opt.filename_seq_in);
        check_files(opt.filename_ref_in);
    }

    if (opt.build){ // Check param. command build

        if (opt.read_chunksize <= 0) {

            cerr << "Error: Chunk size of reads to share among threads cannot be less than or equal to 0." << endl;
            ret = false;
        }

        if (opt.outFilenameBBF.length() != 0){

            FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

            if (fp == NULL) {

                cerr << "Error: Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
                ret = false;
            }
            else {

                fclose(fp);

                if (remove(opt.outFilenameBBF.c_str()) != 0){

                    cerr << "Error: Could not remove temporary file " << opt.outFilenameBBF << "." << endl;
                }
            }
        }

        if (opt.inFilenameBBF.length() != 0){

            if (check_file_exists(opt.inFilenameBBF)){

                FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "Error: Could not read input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                    ret = false;
                }
                else fclose(fp);
            }
            else {

                cerr << "Error: Input Blocked Bloom filter " << opt.inFilenameBBF << " file does not exist." << endl;
                ret = false;
            }
        }
    }

    if (opt.update){

        if (opt.filename_graph_in.length() == 0){

            cerr << "Error: No graph file to update was provided in input." << endl;
            ret = false;
        }
        else if (!check_file_exists(opt.filename_graph_in)){

            cerr << "Error: The graph file to update does not exist." << endl;
            ret = false;
        }
        else {

            FILE* fp = fopen(opt.filename_graph_in.c_str(), "r");

            if (fp == NULL) {

                cerr << "Error: Could not read input graph file " << opt.filename_graph_in << "." << endl;
                ret = false;
            }
            else fclose(fp);
        }

        if (opt.filename_colors_in.length() != 0){

            if (!check_file_exists(opt.filename_colors_in)){

                cerr << "Error: The input color file does not exist." << endl;
                ret = false;
            }
            else {

                FILE* fp = fopen(opt.filename_colors_in.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "Error: Could not read input color file " << opt.filename_colors_in << "." << endl;
                    ret = false;
                }
                else fclose(fp);
            }
        }
    }

    return ret;
}

int main(int argc, char **argv){

    if (argc < 2) PrintUsage();
    else {

        CCDBG_Build_opt opt;

        opt.outputColors = false; // We dont know yet if we want colors or not

        parse_ProgramOptions(argc, argv, opt); // Parse input parameters

        if (!check_ProgramOptions(opt)) return 0; // Check if input parameters are valid
        else if (opt.build){ // We want to build the graph

            if (opt.outputColors){

                ColoredCDBG<> cdbg(opt.k, opt.g);

                cdbg.buildGraph(opt);
                cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
                cdbg.buildColors(opt);
                cdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);
            }
            else {

                CompactedDBG<> cdbg(opt.k, opt.g);

                cdbg.build(opt);
                cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
                cdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.outputGFA, opt.verbose);
            }
        }
        else if (opt.update){

            if (opt.filename_colors_in.length() != 0){ // If colors in or out

                ColoredCDBG<> cdbg1(opt.k, opt.g);
                ColoredCDBG<> cdbg2(opt.k, opt.g);

                cdbg1.read(opt.filename_graph_in, opt.filename_colors_in, opt.nb_threads, opt.verbose);
                cdbg2.buildGraph(opt);
                cdbg2.buildColors(opt);

                const size_t cdbg1_len = cdbg1.length();
                const size_t cdbg2_len = cdbg2.length();

                ColoredCDBG<>& cdbg_a = (cdbg1_len > cdbg2_len) ? cdbg1 : cdbg2;
                ColoredCDBG<>& cdbg_b = (cdbg1_len > cdbg2_len) ? cdbg2 : cdbg1;

                cdbg_a.merge(move(cdbg_b), opt.nb_threads, opt.verbose);

                cdbg_a.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
                cdbg_a.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);
            }
            else {

                CompactedDBG<> cdbg1(opt.k, opt.g);
                CompactedDBG<> cdbg2(opt.k, opt.g);

                cdbg1.read(opt.filename_graph_in, opt.verbose);
                cdbg2.build(opt);

                const size_t cdbg1_len = cdbg1.length();
                const size_t cdbg2_len = cdbg2.length();

                CompactedDBG<>& cdbg_a = (cdbg1_len > cdbg2_len) ? cdbg1 : cdbg2;
                CompactedDBG<>& cdbg_b = (cdbg1_len > cdbg2_len) ? cdbg2 : cdbg1;

                cdbg_a.merge(cdbg_b, opt.nb_threads, opt.verbose);
                cdbg_b.clear();

                cdbg_a.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
                cdbg_a.write(opt.prefixFilenameOut, opt.nb_threads, opt.outputGFA, opt.verbose);
            }
        }
    }
}
