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

#include <jemalloc/jemalloc.h>

#include "CompactedDBG.hpp"

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
    cout << "Highly Parallel and Memory Efficient Compacted de Bruijn Graph Construction" << endl << endl;
    cout << "Usage: Bifrost [Parameters] FAST(A|Q)_file_1 ..." << endl << endl;
    cout << "Parameters with required argument:" << endl << endl <<
    "  -n, --num-kmers          [MANDATORY] Estimated number (upper bound) of different k-mers in the FASTA/FASTQ files" << endl <<
    "  -N, --num-kmer2          [MANDATORY] Estimated number (upper bound) of different k-mers occurring twice or more in the FASTA/FASTQ files" << endl <<
    "  -o, --output             [MANDATORY] Prefix for output GFA file" << endl <<
    "  -t, --threads            Number of threads (default is 1)" << endl <<
    "  -k, --kmer-length        Length of k-mers (default is 31)" << endl <<
    "  -g, --min-length         Length of minimizers (default is 23)" << endl <<
    "  -b, --bloom-bits         Number of Bloom filter bits per k-mer occurring at least once in the FASTA/FASTQ files (default is 14)" << endl <<
    "  -B, --bloom-bits2        Number of Bloom filter bits per k-mer occurring at least twice in the FASTA/FASTQ files (default is 14)" << endl <<
    "  -l, --load               Filename for input Blocked Bloom Filter, skips filtering step (default is no input)" << endl <<
    "  -f, --output2            Filename for output Blocked Bloom Filter (default is no output)" << endl <<
    "  -s, --chunk-size         Read chunksize to split between threads (default is 10000)" << endl <<
    endl << "Parameters with no argument:" << endl << endl <<
    "      --ref                Reference mode, no filtering" << endl <<
    "  -c, --clip-tips          Clip tips shorter than k k-mers in length" << endl <<
    "  -r, --rm-isolated        Delete isolated contigs shorter than k k-mers in length" << endl <<
    "  -v, --verbose            Print information messages during construction" << endl <<
    endl;
}

void parse_ProgramOptions(int argc, char **argv, CDBG_Build_opt& opt) {

    const char* opt_string = "n:N:o:t:k:g:b:B:l:f:s:crv";

    static struct option long_options[] = {

        {"num-kmers",       required_argument,  0, 'n'},
        {"num-kmers2",      required_argument,  0, 'N'},
        {"output",          required_argument,  0, 'o'},
        {"threads",         required_argument,  0, 't'},
        {"kmer-length",     required_argument,  0, 'k'},
        {"min-length",      required_argument,  0, 'g'},
        {"bloom-bits",      required_argument,  0, 'b'},
        {"bloom-bits2",     required_argument,  0, 'B'},
        {"load",            required_argument,  0, 'l'},
        {"output2",         required_argument,  0, 'f'},
        {"chunk-size",      required_argument,  0, 's'},
        {"ref",             no_argument,        0,  0 },
        {"clip-tips",       no_argument,        0, 'c'},
        {"rm-isolated",     no_argument,        0, 'r'},
        {"verbose",         no_argument,        0, 'v'},
        {0,                 0,                  0,  0 }
    };

    int option_index = 0, c;

    while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

        switch (c) {

            case 0:
                if (strcmp(long_options[option_index].name, "ref") == 0) opt.reference_mode = true;
                break;
            case 'v':
                opt.verbose = true;
                break;
            case 'c':
                opt.clipTips = true;
                break;
            case 'r':
                opt.deleteIsolated = true;
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
                opt.prefixFilenameGFA = optarg;
                break;
            case 'f':
                opt.outFilenameBBF = optarg;
                break;
            case 'l':
                opt.inFilenameBBF = optarg;
                break;
            default: break;
        }
    }

    // all other arguments are fast[a/q] files to be read
    while (optind < argc) opt.fastx_filename_in.push_back(argv[optind++]);
}

bool check_ProgramOptions(CDBG_Build_opt& opt) {

    bool ret = true;

    size_t max_threads = std::thread::hardware_concurrency();

    if (opt.nb_threads <= 0){

        cerr << "Error: Number of threads cannot be less than or equal to 0" << endl;
        ret = false;
    }
    else if (opt.nb_threads > max_threads){

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

    if (opt.nb_unique_kmers <= 0){

        cerr << "Error: Number of Bloom filter bits per unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (!opt.reference_mode && (opt.nb_non_unique_kmers <= 0)){

        cerr << "Error: Number of Bloom filter bits per non unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (!opt.reference_mode && (opt.nb_non_unique_kmers > opt.nb_unique_kmers)){

        cerr << "Error: The estimated number of non unique k-mers ";
        cerr << "cannot be greater than the estimated number of unique k-mers" << endl;
        ret = false;
    }

    if (opt.nb_bits_unique_kmers_bf <= 0){

        cerr << "Error: Number of Bloom filter bits per unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (!opt.reference_mode && (opt.nb_bits_non_unique_kmers_bf <= 0)){

        cerr << "Error: Number of Bloom filter bits per non unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (opt.reference_mode) {

        opt.nb_bits_non_unique_kmers_bf = 0;
        opt.nb_non_unique_kmers = 0;
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

    opt.filenameGFA = opt.prefixFilenameGFA + ".gfa";

    FILE* fp = fopen(opt.filenameGFA.c_str(), "w");

    if (fp == NULL) {

        cerr << "Error: Could not open file for writing output graph in GFA format: " << opt.filenameGFA << endl;
        ret = false;
    }
    else fclose(fp);

    if (opt.fastx_filename_in.size() == 0) {

        cerr << "Error: Missing FASTA/FASTQ input files" << endl;
        ret = false;
    }
    else {

        struct stat stFileInfo;
        vector<string>::const_iterator it;
        int intStat;

        for(it = opt.fastx_filename_in.begin(); it != opt.fastx_filename_in.end(); ++it) {

            intStat = stat(it->c_str(), &stFileInfo);

            if (intStat != 0) {
                cerr << "Error: File not found, " << *it << endl;
                ret = false;
            }
        }
    }

    return ret;
}

class myInt : public CDBG_Data_t<myInt> {

    public:

        myInt(int int_init = 1) : Int(int_init) {}

        void join(const myInt& data, CompactedDBG<myInt>& cdbg){

            Int += data.Int;

            cout << "join: " << Int << " += "  << data.Int << endl;
        }

        void split(const size_t pos, const size_t len, myInt& new_data, CompactedDBG<myInt>& cdbg) const {

            new_data.Int = 0;
        }

        int Int;
};

class myColors : public CDBG_Data_t<myColors> {

    public:

        void join(const myColors& data, CompactedDBG<myColors>& cdbg){

            for (auto color : data.colors) add(color);
        }

        void split(const size_t pos, const size_t len, myColors& new_data, CompactedDBG<myColors>& cdbg) const {

        }

        void add(int color){

            if ((colors.size() == 0) || (color > colors[colors.size() - 1])) colors.push_back(color);
            else {

                for (vector<int>::iterator it = colors.begin(); it < colors.end(); it++){

                    if (*it == color) break;
                    else if (*it > color){

                        colors.insert(it, color);
                        break;
                    }
                }
            }
        }

        void erase(int color){

            if ((colors.size() > 0) && (color <= colors[colors.size() - 1])){

                for (vector<int>::iterator it = colors.begin(); it < colors.end(); it++){

                    if (*it == color){

                        colors.erase(it);
                        break;
                    }
                    else if (*it > color) break;
                }
            }
        }

        void print() const {

            for (auto color : colors){

                cout << "Color ID: " << color << endl;
            }
        }

    private:

        vector<int> colors; //color_id, position, len
};

int main(int argc, char **argv){

    if (argc < 2) PrintUsage();
    else {

        CDBG_Build_opt opt;

        parse_ProgramOptions(argc, argv, opt);

        if (check_ProgramOptions(opt)){ //Program options are valid

            CompactedDBG<> cdbg(opt.k, opt.g);

            cdbg.build(opt);

            cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);

            cdbg.write(opt.prefixFilenameGFA, opt.verbose);


            /*CompactedDBG<myColors> cdbg(opt.k, opt.g);
            myColors* data = nullptr;

            string seq1 = "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
            string seq2 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC";
            string seq3 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG";
            string seq4 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
            string seq5 = "GAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

            cdbg.add(seq1, true);
            cdbg.add(seq2, true);

            UnitigMap<myColors> um1 = cdbg.find(Kmer(seq1.c_str()));
            UnitigMap<myColors> um2 = cdbg.find(Kmer(seq2.c_str()));

            data = um1.getData();
            data->add(1);

            data = um2.getData();
            data->add(2);

            for (auto unitig : cdbg){

                cout << unitig.toString() << endl;
                unitig.getData()->print();
            }*/
        }
    }
}
