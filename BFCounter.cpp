

#include <iostream>


#include "CountBF.hpp"
#include "DumpBF.hpp"
#include "Common.hpp"

using namespace std;

void PrintUsage() {
  cerr << "BFCounter " << BFC_VERSION << endl << endl;
  cerr << "A memory efficient k-mer counting program." << endl << endl;
  cerr << "Usage: BFCounter <cmd> [options] ..." << endl << endl;
  cerr << "Where <cmd> can be one of:" << endl;
  cerr << 
    "    count        Counts the occurrences of k-mers in sequence files" << endl <<
    "    dump         Writes occurrences of k-mers into a tab-separated text file" << endl << 
    "    cite         Prints information for citing the paper" << endl <<
    "    version      Displays version number" << endl << endl;
    ;


}

void PrintVersion() {
  cerr <<  BFC_VERSION << endl;
}

void PrintCite() {
  cerr << "The paper describing this software has been published." << endl;
  cerr << "When using this program in your research, please cite" << endl << endl;
  cerr << "Melsted and Pritchard: Efficient counting of k-mers in DNA sequences using a bloom filter.";
  cerr << "BMC Bioinformatics 2011 12:333." << endl << endl;
  cerr << "The full paper is available online with open-access at" << endl << endl;;
  cerr << "http://www.biomedcentral.com/1471-2105/12/333" << endl << endl;
}

int main(int argc, char **argv) {

  if (argc < 2) {
    cout << "Error: too few arguments" << endl;
    PrintUsage();
  } else {
    if (strcmp(argv[1], "cite") == 0) {
      PrintCite();
    } else if (strcmp(argv[1], "version") == 0) {
      PrintVersion();
    } else if (strcmp(argv[1], "count") == 0) {
      CountBF(argc-1,argv+1);
    } else if (strcmp(argv[1], "dump") == 0) {
      DumpBF(argc-1,argv+1);
    } else {
      cout << "Did not understand command " << argv[1] << endl;
      PrintUsage();
    }
  }
}
