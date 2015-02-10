#include <iostream>
#include <cstring>

#include "BuildContigs.hpp"
#include "FilterReads.hpp"
#include "Common.hpp"

using namespace std;


// use:  PrintUsage();
// post: How to run BFGraph has been printed to cout
void PrintUsage() {
  cout << "BFGraph " << BFG_VERSION << endl << endl;
  cout << "A memory efficient de Bruijn graph assembler." << endl << endl;
  cout << "Usage: BFGraph <cmd> [options] ..." << endl << endl;
  cout << "Where <cmd> can be one of:" << endl;
  cout <<
       "    filter       Filters errors from reads" << endl <<
       "    contigs      Builds an initial contig graph" << endl <<
       "    cite         Prints information for citing the paper" << endl <<
       "    version      Displays version number" << endl << endl;
  ;
}


// use:  PrintVersion();
// post: The version of the program has been printed to cout
void PrintVersion() {
  cout <<  BFG_VERSION << endl;
}


// use:  PrintCite();
// post: Information of how to cite this software has been printed to cerr
void PrintCite() {
  cout << "The paper describing this software has not been published." << endl;
  //  cerr << "When using this program in your research, please cite" << endl << endl;
}


int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << "Error: too few arguments" << endl;
    PrintUsage();
  } else {
    if (strcmp(argv[1], "cite") == 0) {
      PrintCite();
    } else if (strcmp(argv[1], "version") == 0) {
      PrintVersion();
    } else if (strcmp(argv[1], "filter") == 0) {
      FilterReads(argc-1,argv+1);
    } else if (strcmp(argv[1], "contigs") == 0) {
      BuildContigs(argc-1,argv+1);
    } else {
      cout << "Did not understand command " << argv[1] << endl;
      PrintUsage();
    }
  }
}
