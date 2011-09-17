

#include <iostream>


#include "CountBF.hpp"
#include "DumpBF.hpp"
#include "Common.hpp"

using namespace std;

void PrintUsage() {
  cerr << "Usage: BFCounter <cmd> [options] ..." << endl;
  cerr << "Where <cmd> can be count, dump, help, version" << endl;
}

void PrintVersion() {
  cerr <<  BFC_VERSION << endl;
}

int main(int argc, char **argv) {

  if (argc < 2) {
    cout << "Error: too few arguments" << endl;
    PrintUsage();
  } else {
    if (strcmp(argv[1], "help") == 0) {
      //PrintHelp();
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
