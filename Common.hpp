#ifndef BFG_COMMON_HPP
#define BFG_COMMON_HPP

#include <cassert>
#include <stdint.h>

#ifndef __clang__ 
  #include <omp.h>
#endif

using namespace std;

#define BFG_VERSION "0.1"
#define BFG_BUG_EMAIL "pmelsted@gmail.com"

static const char alpha[4] = {'A','C','G','T'};

#endif // BFG_COMMON_HPP
