#ifndef BFG_COMMON_HPP
#define BFG_COMMON_HPP

#include <cassert>
#include <stdint.h>

using namespace std;

#define BFG_VERSION "0.2"
#define BFG_BUG_EMAIL "guillaumeholley@gmail.com"

static const char alpha[4] = {'A','C','G','T'};

inline size_t cstrMatch(const char* a, const char* b) {

    const char* a_ = a;

    while ((*a != '\0') && (*a == *b)){ ++a; ++b; }

    return a - a_;
}

inline size_t stringMatch(const string& a, const string& b, const size_t pos) {

    return distance(a.begin(), mismatch(a.begin(), a.end(), b.begin() + pos).first);
}

#endif // BFG_COMMON_HPP
