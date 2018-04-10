#ifndef BFG_COMMON_HPP
#define BFG_COMMON_HPP

#include <cassert>
#include <algorithm>
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

inline string reverse_complement(const string& s){

    string seq(s);

    reverse(seq.begin(), seq.end());

    for (size_t i = 0; i < seq.length(); ++i){

        switch (seq[i]){

            case 'a':
            case 'A':
                seq[i] = 'T';
                break;
            case 'c':
            case 'C':
                seq[i] = 'G';
                break;
            case 'g':
            case 'G':
                seq[i] = 'C';
                break;
            case 't':
            case 'T':
                seq[i] = 'A';
                break;
            default:
                return string();
        }
    }

    return seq;
}

inline size_t rndup(size_t v) {

    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;

    return v;
}

inline uint32_t rndup(uint32_t v) {

    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;

    return v;
}

template<typename T> class wrapperData {

    public:

        inline const T* getData() const { return &data; }
        inline T* getData() { return &data; }

    private:

        T data;
};

template<> class wrapperData<void> {

    public:

        inline const void* getData() const { return nullptr; }
        inline void* getData() { return nullptr; }
};

#endif // BFG_COMMON_HPP
