#include <string>
#include <cstdio>
#include <stdint.h>
#include <assert.h>
#include <sstream>
#include <algorithm>
#include "CompressedCoverage.hpp"

using namespace std;

// use:  cc = CompressedCoverage(sz, full);
// post: if (sz > 0) then initialize the instance else skip initializing, if full is true, initialize as full regardlesss of sz.
CompressedCoverage::CompressedCoverage(size_t sz, bool full) {

    if (sz > 0) initialize(sz, full);
    else asBits = full ? fullMask : tagMask;
}

// use:  delete cc;
// post:
CompressedCoverage::~CompressedCoverage() {

    releasePointer();
}

CompressedCoverage::CompressedCoverage(const CompressedCoverage& o) {

    if (((o.asBits & fullMask) == fullMask) || ((o.asBits & tagMask) == tagMask)) asBits = o.asBits;
    else {

        size_t sz = o.size();

        asPointer = new uint8_t[8+round_to_bytes(sz)];

        *(get32Pointer()) = sz;
        *(get32Pointer() + 1) = sz;

        memcpy(asPointer + 8, o.asPointer + 8, round_to_bytes(sz)); // 0 out array allocated
    }
}

CompressedCoverage::CompressedCoverage(CompressedCoverage&& o){

    asBits = o.asBits;
    o.asBits = o.isFull() ? fullMask : tagMask;
}

CompressedCoverage& CompressedCoverage::operator=(const CompressedCoverage& o){

    if (((o.asBits & fullMask) == fullMask) || ((o.asBits & tagMask) == tagMask)) asBits = o.asBits;
    else {

        size_t sz = o.size();

        releasePointer();
        asPointer = new uint8_t[8+round_to_bytes(sz)];

        *(get32Pointer()) = sz;
        *(get32Pointer() + 1) = sz;

        memcpy(asPointer + 8, o.asPointer + 8, round_to_bytes(sz)); // 0 out array allocated
    }

    return *this;
}

CompressedCoverage& CompressedCoverage::operator=(CompressedCoverage&& o){

    if (this != &o) {

        releasePointer();

        asBits = o.asBits;
        o.asBits = o.isFull() ? fullMask : tagMask;
    }

    return *this;
}


// use:  cc.initialize(sz, full);
// post: the data structure has been initialized either as a small local array on the stack
//       or a bigger array on the heap if sz > size_limit
void CompressedCoverage::initialize(const size_t sz, const bool full) {

    if (full) asBits = fullMask | (sz << 32);
    else if (sz <= size_limit) asBits = tagMask | (sizeMask & (sz << 2));
    else {

        asPointer = new uint8_t[8+round_to_bytes(sz)];

        *(get32Pointer()) = sz;
        *(get32Pointer() + 1) = sz;

        memset(asPointer + 8, 0, round_to_bytes(sz)); // 0 out array allocated
    }

    assert(sz == size());

}

// use:  cc.releasePointer();
// post: if there was data on the heap then it has been freed
void CompressedCoverage::releasePointer() {

    if (((asBits & tagMask) != tagMask) && ((asBits & fullMask) != fullMask)) {
        // release pointer
        uint8_t* ptr = get8Pointer();

        uintptr_t* change = &asBits;

        uintptr_t oldval = asBits;
        uintptr_t newval = fullMask | (size() << 32);

        bool replaced = false;

        while (oldval != newval){

            replaced = __sync_bool_compare_and_swap(change, oldval, newval);
            oldval = *change;
        }

        if (replaced) delete[] ptr;
    }
}


// use:  i = cc.size();
// post: i is the number of kmers that cc can hold coverage for
size_t CompressedCoverage::size() const {

    if ((asBits & tagMask) == tagMask) return ((asBits & sizeMask) >> 2);
    if ((asBits & fullMask) == fullMask) return asBits >> 32;
    return *(get32Pointer());
}


// use:  s = cc.toString();
// post: s contains all important information about cc
string CompressedCoverage::toString() const {

    bool isPtr = ((asBits & tagMask) == 0);
    size_t sz = size();
    bool full = isFull();
    uintptr_t one(1);

    string bits(64, '0');

    for (int i = 0; i < 64; i++) {

        if (asBits & (one << (63-i))) bits[i] = '1';
    }

    ostringstream info;

    if (isPtr) {

        info << "Pointer: ";

        if (full) {

            info << "Full, size = ";
            info << sz;
            info << endl;
        }
        else {

            const uint32_t filled = *(getConst32Pointer() + 1);

            info << "Non-full, size = " << sz << ", not-filled = ";
            info << filled << endl;
            info << "[";

            for (int i = 0; i < sz; i++) {

                if (i > 0) info << ", ";

                info << (int)covAt(i);
            }

            info << "] " << endl;
        }
    }
    else {

        info << "Local array:";

        if (full) info << ", Full,";

        info <<" size = " << sz;

        info << endl <<  "[";

        for (int i = 0; i < sz; i++) {

            if (i > 0) info << ", ";

            info << (int)covAt(i);
        }

        info << "] " << endl;
    }

    return bits + "\n" + info.str();
}


// use:  cc.cover(start, end);
// pre:  0 <= start , end < cc.size()
// post: the coverage of kmers: start,...,end (or reverse)
//       has been increased by one
void CompressedCoverage::cover(size_t start, size_t end) {

    if (end < start) std::swap(start, end);

    assert(end < size());

    if (isFull()) return;
    else if ((asBits & tagMask) == tagMask) { // local array

        uintptr_t s = 0x3;
        uintptr_t val = asBits;

        start = 8 + 2 * start;
        end = 8 + 2 * end;

        s <<= start;
        val >>= start;

        for (uintptr_t val_tmp; start <= end; start += 2, s <<= 2, val >>= 2) {

            val_tmp = val & 0x3;

            if (val_tmp < cov_full) {

                val_tmp = (val_tmp + 1) << start;

                uintptr_t* change = &asBits;
                uintptr_t oldval = asBits;
                uintptr_t newval = (oldval & ~s) | val_tmp;

                while (!__sync_bool_compare_and_swap(change, oldval, newval)) {

                    oldval = *change;
                    val_tmp = (oldval >> start) & 0x3;
                    val_tmp = (val_tmp < cov_full ? val_tmp + 1 : val_tmp) << start;
                    newval = (oldval & ~s) | val_tmp;
                }
            }
        }

        if (isFull()) asBits = fullMask | (size() << 32);
    }
    else {

        const uint8_t s = 0x3;
        uint8_t* ptr = get8Pointer() + 8;
        size_t fillednow = 0;
        size_t /*val,*/ index, pos;

        uint8_t val;

        for (; start <= end; ++start) {

            index = start >> 2; // start / 4
            pos = 2 * (start & 0x3); // 2 * (start % 4)
            val = (ptr[index] >> pos) & 0x3; //(ptr[index] & (s << pos)) >> pos;

            if (val < cov_full) {

                ++val;
                if (val == cov_full) ++fillednow;
                val <<= pos;

                uint8_t* change = &ptr[index];
                uint8_t oldval = ptr[index];
                uint8_t newval = (oldval & ~(s << pos)) | val;

                while (!__sync_bool_compare_and_swap(change, oldval, newval)) {

                    oldval = *change;
                    val = (oldval >> pos) & 0x3;

                    if (val < cov_full) val = (val + 1) << pos;
                    else {

                        val <<= pos;
                        --fillednow;
                    }

                    newval = (oldval & ~(s << pos)) | val;
                }
            }

            /*if (val < cov_full) {

                ++val;

                uint8_t* change = &ptr[index];

                while (1) {

                    uintptr_t oldval = *change;

                    if ((oldval & (s << pos)) >> pos == 2) break;

                    uintptr_t newval = oldval & ~(s << pos);

                    newval |= ( val << pos);

                    if (__sync_bool_compare_and_swap(change, oldval, newval)) {

                        if (((newval & (s << pos)) >> pos)  == 2) ++fillednow;

                        break;
                    }
                }
            }*/
        }

        // Decrease filledcounter
        if (fillednow > 0) __sync_add_and_fetch(get32Pointer() + 1, -fillednow);
        if (isFull()) releasePointer();
    }
}


// use:  k = cc.covat(index);
// pre:  0 <= index < size(), cc is not full
// post: k is the coverage at index
uint8_t CompressedCoverage::covAt(const size_t index) const {

    if ((asBits & fullMask) == fullMask) return cov_full;
    else if ((asBits & tagMask) == tagMask) return ((asBits >> (8 + 2*index)) & 0x03);
    else {

        const uint8_t* ptr = get8Pointer() + 8;
        const size_t pos = 2 * (index & 0x03);
        return (ptr[index >> 2] & (0x03 << pos)) >> pos;
    }
}


// use:  (low, sum) = ccov.lowCoverage();
// pre:
// post: low is the number of kmers under coverage limtis
//       sum is the sum of these low coverages
pair<size_t, size_t> CompressedCoverage::lowCoverageInfo() const {

    if (isFull()) return make_pair(0,0);

    size_t sz = size();
    size_t low = 0;
    size_t sum = 0;

    for(size_t i=0; i<sz; ++i) {

        size_t cov = covAt(i);

        if (cov < cov_full) {

            ++low;
            sum += cov;
        }
    }

    return make_pair(low, sum);
}


// use:  v = ccov.splittingVector();
// pre:  ccov.isFull() == false
// post: v is a vector of pairs (a1,b1), (a2,b2),...,(an,bn)
//       where ai < aj if i < j
//         and bi < bj if i < j
//       these pairs are all the fully covered subintervals of the corresponding contig
//       i.e. [ai,...,bi-1] is fully covered
vector<pair<int, int>> CompressedCoverage::splittingVector() const {

    size_t a = 0, b = 0;
    const size_t sz = size();

    vector<pair<int, int>> v;

    while (b != sz) {
        // [a,...,b-1] is a fully covered subinterval and (a,b) has been added to v
        while ((a < sz) && (covAt(a) < cov_full)) ++a;

        if (a == sz) break;

        b = a;

        while ((b < sz) && (covAt(b) >= cov_full)) ++b;

        v.push_back(make_pair(a,b));

        a = b;
    }

    return v;
}


// use:  b = cc.isFull();
// post: (b == true) <==> cc is full
bool CompressedCoverage::isFull() const {

    if ((asBits & fullMask) == fullMask) return true;
    if ((asBits & tagMask) == tagMask) return (asBits >> 8) == (localCoverageMask >> 2*(28 - size()));

    return *(getConst32Pointer() + 1) == 0;
}

// use: cc.setFull()
// pre:
// post: cc is full and any memory is released
void CompressedCoverage::setFull() {

    if ((asBits & fullMask) != fullMask){

        if ((asBits & tagMask) == tagMask) asBits = fullMask | (size() << 32);
        else releasePointer();
    }

    return;
}

size_t CompressedCoverage::cov_full = 2;

uintptr_t CompressedCoverage::localCoverageMask = 0xAAAAAAAAAAAAAA; // 0b10101010101010101010101010101010101010101010101010101010
