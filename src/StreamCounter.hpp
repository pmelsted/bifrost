#ifndef BIFROST_STREAMCOUNTER_HPP
#define BIFROST_STREAMCOUNTER_HPP

#include <sstream>
#include <stdint.h>
#include <string>
#include <cmath>
#include <fstream>
#include <assert.h>

class StreamCounter {

    public:

        StreamCounter(double e_, int seed_ = 0) : MAX_TABLE(32), maxVal(3ULL), countWidth(2), countsPerLong(32), e(e_), seed(seed_), sumCount(0) {

            const size_t numcounts = max(static_cast<size_t>(48.0/(e*e) + 1), static_cast<size_t>(8192)); // approx 3 std-dev, true with 0.001 prob.

            size = (numcounts + countsPerLong - 1) / countsPerLong; // size is number of uint64_t's use
            size = roundUpPowerOfTwo(size);
            mask = (size * countsPerLong) - 1;
            maxCount = size * countsPerLong * maxVal;

            M = new size_t[MAX_TABLE]();
            table = new uint64_t[size * MAX_TABLE]();

        }

        StreamCounter(const StreamCounter& o) : seed(o.seed), e(o.e), size(o.size), maxCount(o.maxCount), mask(o.mask),
                                                MAX_TABLE(o.MAX_TABLE), countWidth(o.countWidth), countsPerLong(o.countsPerLong),
                                                maxVal(o.maxVal), sumCount(o.sumCount) {
            // copy constructor, creates object of same size with same seed, but empty data
            M = new size_t[MAX_TABLE]();
            table = new uint64_t[size * MAX_TABLE]();

            memcpy(M, o.M, MAX_TABLE * sizeof(size_t));
            memcpy(table, o.table, size * MAX_TABLE * sizeof(uint64_t));
        }

        ~StreamCounter() {

            delete[] table;
            delete[] M;
        }

        inline int getSeed() const { return seed; }

        void operator()(const uint64_t hashval) {

            ++sumCount;

            // hashval is XXX .. XXX1000.. (w times) ..00 0
            size_t w = __builtin_ctzll(hashval);

            if (w >= MAX_TABLE) w = MAX_TABLE - 1;
            if (M[w] == size * countsPerLong * maxVal) return;

            // hashval is now XXXX...XX, random
            const uint64_t hval = hashval >> (w + 1); // shift away pattern of XXX10000
            const uint64_t index = hval & mask;
            const uint64_t val = getVal(index,w);

            if (val != maxVal) { // max count

                setVal(index, w, val+1);
                ++M[w];
            }
        }

        /*void update(const uint64_t h) {

            ++sumCount;

            // hashval is XXX .. XXX1000.. (w times) ..00 0
            const uint64_t w = __builtin_ctzll(h);

            if (w >= MAX_TABLE) w = MAX_TABLE - 1;
            if (M[w] == size * countsPerLong * maxVal) return;

            // hashval is now XXXX...XX, random
            const uint64_t hval = h >> (w + 1); // shift away pattern of XXX10000
            const uint64_t index = hval & mask;
            const uint64_t wordindex = w * size + (index / countsPerLong);
            const uint64_t bitshift = countWidth * (index & (countsPerLong - 1));
            const uint64_t bitmask = maxVal << bitshift;
            const uint64_t val = (table[wordindex] & bitmask) >> bitshift;

            if (val != maxVal){

                table[wordindex] = ((((val + 1) & maxVal) << bitshift) & bitmask) | (table[wordindex] & ~bitmask);
                ++M[w];
            }
        }

        void update_p(const uint64_t h) {

            ++sumCount;

            // hashval is XXX .. XXX1000.. (w times) ..00 0
            const uint64_t w = __builtin_ctzll(h);

            if (w >= MAX_TABLE) w = MAX_TABLE - 1;
            if (M[w] == size * countsPerLong * maxVal) return;

            // hashval is now XXXX...XX, random
            const uint64_t hval = h >> (w + 1); // shift away pattern of XXX10000
            const uint64_t index = hval & mask;

            const uint64_t wordindex = w * size + (index / countsPerLong);
            const uint64_t bitshift = countWidth * (index & (countsPerLong - 1));
            const uint64_t bitmask = maxVal << bitshift;

            uint64_t old_val = table[wordindex];
            uint64_t val = ((old_val & bitmask) >> bitshift);
            uint64_t new_val = ((((val + 1) & maxVal) << bitshift) & bitmask) | (old_val & ~bitmask);

            while ((val != maxVal) && !table[wordindex].compare_exchange_weak(old_val, new_val)){

                old_val = table[wordindex];
                val = ((old_val & bitmask) >> bitshift);
                new_val = ((((val + 1) & maxVal) << bitshift) & bitmask) | (old_val & ~bitmask);
            }

            if (val != maxVal) {

                old_val = M[w];
                new_val = old_val + 1;

                while (!M[w].compare_exchange_weak(old_val, new_val)) {

                    old_val = M[w];
                    new_val = old_val + 1;
                }
            }
        }*/

        bool join(const StreamCounter& o) {

            if (size != o.size || seed != o.seed) return false;

            for (size_t i = 0; i < MAX_TABLE; ++i) {

                M[i] = 0;

                for (size_t j = 0; j < size * countsPerLong; ++j) {

                    setVal(j, i, getVal(j,i) + o.getVal(j,i));

                    M[i] += getVal(j,i);
                }
            }

            sumCount += o.sumCount;

            return true;
        }

        size_t F0() const {

            int n = 0;
            double sum = 0;
            double limit = 0.2;

            const size_t R = size * countsPerLong;

            while (n == 0 && limit > 1e-8) {

                for (size_t i = 0; i < MAX_TABLE; ++i) {

                    size_t ts = 0;

                    for (size_t j = 0; j < R; ++j) {

                        if (getVal(j,i) > 0) ++ts;
                    }

                    if (ts == 0) { // nothing in this level, hack to break out of loop

                        limit = 0;
                        break;
                    }

                    if ((ts <= (1 - limit) * R) && ((ts >= limit * R) || (i == 0))) {

                        sum += (log(1.0 - ts/((double) R)) / log(1.0-1.0 / R)) * pow(2.0, i+1);
                        ++n;
                        break;
                    }
                }

                limit = limit / 1.5;
            }

            if (n == 0) return 0;
            return (size_t)(sum/n);
        }

        inline size_t F1() const { return sumCount; }

        size_t f1() const {

            int n = 0;

            double sum = 0;
            double limit = 0.2;

            const size_t R = size * countsPerLong;

            while (n == 0 && limit > 1e-8) {

                for (size_t i = 0; i < MAX_TABLE; ++i) {

                    size_t r1 = 0, r0 = 0;

                    for (size_t j = 0; j < R; ++j) {

                        const uint64_t val = getVal(j,i);

                        if (val == 0) ++r0;
                        if (val == 1) ++r1;

                    }

                    if (r0 == R) { // empty level

                        limit = 0;
                        break;
                    }

                    if (((r0 <= ((1 - limit) * R)) || (i == 0)) && (r0 >= (limit * R))) {
                        // i==0 takes care of small first levels where R is too large
                        sum += (R-1) * (r1/((double) r0)) * pow(2.0,i+1);
                        ++n;
                        break;
                    }
                }

                limit = limit/1.5;
            }

            if (n == 0) return 0;
            return (size_t) (sum/n);
        }

        std::string humanReport() const {

            std::stringstream s;

            const size_t eF0 = F0();
            const size_t ef1 = f1();
            const size_t eF1 = sumCount;

            s << readable(eF0-ef1) << " repeated, " << readable(eF0) << " distinct, " <<
            readable(ef1) << " singletons, " << readable(eF1) << " total k-mers processed";

            return s.str();
        }

        std::string readable(size_t x) const {

            std::stringstream s;

            if (x < (1ULL<<10)) s << x;
            else if (x < (1ULL<<20)) s << (x/1024) << "K";
            else if (x < (1ULL<<30)) s << (x/(1ULL<<20)) << "M";
            else s << (x/(1ULL<<30)) << "G";

            return s.str();

        }

        std::string report(bool useTSV = false) const {

            std::stringstream s;

            if (useTSV) s << F0() << "\t" << f1() << "\t" << F1() << std::endl;
            else {

                s << "F0 = " << F0() << std::endl;
                s << "f1 = " << f1() << std::endl;
                s << "F1 = " << sumCount << std::endl;
            }

            return s.str();
        }

        bool writeBinary(const std::string& fn) {

            std::ofstream out;

            out.open(fn.c_str(), std::ios::out | std::ios::binary);

            if (!out.is_open()) {

                std::cerr << "Error: could not write to file " << fn << std::endl;
                return false;
            }

            // 1. write out seed
            out.write((char*)&seed, sizeof(seed));
            // 2. write out size
            out.write((char*)&size, sizeof(size));
            // 3. write out e
            out.write((char*)&e, sizeof(e));
            // 3.5 write out sumCount
            out.write((char*)&sumCount, sizeof(sumCount));
            // 4. write out MAX_TABLE
            out.write((char*)&MAX_TABLE, sizeof(MAX_TABLE));
            // 5. write out M
            out.write((char*)M, MAX_TABLE*sizeof(M[0]));
            // 6. write out table
            out.write((char*)table, size*MAX_TABLE*sizeof(table[0]));

            out.flush();
            out.close();

            return true;
        }

        bool loadBinary(const std::string& fn) {

            std::ifstream in;
            const size_t oldsize = size;

            in.open(fn.c_str(), std::ios::in | std::ios::binary);

            // 1. read seed
            in.read((char*)&seed,sizeof(seed));
            // 2. read size
            in.read((char*)&size,sizeof(size));
            // 3. read e
            in.read((char*)&e, sizeof(e));
            // 3.5 read sumCount
            in.read((char*)&sumCount, sizeof(sumCount));

            size_t max_table;
            // 4. read MAX_TABLE
            in.read((char*)&max_table, sizeof(max_table));

            if(MAX_TABLE != max_table) {

                std::cerr <<"Error: Max table size doesn't match" << std::endl;
                std::cerr << "MAX_TABLE = " << MAX_TABLE << std::endl << "max_table = " << max_table << std::endl;
                exit(1);
            }

            // fill in other variables
            mask = (size * countsPerLong) - 1;
            maxCount = size * countsPerLong * maxVal;

            // allocate space
            if (oldsize != size) {

                if (M != 0) {

                    delete[] M;
                    M = new size_t[MAX_TABLE];
                }

                if (table != 0) {

                    delete[] table;
                    table = new uint64_t[size * MAX_TABLE];
                }
            }

            // 5. read M
            in.read((char*)M, MAX_TABLE * sizeof(M[0]));

            // 6. read T
            in.read((char*)table, size * MAX_TABLE * sizeof(table[0]));
            in.close();

            return true;
        }

    private:

        static size_t roundUpPowerOfTwo(size_t size) {

            --size;

            size |= size >> 1;
            size |= size >> 2;
            size |= size >> 4;
            size |= size >> 8;
            size |= size >> 16;
            size |= size >> 32;

            ++size;

            return size;
        }

        uint64_t getVal(const size_t index, const size_t w) const {

            const size_t wordindex = w * size + (index / countsPerLong);
            const size_t bitindex = index & (countsPerLong - 1) ;
            const uint64_t bitmask = maxVal << (countWidth * bitindex);

            return (table[wordindex] & bitmask) >> (countWidth * bitindex);
        }

        void setVal(const size_t index, const size_t w, uint64_t val) {

            if (val > maxVal) val = maxVal;

            const size_t wordindex = w * size + (index / countsPerLong);
            const size_t bitindex = index & (countsPerLong - 1);
            const uint64_t bitmask = maxVal << (countWidth * bitindex);

            table[wordindex] = (((val & maxVal) << (countWidth * bitindex)) & bitmask) | (table[wordindex] & ~bitmask);
        }

        /*void incVal_p(const size_t index, const size_t w) {

            const size_t wordindex = w * size + (index / countsPerLong);
            const size_t bitshift = countWidth * (index & (countsPerLong - 1));
            const uint64_t bitmask = maxVal << bitshift;

            uint64_t old_val = table[wordindex];
            uint64_t val = ((old_val & bitmask) >> bitshift);
            uint64_t new_val = ((((val + 1) & maxVal) << bitshift) & bitmask) | (old_val & ~bitmask);

            while ((val != maxVal) && !table[wordindex].compare_exchange_weak(old_val, new_val)){

                old_val = table[wordindex];
                val = ((old_val & bitmask) >> bitshift);
                new_val = ((((val + 1) & maxVal) << bitshift) & bitmask) | (old_val & ~bitmask);
            }

            if (val != maxVal) {

                old_val = M[w];
                new_val = old_val + 1;

                while (!M[w].compare_exchange_weak(old_val, new_val)) {

                    old_val = M[w];
                    new_val = old_val + 1;
                }
            }
        }

        void incVal(const size_t index, const size_t w) {

            const size_t wordindex = w * size + (index / countsPerLong);
            const size_t bitshift = countWidth * (index & (countsPerLong - 1));
            const uint64_t bitmask = maxVal << bitshift;
            const uint64_t val = (table[wordindex] & bitmask) >> bitshift;

            if (val != maxVal){

                table[wordindex] = ((((val + 1) & maxVal) << bitshift) & bitmask) | (table[wordindex] & ~bitmask);
                ++M[w];
            }
        }*/

        int seed;
        double e;
        uint64_t* table;
        size_t sumCount;
        size_t *M;
        size_t size;
        size_t maxCount;
        uint64_t mask;
        const size_t MAX_TABLE;
        const size_t countWidth; // number of bits per count, even number
        const size_t countsPerLong;  // fix this
        const uint64_t maxVal; // has to be a power of 2-1
};

#endif
