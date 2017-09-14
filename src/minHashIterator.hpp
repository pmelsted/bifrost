#ifndef MINHASHITERATOR_H
#define MINHASHITERATOR_H

#include <stdint.h>
#include <deque>
#include <vector>

#include <iostream>
using namespace std;

//using std::vector;

struct minHashResult {

	minHashResult() : hash((uint64_t) -1),pos(-1) {}
	minHashResult(const uint64_t h, const int p) : hash(h), pos(p) {}
	minHashResult(const minHashResult& o) : hash(o.hash), pos(o.pos) {}
	uint64_t hash;
	int pos;
};

template<class HF>
struct minHashResultIterator;

// TODO: return current position as well as minimum
// templated?
template<class HF>
class minHashIterator {

    public:

        minHashIterator(int _k, int _g, HF _h) : s(NULL), n(0), k(_k), g(_g), hf(_h), v(k-g+1), p(-1), invalid(true), nh(false) {
            hf.setK(g);
        }

        minHashIterator(const char* _s, int _length, int _k, int _g, HF _h, bool _nh) : k(_k), g(_g), hf(_h), v(k-g+1), p(-1), invalid(true), nh(_nh) {
            hf.setK(g);
            initString(_s,_length);
        }

        minHashIterator() : s(NULL), n(0), k(0), g(0), hf(HF(0)), invalid(true), nh(false) {}

        void seed(HF &ohf) {
            hf = ohf;
            invalid=true;
        }

        void initString(const char* _s, int _length) {

            n = _length;
            s = _s;
            p = -1;// reset position
            invalid = false;
            v.clear();

            // advance to first position or set to invalid
            if (n < k || k < g) invalid = true;
            else operator++();
        }


        bool operator==(const minHashIterator& o) {

            if (invalid || o.invalid) return invalid && o.invalid;
            return s==o.s && n==o.n && g==o.g && k==o.k && nh==o.nh;
        }

        bool operator!=(const minHashIterator& o) { return !this->operator==(o); }

        // invariant:
        //   v[0..sz] contains only ascending elements in s from [p,p+k-g+1)
        //
        //   there exists no a<b s.t. v[a] > v[b].
        minHashIterator& operator++() {

            if (invalid) return *this;

            ++p; // advance to next k-mer

            if (p >= n-k+1 || s[p+k-1] == 0) {
                // out of bounds
                invalid = true;
                return *this;
            }

            const int shift = nh ? 1 : 0;

            if (p==0) {

                hf.init(&s[shift]);

                v.push_back(minHashResult(hf.hash(), shift));

                for (int j = shift; j < k-g-shift;) {

                    hf.update(s[j],s[j+g]);

                    uint64_t h = hf.hash();
                    int t = ((int)v.size())-1;

                    while (t >= 0 && v[t].hash > h) {

                        v.pop_back();
                        t--;
                    }

                    j++;
                    v.push_back(minHashResult(h,j));
                }
            }
            else {

                if (v[0].pos < p + shift) v.pop_front(); // remove first element, fell outside of window

                hf.update(s[p+k-g-1-shift],s[p+k-1-shift]);

                uint64_t h = hf.hash();
                int t = ((int) v.size())-1;

                while (t >= 0 && v[t].hash > h) {

                    v.pop_back();
                    t--;
                }

                v.push_back(minHashResult(h,p+k-g-shift));
            }

            return *this;
        }

        minHashIterator operator++(int) {

            minHashIterator tmp(*this);
            operator++();

            return tmp;
        }

        minHashIterator& operator+=(int i) {

            for (; i > 0; i--) operator++();
            return *this;
        }

        minHashResult getNewMin(const minHashResult& mhr_discard) const {

            if (invalid) return minHashResult();

            const int shift = nh ? 1 : 0;
            const int end = p+k-g-shift;

            int j = p + shift;

            uint64_t h;

            HF hf_tmp;

            hf_tmp.setK(g);
            hf_tmp.init(&s[j]);

            while ((hf_tmp.hash() <= mhr_discard.hash) && (j < end)){

                hf_tmp.update(s[j],s[j+g]);
                j++;
            }

            if ((j == end) && (hf_tmp.hash() <= mhr_discard.hash)) return /*v[0]*/mhr_discard; //No other minimizer can be found

            minHashResult mhr = minHashResult(hf_tmp.hash(), j);

            while (j < end) {

                hf_tmp.update(s[j],s[j+g]);
                j++;

                h = hf_tmp.hash();

                if ((h <= mhr.hash) && (h > mhr_discard.hash)){

                    if (h == mhr.hash){ // break ties if two minimizers have same hash

                        if (Minimizer(&s[j]).rep() < Minimizer(&s[mhr.pos]).rep()){

                            mhr.hash = hf_tmp.hash();
                            mhr.pos = j;
                        }
                    }
                    else {

                        mhr.hash = hf_tmp.hash();
                        mhr.pos = j;
                    }
                }
            }

            return mhr;
        }

        minHashResultIterator<HF> operator*() const {return minHashResultIterator<HF>(this);}

        inline uint64_t getHash() const { return invalid ? 0 : v[0].hash; }

        inline int getPosition() const { return invalid ? 0 : v[0].pos; }

        inline int getKmerPosition() const { return p; }

        const char *s; //Minimizers are from k-mers, k-mers are from a sequence s
        int n; // Length of sequence s
        int k; // Length of k-mers
        int g; // Length of minimizers
        HF hf; // Rolling hash function
        deque<minHashResult> v; //Hash values and positions of a same minimizer with k-mer at position p
        int p; // Position of current k-mer traversed in the sequence
        bool invalid; // If sequence is invalid (iterating on k-mers out of bounds, etc.)
        bool nh; // If true, minimizer of k-mers km cannot start at position 0 or k-g

        // private copy constructor
        minHashIterator(const minHashIterator& o) : s(o.s), n(o.n), k(o.k), g(o.g), hf(o.hf), v(o.v), p(o.p), invalid(o.invalid), nh(o.nh) {}
};

template <class HF>
struct minHashResultIterator {

    minHashResultIterator(const minHashIterator<HF> *p) : p(p), invalid(false), pos(0), p_pos(p->p), p_s(p->s) {}
    minHashResultIterator() : invalid(true), p_pos(-1), p_s(NULL) {}

    bool operator==(const minHashResultIterator& o) {

        if (o.invalid || invalid) return o.invalid && invalid;
        return p_pos == o.p_pos && p_s == o.p_s && pos==o.pos;

    }

    bool operator!=(const minHashResultIterator& o) {
        return !this->operator==(o);
    }

    minHashResultIterator& operator++() {

        if (invalid) return *this;

        if ((p_pos != p->p || p_s != p->s) // check if parent iterator has moved
            || (pos>=p->v.size()-1) // or if we advance past the end
            || (p->v[pos+1].hash != p->v[pos].hash))    { // or if the next position doesn't share the hash value
            // invalidate the iterator
            invalid = true;
            return *this;
        }
        else pos++; // advance to next equal hash value

        return *this;
    }


    minHashResultIterator operator++(int) {

        minHashResultIterator tmp(*this);
        operator++();

        return tmp;
    }

    const minHashResult& operator*() const { return p->v[pos]; }
    const minHashResult* operator->() const { return &(p->v[pos]); }

    // pos points to a minHashIterator, all the values from p.v[0] to p.v[pos] have the
    // same (minimum) hash value or the this.invalid is true
    // at the time this was created p.s==p_s and p.p=p_pos
    const minHashIterator<HF> *p;
    bool invalid;
    int pos;
    const int p_pos;
    const char *p_s;
};


template<class HF>
struct preAllocMinHashResultIterator;

template<class HF>
class preAllocMinHashIterator {

    public:

        preAllocMinHashIterator() : s(NULL), n(0), k(0), g(0), hf(HF(0)), p(-1), v(0), p_cur_start(0), p_cur_end(0), invalid(true), nh(false) {}

        preAllocMinHashIterator(const char* _s, int _n, int _k, int _g, HF _h, bool _nh) :
                                s(_s), n(_n), k(_k), g(_g), hf(_h), p(-1), p_cur_start(0), p_cur_end(0), invalid(true), nh(_nh) {

            if ((s != NULL) && (n >= k) && (k >= g)){

                invalid = false;

                v = vector<minHashResult>(n-g+1);

                hf.setK(g);
                operator++();
            }
        }

        bool operator==(const preAllocMinHashIterator& o) {

            if (invalid || o.invalid) return invalid && o.invalid;
            return s==o.s && n==o.n && g==o.g && k==o.k && nh==o.nh;
        }

        bool operator!=(const preAllocMinHashIterator& o) { return !this->operator==(o); }

        // invariant:
        //   v[0..sz] contains only ascending elements in s from [p,p+k-g+1)
        //
        //   there exists no a<b s.t. v[a] > v[b].
        preAllocMinHashIterator& operator++() {
            //cerr << "operator++";
            if (invalid) return *this;

            ++p; // advance to next k-mer

            if (p >= n-k+1 || s[p+k-1] == 0) {
                // out of bounds
                invalid = true;
                return *this;
            }

            const int shift = nh ? 1 : 0;

            if (p == 0) {

                hf.init(&s[shift]);

                v[p_cur_end] = minHashResult(hf.hash(), shift);
                p_cur_end++;

                for (int j = shift; j < k-g-shift;) {

                    hf.update(s[j], s[j+g]);

                    uint64_t h = hf.hash();

                    while (p_cur_end > p_cur_start && v[p_cur_end-1].hash > h) p_cur_end--;

                    j++;

                    v[p_cur_end] = minHashResult(h,j);
                    p_cur_end++;
                }
            }
            else {

                if (v[p_cur_start].pos < p + shift) p_cur_start++; // remove first element, fell outside of window

                hf.update(s[p+k-g-1-shift], s[p+k-1-shift]);

                uint64_t h = hf.hash();

                while (p_cur_end > p_cur_start && v[p_cur_end-1].hash > h) p_cur_end--;

                v[p_cur_end] = minHashResult(h,p+k-g-shift);
                p_cur_end++;
            }

            return *this;
        }

        preAllocMinHashIterator operator++(int) {

            preAllocMinHashIterator tmp(*this);
            operator++();

            return tmp;
        }

        preAllocMinHashIterator& operator+=(int i) {

            for (; i > 0; i--) operator++();
            return *this;
        }

        preAllocMinHashResultIterator<HF> operator*() const {return preAllocMinHashResultIterator<HF>(*this);}

        inline uint64_t getHash() const { return invalid ? 0 : v[p_cur_start].hash; }

        inline int getPosition() const { return invalid ? 0 : v[p_cur_start].pos; }

        inline int getNbMin() const { return p_cur_end - p_cur_start; }

        inline int getKmerPosition() const { return p; }

        minHashResult getNewMin(const minHashResult& mhr_discard) const {

            if (invalid) return minHashResult();

            const int shift = nh ? 1 : 0;
            const int end = p+k-g-shift;

            int j = p + shift;

            uint64_t h;

            HF hf_tmp;

            hf_tmp.setK(g);
            hf_tmp.init(&s[j]);

            while ((hf_tmp.hash() <= mhr_discard.hash) && (j < end)){

                hf_tmp.update(s[j],s[j+g]);
                j++;
            }

            if ((j == end) && (hf_tmp.hash() <= mhr_discard.hash)) return /*v[p_cur_start]*/ mhr_discard;

            minHashResult mhr = minHashResult(hf_tmp.hash(), j);

            while (j < end) {

                hf_tmp.update(s[j],s[j+g]);
                j++;

                h = hf_tmp.hash();

                if ((h <= mhr.hash) && (h > mhr_discard.hash)){

                    if (h == mhr.hash){ // break ties if two minimizers have same hash

                        if (Minimizer(&s[j]).rep() < Minimizer(&s[mhr.pos]).rep()){

                            mhr.hash = hf_tmp.hash();
                            mhr.pos = j;
                        }
                    }
                    else {

                        mhr.hash = hf_tmp.hash();
                        mhr.pos = j;
                    }
                }
            }

            if (mhr.hash <= mhr_discard.hash) cerr << "Problem here" << endl;

            return mhr;
        }

        const char *s; //Minimizers are from k-mers, k-mers are from a sequence s
        int n; // Length of sequence s
        int k; // Length of k-mers
        int g; // Length of minimizers
        HF hf; // Rolling hash function
        vector<minHashResult> v; //Hash values and positions of a same minimizer with k-mer at position p
        size_t p_cur_start;
        size_t p_cur_end;
        int p; // Position of current k-mer traversed in the sequence
        bool invalid; // If sequence is invalid (iterating on k-mers out of bounds, etc.)
        bool nh;

        // private copy constructor
        preAllocMinHashIterator(const preAllocMinHashIterator& o) : s(o.s), n(o.n), k(o.k), g(o.g), hf(o.hf), v(o.v), p(o.p),
                                                                    p_cur_start(o.p_cur_start), p_cur_end(o.p_cur_end), invalid(o.invalid), nh(o.nh) {}

        preAllocMinHashIterator(const preAllocMinHashIterator& o, int len) : s(o.s + o.p), n(len), k(o.k), g(o.g), hf(o.hf), p(0), p_cur_start(0),
                                                                                p_cur_end(o.p_cur_end - o.p_cur_start), invalid(o.invalid), nh(o.nh) {

            if (!invalid && (o.p + n <= o.n)){

                v = std::vector<minHashResult>(o.v.begin() + o.p_cur_start, o.v.begin() + o.p_cur_end);

                for (auto& min_h : v) min_h.pos -= o.p;
            }
            else invalid = true;
        }
};

template <class HF>
struct preAllocMinHashResultIterator {

    preAllocMinHashResultIterator(const preAllocMinHashIterator<HF>& _p) : p(&_p), p_pos(_p.p), p_it(_p.p_cur_start), p_it_end(_p.p_cur_end), p_s(_p.s), invalid(false) {}
    preAllocMinHashResultIterator() : invalid(true), p(NULL), p_pos(-1), p_s(NULL), p_it(0), p_it_end(0) {}

    bool operator==(const preAllocMinHashResultIterator& o) const {

        if (o.invalid || invalid) return o.invalid && invalid;
        return p_pos == o.p_pos && p_s == o.p_s && p_it==o.p_it && p_it_end==o.p_it_end;
    }

    bool operator!=(const preAllocMinHashResultIterator& o) const {
        return !this->operator==(o);
    }

    preAllocMinHashResultIterator& operator++() {

        if (invalid) return *this;

        if ((p_s != p->s || p_pos != p->p || p_it_end != p->p_cur_end) // check if parent iterator has moved
            || (p_it >= p_it_end - 1) // or if we advance past the end
            || (p->v[p_it + 1].hash != p->v[p_it].hash)) { // or if the next position doesn't share the hash value
            // invalidate the iterator
            invalid = true;
            return *this;
        }
        else p_it++; // advance to next equal hash value

        return *this;
    }

    preAllocMinHashResultIterator operator++(int) {

        preAllocMinHashResultIterator tmp(*this);
        operator++();

        return tmp;
    }

    preAllocMinHashResultIterator& operator=(const preAllocMinHashResultIterator &o){

        p_it = o.p_it;
        p_it_end = o.p_it_end;
        invalid = o.invalid;

        if (operator!=(o)) invalid = true;

        return *this;
    }

    const minHashResult& operator*() const { return p->v[p_it]; }
    const minHashResult* operator->() const { return &(p->v[p_it]); }

    // pos points to a minHashIterator, all the values from p.v[0] to p.v[pos] have the
    // same (minimum) hash value or the this.invalid is true
    // at the time this was created p.s==p_s and p.p=p_pos

    const preAllocMinHashIterator<HF>* p;
    bool invalid;
    size_t p_it;
    size_t p_it_end;
    const int p_pos;
    const char *p_s;
};










template <class HF>
struct minHashKmer {

    public:

        minHashKmer(const char* _s, int _k, int _g, HF _h, bool neighbor_hash) : s(_s), k(_k), g(_g), hf(_h), h(0), p(-1), invalid(true), nh(neighbor_hash) {

            if ((s != NULL) && ((n = strlen(s)) >= k) && (k >= g)){

                invalid = false;

                hf.setK(g);
                compute_min();
            }
        }

        minHashKmer() : s(NULL), n(0), k(0), g(0), hf(HF(0)), h(0), p(-1), nb(0), invalid(true), nh(false) {}

        minHashKmer(const preAllocMinHashIterator<HF>& o) : s(o.s), k(o.k), g(o.g), hf(o.hf), nh(o.nh), h(0), p(-1), nb(o.nb), invalid(o.invalid){

            if (!invalid){

                h = o.v[o.p_cur_start].hash;
                p = o.v[o.p_cur_start].pos - o.p;
            }
        }

        bool operator==(const minHashKmer& o) {

            if(invalid || o.invalid) return invalid && o.invalid;
            return s==o.s && n==o.n && g==o.g && k==o.k && nh==o.nh && h==o.h && p==o.p && nb==o.nb;
        }

        bool operator!=(const minHashKmer& o) { return !this->operator==(o); }

        uint64_t getHash() const { return h; }

        int getPosition() const { return p; }

        int getNbMin() const { return nb; }

    private:

        void compute_min(){

            if (invalid) return;

            const int shift = nh ? 1 : 0;

            hf.init(&s[shift]);

            h = hf.hash();
            p = shift;
            nb = 1;

            for (int j = shift; j < k-g-shift; j++) {

                hf.update(s[j], s[j+g]);

                if (hf.hash() <= h){

                    if (hf.hash() == h) nb++;
                    else {

                        h = hf.hash();
                        p = j + 1;
                        nb = 1;
                    }
                }
            }
        }

        const char* s;
        HF hf;
        uint64_t h;
        int n;
        int k;
        int g;
        int p;
        int nb;
        bool invalid;
        bool nh;
};

#endif // MINHASHITERATOR_H
