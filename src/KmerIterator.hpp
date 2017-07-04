#ifndef BFG_KMER_ITERATOR_HPP
#define BFG_KMER_ITERATOR_HPP

#include <iostream>
#include <iterator>
#include "Kmer.hpp"


/* Short description:
 *  - Easily iterate through kmers in a read
 *  - If the read contains any N, then the N is skipped and checked whether
 *    there is a kmer to the right of the N
 *  - iter->first gives the kmer, iter->second gives the position within the reads
 * */
class KmerIterator : public std::iterator<std::input_iterator_tag, std::pair<Kmer, int>, int> {

    public:

        KmerIterator() : s_(NULL), p_(), invalid_(true) {}
        KmerIterator(const char *s) : s_(s), p_(), invalid_(false) { find_next(-1,-1,false);}
        KmerIterator(const KmerIterator& o) : s_(o.s_), p_(o.p_), invalid_(o.invalid_) {}

        KmerIterator& operator++();
        KmerIterator operator++(int);
        void raise(Kmer& km, Kmer& rep);

        bool operator==(const KmerIterator& o);
        bool operator!=(const KmerIterator& o) { return !this->operator==(o);}

        std::pair<Kmer, int>& operator*();
        std::pair<Kmer, int> *operator->();

    private:

        void find_next(size_t i, size_t j, bool last_valid);

        const char *s_;
        std::pair<Kmer, int> p_;
        bool invalid_;
};

template<class HF>
class KmerHashIterator {

    public:

        KmerHashIterator(const char* _s, const int _length, const int _k) : s(_s), n(_length), k(_k), hf(HF(_k)), p_(0,-1), invalid(true) {

            if ((_s != NULL) || (n >= k)) {

                invalid = false;
                operator++();
            }
        }

        KmerHashIterator() : s(NULL), n(0), k(0), hf(HF(0)), p_(0,-1), invalid(true) {}

        KmerHashIterator(const KmerHashIterator& o) : s(o.s), n(o.n), k(o.k), hf(o.hf), p_(o.p_), invalid(o.invalid) {}

        bool operator==(const KmerHashIterator& o) {

            if (invalid || o.invalid) return invalid && o.invalid;
            return s==o.s && n==o.n && k==o.k && p_.first==o.p_.first && p_.second==o.p_.second;
        }

        bool operator!=(const KmerHashIterator& o) { return !this->operator==(o); }

        KmerHashIterator& operator++() {

            if (invalid) return *this;

            ++(p_.second); // advance to next k-mer

            if (p_.second >= n - k + 1 || s[p_.second + k - 1] == '\0') { // out of bounds

                invalid = true;
                p_ = std::make_pair(0,-1);
                return *this;
            }

            char c, c_twin;

            int j = p_.second + k - 1;

            if (p_.second == 0){

                while (j >= p_.second) {

                    c = s[j] & 0xDF; // mask lowercase bit

                    if ((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) j--;
                    else {

                        p_.second += j - p_.second + 1;

                        if (p_.second >= n - k + 1 || s[p_.second + k - 1] == '\0') { // out of bounds

                            invalid = true;
                            p_ = std::make_pair(0,-1);
                            return *this;
                        }

                        j = p_.second + k - 1;
                    }
                }

                hf.init(&s[p_.second]);
            }
            else {

                c = s[j] & 0xDF; // mask lowercase bit

                if ((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) hf.update(s[j-k], s[j]);
                else {

                    p_.second += k;

                    if (p_.second >= n - k + 1 || s[p_.second + k - 1] == '\0') { // out of bounds

                        invalid = true;
                        p_ = std::make_pair(0,-1);
                        return *this;
                    }

                    j = p_.second + k - 1;

                    while (j >= p_.second) {

                        c = s[j] & 0xDF; // mask lowercase bit

                        if ((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) j--;
                        else {

                            p_.second += j - p_.second + 1;

                            if (p_.second >= n - k + 1 || s[p_.second + k - 1] == '\0') { // out of bounds

                                invalid = true;
                                p_ = std::make_pair(0,-1);
                                return *this;
                            }

                            j = p_.second + k - 1;
                        }
                    }

                    hf.init(&s[p_.second]);
                }
            }

            p_.first = hf.hash();

            return *this;
        }

        KmerHashIterator operator++(int) {

            KmerHashIterator tmp(*this);
            operator++();

            return tmp;
        }

        //Move iterator to next VALID position >= p_.second + length
        KmerHashIterator& operator+=(const int length){

            size_t next_pos = p_.second + length;

            while (!invalid && p_.second < next_pos) operator++();

            return *this;
        }

        std::pair<uint64_t, int>& operator*() { return p_; }

        std::pair<uint64_t, int>* operator->() { return &(operator*()); }

        const char *s; // K-mers are from a sequence s
        int n; // Length of sequence s
        int k; // Length of k-mers
        HF hf; // Rolling hash function for k-mers of s
        std::pair<uint64_t, int> p_; // <hash, position> current k-mer
        bool invalid; // If sequence is invalid (iterating on k-mers out of bounds, etc.)
};

#endif // BFG_KMER_ITERATOR_HPP
