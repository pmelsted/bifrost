#ifndef BIFROST_MINIMIZER_IDX_HPP
#define BIFROST_MINIMIZER_IDX_HPP

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iterator>
#include <mutex>
#include <thread>
#include <string>
#include <utility>

#include "fastmod.h"

#include "Kmer.hpp"
#include "Lock.hpp"
#include "TinyVector.hpp"

#define BIFROST_MI_MAX_OCCUPANCY 0.95
#define BIFROST_MI_INIT_SZ 128

class MinimizerIndex {

    template<bool is_const = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, packed_tiny_vector> {

        public:

            typedef typename std::conditional<is_const, const MinimizerIndex*, MinimizerIndex*>::type MI_ptr_t;
            typedef typename std::conditional<is_const, const packed_tiny_vector&, packed_tiny_vector&>::type MI_tinyv_ref_t;
            typedef typename std::conditional<is_const, const packed_tiny_vector*, packed_tiny_vector*>::type MI_tinyv_ptr_t;
            typedef typename std::conditional<is_const, const uint8_t&, uint8_t&>::type MI_tinyv_sz_ref_t;
            typedef typename std::conditional<is_const, const uint8_t*, uint8_t*>::type MI_tinyv_sz_ptr_t;

            iterator_() : ht(nullptr), h(0xffffffffffffffffULL), psl(0xffffffffffffffffULL) {}
            iterator_(MI_ptr_t ht_) : ht(ht_), h(0xffffffffffffffffULL), psl(0xffffffffffffffffULL) {}
            iterator_(MI_ptr_t ht_, size_t h_) :  ht(ht_), h(h_), psl(0xffffffffffffffffULL) {}
            iterator_(MI_ptr_t ht_, size_t h_, size_t psl_) :  ht(ht_), h(h_), psl(psl_) {}
            iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h), psl(o.psl) {}

            iterator_& operator=(const iterator_& o) {

                if (this != &o) {

                    ht=o.ht;
                    h=o.h;
                    psl=o.psl;
                }

                return *this;
            }

            BFG_INLINE Minimizer getKey() const {

                return ht->table_keys[h];
            }

            BFG_INLINE size_t getHash() const {

                return h;
            }

            BFG_INLINE size_t getPSL() const {

                return psl;
            }

            BFG_INLINE MI_tinyv_sz_ref_t getVectorSize() const {

                return ht->table_tinyv_sz[h];
            }

            BFG_INLINE MI_tinyv_ref_t getVector() const {

                return ht->table_tinyv[h];
            }

            MI_tinyv_ref_t operator*() const {

                return ht->table_tinyv[h];
            }

            MI_tinyv_ptr_t operator->() const {

                return &(ht->table_tinyv[h]);
            }

            iterator_ operator++(int) {

                const iterator_ tmp(*this);
                operator++();
                return tmp;
            }

            iterator_& operator++() {

                h += static_cast<size_t>(h < ht->size_);

                while ((h < ht->size_) && ht->table_keys[h].isEmpty()) ++h;

                h |= static_cast<size_t>(h < ht->size_) - 1;
                psl = 0xffffffffffffffffULL;

                return *this;
            }

            BFG_INLINE bool operator==(const iterator_ &o) const {

                return (ht == o.ht) && (h == o.h);
            }

            BFG_INLINE bool operator!=(const iterator_ &o) const {

                return (ht != o.ht) || (h != o.h);
            }

            void get_to_first() {

                h = 0xffffffffffffffffULL;
                psl = 0xffffffffffffffffULL;

                if ((ht != nullptr) && (ht->size_ != 0)) {

                    h = 0;

                    while ((h < ht->size_) && ht->table_keys[h].isEmpty()) ++h;

                    h |= static_cast<size_t>(h < ht->size_) - 1;
                }
            }

            friend class iterator_<true>;

        //private:

            MI_ptr_t ht;

            size_t h;
            size_t psl;
    };

    public:

        template<bool is_const> friend class iterator_;

        typedef iterator_<true> const_iterator;
        typedef iterator_<false> iterator;

        MinimizerIndex();
        MinimizerIndex(const size_t sz, const double max_ratio_occupancy = BIFROST_MI_MAX_OCCUPANCY);

        MinimizerIndex(const MinimizerIndex& o);
        MinimizerIndex(MinimizerIndex&& o);

        MinimizerIndex& operator=(const MinimizerIndex& o);
        MinimizerIndex& operator=(MinimizerIndex&& o);

        ~MinimizerIndex();

        BFG_INLINE size_t size() const {

            return pop;
        }

        BFG_INLINE size_t capacity() const {

            return size_;
        }

        BFG_INLINE bool empty() const {

            return (pop == 0);
        }

        void clear();

        iterator find(const Minimizer& key);
        const_iterator find(const Minimizer& key) const;

        iterator find(const size_t h);
        const_iterator find(const size_t h) const;

        size_t erase(const_iterator it);

        BFG_INLINE size_t erase(const Minimizer& minz) {

            const const_iterator it = find(minz);

            return erase(it);
        }

        pair<iterator, bool> insert(const Minimizer& key, const packed_tiny_vector& v, const uint8_t& flag);

        iterator begin();
        const_iterator begin() const;

        iterator end();
        const_iterator end() const;

        void recomputeMaxPSL(const size_t nb_threads = 1);
        void recomputeMaxStdPSL(const size_t nb_threads = 1);

        BFG_INLINE size_t get_mean_psl() const {

            return ceil(static_cast<double>(sum_psl) / static_cast<double>(pop + 1)); // Slightly biased computation but avoids to check for (psl == 0). Fine since we just need an approximate result.
        }

        BFG_INLINE size_t get_std_psl() const {

            return std_psl;
        }

        BFG_INLINE size_t get_max_psl() const {

            return max_psl;
        }

    private:

        void clear_tables();
        void init_tables(const size_t sz);
        void reserve(const size_t sz);
        void swap(const size_t i, const size_t j);

        double max_ratio_occupancy;

        __uint128_t M_u64;

        size_t size_, pop;
        size_t max_psl, sum_psl, std_psl;

        Minimizer* table_keys;

        packed_tiny_vector* table_tinyv;
        uint8_t* table_tinyv_sz;

        uint64_t* table_outliers_psl;
};

#endif
