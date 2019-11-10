#ifndef BFG_MINIMIZER_IDX_HPP
#define BFG_MINIMIZER_IDX_HPP

#include <utility>
#include <string>
#include <iterator>
#include <algorithm>

#include "Kmer.hpp"
#include "Lock.hpp"
#include "TinyVector.hpp"

class MinimizerIndex {

    template<bool is_const = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, packed_tiny_vector> {

        public:

            typedef typename std::conditional<is_const, const MinimizerIndex*, MinimizerIndex*>::type MI_ptr_t;
            typedef typename std::conditional<is_const, const packed_tiny_vector&, packed_tiny_vector&>::type MI_tinyv_ref_t;
            typedef typename std::conditional<is_const, const packed_tiny_vector*, packed_tiny_vector*>::type MI_tinyv_ptr_t;
            typedef typename std::conditional<is_const, const uint8_t&, uint8_t&>::type MI_tinyv_sz_ref_t;
            typedef typename std::conditional<is_const, const uint8_t*, uint8_t*>::type MI_tinyv_sz_ptr_t;

            iterator_() : ht(nullptr), h(0xffffffffffffffffULL) {}
            iterator_(MI_ptr_t ht_) : ht(ht_), h(ht_->size_) {}
            iterator_(MI_ptr_t ht_, size_t h_) :  ht(ht_), h(h_) {}
            iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h) {}
            iterator_& operator=(const iterator_& o) { ht=o.ht; h=o.h; return *this; }

            const Minimizer getKey() const {

                /*ht->lck_edit_keys.acquire();

                const Minimizer minz(ht->table_keys[h]);

                ht->lck_edit_keys.release();

                return minz;*/

                return ht->table_keys[h];
            }

            inline size_t getHash() const {

                return h;
            }

            inline MI_tinyv_sz_ref_t getVectorSize() const {

                return ht->table_tinyv_sz[h];
            }

            inline MI_tinyv_ref_t getVector() const {

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

                if (h == ht->size_) return *this;

                ++h;

                for (; h < ht->size_; ++h) {

                    if ((ht->table_keys[h] != ht->empty_key) && (ht->table_keys[h] != ht->deleted_key)) break;
                }

                return *this;
            }

            bool operator==(const iterator_ &o) const { return (ht == o.ht) && (h == o.h); }
            bool operator!=(const iterator_ &o) const { return (ht != o.ht) || (h != o.h); }

            friend class iterator_<true>;

        //private:

            MI_ptr_t ht;
            size_t h;
    };

    public:

        template<bool is_const> friend class iterator_;

        typedef iterator_<true> const_iterator;
        typedef iterator_<false> iterator;

        MinimizerIndex() :  table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                            size_(0), pop(0), num_empty(0), lck_min_block_sz(64) {

            empty_key.set_empty();
            deleted_key.set_deleted();

            init_tables(1024);
        }

        MinimizerIndex(const size_t sz) :   table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                            size_(0), pop(0), num_empty(0), lck_min_block_sz(64) {

            empty_key.set_empty();
            deleted_key.set_deleted();

            init_tables(std::max(static_cast<size_t>(1.2 * sz), static_cast<size_t>(lck_min_block_sz)));
        }

        MinimizerIndex(const MinimizerIndex& o) :   size_(o.size_), pop(o.pop), num_empty(o.num_empty),
                                                    empty_key(o.empty_key), deleted_key(o.deleted_key) {

            table_keys = new Minimizer[size_];
            table_tinyv = new packed_tiny_vector[size_];
            table_tinyv_sz = new uint8_t[size_];

            lck_min = vector<SpinLock>(o.lck_min.size());
            lck_min_block_sz = o.lck_min_block_sz;

            std::copy(o.table_keys, o.table_keys + size_, table_keys);

            for (size_t i = 0; i < size_; ++i){

                table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;
                table_tinyv[i].copy(table_tinyv_sz[i], o.table_tinyv[i], o.table_tinyv_sz[i]);
            }
        }

        MinimizerIndex(MinimizerIndex&& o){

            size_ = o.size_;
            pop = o.pop;
            num_empty = o.num_empty;

            empty_key = o.empty_key;
            deleted_key = o.deleted_key;

            table_keys = o.table_keys;
            table_tinyv = o.table_tinyv;
            table_tinyv_sz = o.table_tinyv_sz;

            lck_min = vector<SpinLock>(o.lck_min.size());
            lck_min_block_sz = o.lck_min_block_sz;

            o.table_keys = nullptr;
            o.table_tinyv = nullptr;
            o.table_tinyv_sz = nullptr;

            o.clear();
        }

        MinimizerIndex& operator=(const MinimizerIndex& o) {

            if (this != &o) {

                clear();

                size_ = o.size_;
                pop = o.pop;
                num_empty = o.num_empty;

                empty_key = o.empty_key;
                deleted_key = o.deleted_key;

                table_keys = new Minimizer[size_];
                table_tinyv = new packed_tiny_vector[size_];
                table_tinyv_sz = new uint8_t[size_];

                lck_min = vector<SpinLock>(o.lck_min.size());
                lck_min_block_sz = o.lck_min_block_sz;

                std::copy(o.table_keys, o.table_keys + size_, table_keys);

                for (size_t i = 0; i < size_; ++i){

                    table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;
                    table_tinyv[i].copy(table_tinyv_sz[i], o.table_tinyv[i], o.table_tinyv_sz[i]);
                }
            }

            return *this;
        }

        MinimizerIndex& operator=(MinimizerIndex&& o){

            if (this != &o) {

                clear();

                size_ = o.size_;
                pop = o.pop;
                num_empty = o.num_empty;

                empty_key = o.empty_key;
                deleted_key = o.deleted_key;

                table_keys = o.table_keys;
                table_tinyv = o.table_tinyv;
                table_tinyv_sz = o.table_tinyv_sz;

                lck_min = vector<SpinLock>(o.lck_min.size());
                lck_min_block_sz = o.lck_min_block_sz;

                o.table_keys = nullptr;
                o.table_tinyv = nullptr;
                o.table_tinyv_sz = nullptr;

                o.clear();
            }

            return *this;
        }

        ~MinimizerIndex() {

            clear();
        }

        size_t size() const {

            lck_edit_table.acquire_reader();

            const size_t pop_ret = pop;

            lck_edit_table.release_reader();

            return pop_ret;
        }

        bool empty() const {

            lck_edit_table.acquire_reader();

            const bool pop_ret = (pop == 0);

            lck_edit_table.release_reader();

            return pop_ret;
        }

        void clear() {

            if (table_tinyv != nullptr){

                for (size_t i = 0; i < size_; ++i) table_tinyv[i].destruct(table_tinyv_sz[i]);
            }

            clear_tables();

            lck_min.clear();
            lck_edit_table.release_all();

            lck_min_block_sz = 0;
        }

        iterator find(const Minimizer& key) {

            const size_t end_table = size_-1;
            size_t h = key.hash() & end_table;
            const size_t end_h = (h-1) & end_table;

            for (; h != end_h; h = (h+1) & end_table) {

                if ((table_keys[h] == empty_key) || (table_keys[h] == key)) break;
            }

            if ((h != end_h) && (table_keys[h] == key)) return iterator(this, h);

            return iterator(this);
        }

        const_iterator find(const Minimizer& key) const {

            const size_t end_table = size_-1;
            size_t h = key.hash() & end_table;
            const size_t end_h = (h-1) & end_table;

            for (; h != end_h; h = (h+1) & end_table) {

                if ((table_keys[h] == empty_key) || (table_keys[h] == key)) break;
            }

            if ((h != end_h) && (table_keys[h] == key)) return const_iterator(this, h);

            return const_iterator(this);
        }

        iterator find(const size_t h) {

            if ((h < size_) && (table_keys[h] != empty_key) && (table_keys[h] != deleted_key)) return iterator(this, h);

            return iterator(this);
        }

        const_iterator find(const size_t h) const {

            if ((h < size_) && (table_keys[h] != empty_key) && (table_keys[h] != deleted_key)) return const_iterator(this, h);

            return const_iterator(this);
        }

        iterator erase(const_iterator it) {

            if (it == end()) return end();

            table_keys[it.h] = deleted_key;
            table_tinyv[it.h].destruct(table_tinyv_sz[it.h]);
            table_tinyv_sz[it.h] = packed_tiny_vector::FLAG_EMPTY;

            --pop;

            return ++iterator(this, it.h); // return pointer to next element
        }

        size_t erase(const Minimizer& minz) {

            const size_t end_table = size_-1;
            size_t h = minz.hash() & end_table;
            const size_t end_h = (h-1) & end_table;
            const size_t oldpop = pop;

            for (; h != end_h; h = (h+1) & end_table) {

                if ((table_keys[h] == empty_key) || (table_keys[h] == minz)) break;
            }

            if ((h != end_h) && (table_keys[h] == minz)){

                table_keys[h] = deleted_key;
                table_tinyv[h].destruct(table_tinyv_sz[h]);
                table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;

                --pop;
            }

            return oldpop - pop;
        }

        pair<iterator, bool> insert(const Minimizer& key, const packed_tiny_vector& v, const uint8_t& flag) {

            if ((5 * num_empty) < size_) reserve(2 * size_); // if more than 80% full, resize

            bool is_deleted = false;

            const size_t end_table = size_-1;

            for (size_t h = key.hash() & (size_-1), h_tmp;; h = (h+1) & end_table) {

                if (table_keys[h] == empty_key) {

                    is_deleted ? h = h_tmp : --num_empty;

                    table_keys[h] = key;
                    table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;

                    table_tinyv[h].copy(table_tinyv_sz[h], v, flag);

                    ++pop;

                    return {iterator(this, h), true};
                }
                else if (table_keys[h] == key){

                    return {iterator(this, h), false};
                }
                else if (!is_deleted && (table_keys[h] == deleted_key)) {

                    is_deleted = true;
                    h_tmp = h;
                }
            }
        }

        //-------------------------------------
        void init_threads() const {

            lck_edit_table.acquire_reader();

            lck_min = vector<SpinLock>((size_ + lck_min_block_sz - 1) / lck_min_block_sz);

            lck_edit_table.release_reader();
        }

        void release_threads() const {

            lck_min.clear();
        }

        iterator find_p(const Minimizer& key) {

            size_t h = key.hash();

            lck_edit_table.acquire_reader();

            const size_t end_table = size_-1;

            h &= end_table;

            size_t id_block = h / lck_min_block_sz;

            const size_t end_h = (h-1) & end_table;

            lck_min[id_block].acquire();

            for (; h != end_h; h = (h+1) & end_table) {

                if ((h / lck_min_block_sz) != id_block){

                    lck_min[id_block].release();
                    id_block = h / lck_min_block_sz;
                    lck_min[id_block].acquire();
                }

                if (table_keys[h] == empty_key){

                    lck_min[id_block].release();
                    lck_edit_table.release_reader();

                    iterator(this);
                }
                else if (table_keys[h] == key){

                    lck_min[id_block].release();
                    lck_edit_table.release_reader();

                    return iterator(this, h);
                }
            }

            lck_min[id_block].release();
            lck_edit_table.release_reader();

            return iterator(this);
        }

        const_iterator find_p(const Minimizer& key) const {

            size_t h = key.hash();

            lck_edit_table.acquire_reader();

            const size_t end_table = size_-1;

            h &= end_table;

            size_t id_block = h / lck_min_block_sz;

            const size_t end_h = (h-1) & end_table;

            lck_min[id_block].acquire();

            for (; h != end_h; h = (h+1) & end_table) {

                if ((h / lck_min_block_sz) != id_block){

                    lck_min[id_block].release();
                    id_block = h / lck_min_block_sz;
                    lck_min[id_block].acquire();
                }

                if (table_keys[h] == empty_key){

                    lck_min[id_block].release();
                    lck_edit_table.release_reader();

                    const_iterator(this);
                }
                else if (table_keys[h] == key){

                    lck_min[id_block].release();
                    lck_edit_table.release_reader();

                    return const_iterator(this, h);
                }
            }

            lck_min[id_block].release();
            lck_edit_table.release_reader();

            return const_iterator(this);
        }

        iterator find_p(const size_t h) {

            lck_edit_table.acquire_reader();

            if (h < size_){

                const size_t id_block = h / lck_min_block_sz;

                lck_min[id_block].acquire();

                if ((table_keys[h] != empty_key) && (table_keys[h] != deleted_key)) {

                    lck_min[id_block].release();
                    lck_edit_table.release_reader();

                    return iterator(this, h);
                }

                lck_min[id_block].release();
            }

            lck_edit_table.release_reader();

            return iterator(this);
        }

        const_iterator find_p(const size_t h) const {

            lck_edit_table.acquire_reader();

            if (h < size_){

                const size_t id_block = h / lck_min_block_sz;

                lck_min[id_block].acquire();

                if ((table_keys[h] != empty_key) && (table_keys[h] != deleted_key)) {

                    lck_min[id_block].release();
                    lck_edit_table.release_reader();

                    return const_iterator(this, h);
                }

                lck_min[id_block].release();
            }

            lck_edit_table.release_reader();

            return const_iterator(this);
        }

        void release_p(const_iterator it) const {

            if (it != end()){

                lck_min[it.h / lck_min_block_sz].release();
                lck_edit_table.release_reader();
            }
        }

        void release_p(iterator it) {

            if (it != end()){

                lck_min[it.h / lck_min_block_sz].release();
                lck_edit_table.release_reader();
            }
        }

        iterator erase_p(const_iterator it) {

            if (it == end()) return end();

            const size_t id_block = it.h / lck_min_block_sz;

            lck_edit_table.acquire_reader();
            lck_min[id_block].acquire();

            table_keys[it.h] = deleted_key;
            table_tinyv[it.h].destruct(table_tinyv_sz[it.h]);
            table_tinyv_sz[it.h] = packed_tiny_vector::FLAG_EMPTY;

            lck_min[id_block].release();

            --pop;

            lck_edit_table.release_reader();

            return ++iterator(this, it.h); // return pointer to next element
        }

        size_t erase_p(const Minimizer& minz) {

            size_t h = minz.hash();

            lck_edit_table.acquire_reader();

            const size_t end_table = size_-1;

            h &= end_table;

            size_t id_block = h / lck_min_block_sz;
            size_t l_pop = pop;

            const size_t end_h = (h-1) & end_table;

            lck_min[id_block].acquire();

            for (; h != end_h; h = (h+1) & end_table) {

                if ((h / lck_min_block_sz) != id_block){

                    lck_min[id_block].release();
                    id_block = h / lck_min_block_sz;
                    lck_min[id_block].acquire();
                }

                if (table_keys[h] == empty_key){

                    lck_min[id_block].release();
                    lck_edit_table.release_reader();

                    return 0;
                }
                else if (table_keys[h] == minz) break;
            }

            if ((h != end_h) && (table_keys[h] == minz)){

                table_keys[h] = deleted_key;
                table_tinyv[h].destruct(table_tinyv_sz[h]);
                table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;

                lck_min[id_block].release();

                --pop;
            }

            l_pop -= pop;

            lck_edit_table.release_reader();

            return l_pop;
        }

        pair<iterator, bool> insert_p(const Minimizer& key, const packed_tiny_vector& v, const uint8_t& flag) {

            size_t h = key.hash();

            lck_edit_table.acquire_reader();

            if ((5 * num_empty) < size_){

                lck_edit_table.release_reader();
                lck_edit_table.acquire_writer();

                reserve(2 * size_); // if more than 80% full, resize

                lck_edit_table.release_writer_acquire_reader();
            }

            bool is_deleted = false;

            const size_t end_table = size_-1;

            h &= end_table;

            size_t id_block = h / lck_min_block_sz;

            lck_min[id_block].acquire();

            for (size_t h_tmp;; h = (h+1) & end_table) {

                if ((h / lck_min_block_sz) != id_block){

                    lck_min[id_block].release();
                    id_block = h / lck_min_block_sz;
                    lck_min[id_block].acquire();
                }

                if (table_keys[h] == empty_key) {

                    is_deleted ? h = h_tmp : --num_empty;

                    table_keys[h] = key;
                    table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;

                    table_tinyv[h].copy(table_tinyv_sz[h], v, flag);

                    ++pop;

                    return {iterator(this, h), true};
                }
                else if (table_keys[h] == key){

                    return {iterator(this, h), false};
                }
                else if (!is_deleted && (table_keys[h] == deleted_key)) {

                    is_deleted = true;
                    h_tmp = h;
                }
            }

            lck_min[id_block].release(); // Just for safety
            lck_edit_table.release_reader(); // Just for safety
        }

        iterator begin() {

            iterator it(this);
            it.operator++();
            return it;
        }

        const_iterator begin() const {

            const_iterator it(this);
            it.operator++();
            return it;
        }

        iterator end() { return iterator(this); }

        const_iterator end() const { return const_iterator(this); }

    private:

        void clear_tables() {

            if (table_keys != nullptr) {

                delete[] table_keys;
                table_keys = nullptr;
            }

            if (table_tinyv != nullptr) {

                delete[] table_tinyv;
                table_tinyv = nullptr;
            }

            if (table_tinyv_sz != nullptr) {

                delete[] table_tinyv_sz;
                table_tinyv_sz = nullptr;
            }

            size_ = 0;
            pop  = 0;
            num_empty = 0;
        }

        void init_tables(const size_t sz) {

            clear_tables();

            pop = 0;
            size_ = rndup(sz);
            num_empty = size_;

            table_keys = new Minimizer[size_];
            table_tinyv = new packed_tiny_vector[size_];
            table_tinyv_sz = new uint8_t[size_];

            std::fill(table_keys, table_keys + size_, empty_key);

            memset(table_tinyv_sz, packed_tiny_vector::FLAG_EMPTY, size_ * sizeof(uint8_t));
        }

        void reserve(const size_t sz) {

            if (sz <= size_) return;

            Minimizer* old_table_keys = table_keys;
            packed_tiny_vector* old_table_tinyv = table_tinyv;
            uint8_t* old_table_tinyv_sz = table_tinyv_sz;

            const size_t old_size_ = size_;

            size_ = rndup(sz);
            pop = 0;
            num_empty = size_;

            table_keys = new Minimizer[size_];
            table_tinyv = new packed_tiny_vector[size_];
            table_tinyv_sz = new uint8_t[size_];

            std::fill(table_keys, table_keys + size_, empty_key);

            memset(table_tinyv_sz, packed_tiny_vector::FLAG_EMPTY, size_ * sizeof(uint8_t));

            for (size_t i = 0; i < old_size_; ++i) {

                if (old_table_keys[i] != empty_key && old_table_keys[i] != deleted_key){

                    insert(old_table_keys[i], old_table_tinyv[i], old_table_tinyv_sz[i]);

                    old_table_tinyv[i].destruct(old_table_tinyv_sz[i]);
                }
            }

            delete[] old_table_keys;
            delete[] old_table_tinyv;
            delete[] old_table_tinyv_sz;
        }

        size_t size_, pop, num_empty;

        Minimizer* table_keys;
        packed_tiny_vector* table_tinyv;
        uint8_t* table_tinyv_sz;

        Minimizer empty_key;
        Minimizer deleted_key;

        mutable vector<SpinLock> lck_min;
        mutable SpinLockRW lck_edit_table;

        int64_t lck_min_block_sz;
};

#endif
