#ifndef KALLISTO_KMERHASHTABLE_H
#define KALLISTO_KMERHASHTABLE_H

#include <utility>
#include <string>
#include <iterator>

#include "Kmer.hpp"

template<typename T, typename Hash = KmerHash>
struct KmerHashTable {

    using value_type = std::pair<Kmer, T>;
    using key_type = Kmer;
    using mapped_type = T;

    Hash hasher;
    value_type *table;
    size_t size_, pop, num_empty;
    value_type empty_val;
    value_type deleted;

// ---- iterator ----

    template<bool is_const_iterator = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, value_type> {

        public:

            typedef typename std::conditional<is_const_iterator, const KmerHashTable *, KmerHashTable *>::type DataStructurePointerType;
            typedef typename std::conditional<is_const_iterator, const value_type&, value_type&>::type ValueReferenceType;
            typedef typename std::conditional<is_const_iterator, const value_type *, value_type *>::type ValuePointerType;


            DataStructurePointerType ht;
            size_t h;

            iterator_() : ht(nullptr), h(0) {}
            iterator_(DataStructurePointerType ht_) : ht(ht_), h(ht_->size_) {}
            iterator_(DataStructurePointerType ht_, size_t h_) :  ht(ht_), h(h_) {}
            iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h) {}
            iterator_& operator=(const iterator_& o) {ht=o.ht; h=o.h; return *this;}

            ValueReferenceType operator*() const {return ht->table[h];}
            ValuePointerType operator->() const {return &(ht->table[h]);}

            size_t getHash() const { return h; }

            void find_first() {

                h = 0;

                if (ht->table != nullptr && ht->size_>0) {

                    Kmer& km = ht->table[h].first;
                    if (km == ht->empty_val.first || km == ht->deleted.first) operator++();
                }
            }

            iterator_ operator++(int) {

                const iterator_ old(*this);
                ++(*this);
                return old;
            }

            iterator_& operator++() {

                if (h == ht->size_) return *this;

                ++h;

                for (; h < ht->size_; ++h) {

                    Kmer& km = ht->table[h].first;

                    if (km != ht->empty_val.first && km != ht->deleted.first) break;
                }

                return *this;
            }

            bool operator==(const iterator_ &o) const {return (ht->table == o.ht->table) && (h == o.h);}
            bool operator!=(const iterator_ &o) const {return !(this->operator==(o));}

            friend class iterator_<true>;
    };

    typedef iterator_<true> const_iterator;
    typedef iterator_<false> iterator;


  // --- hash table
    KmerHashTable(const Hash& h = Hash() ) : hasher(h), table(nullptr), size_(0), pop(0), num_empty(0) {

        empty_val.first.set_empty();
        deleted.first.set_deleted();
        init_table(1024);
    }

    KmerHashTable(size_t sz, const Hash& h = Hash() ) : hasher(h), table(nullptr), size_(0), pop(0), num_empty(0) {
        empty_val.first.set_empty();
        deleted.first.set_deleted();
        init_table((size_t) (1.2*sz));
    }

    KmerHashTable(KmerHashTable&& o){

        hasher = o.hasher;
        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;
        table = o.table;
        empty_val = o.empty_val;
        deleted = o.deleted;

        o.table = nullptr;

        o.clear_table();
    }

    KmerHashTable& operator=(KmerHashTable&& o){

        if (this != &o) {

            clear_table();

            hasher = o.hasher;
            size_ = o.size_;
            pop = o.pop;
            num_empty = o.num_empty;
            table = o.table;
            empty_val = o.empty_val;
            deleted = o.deleted;

            o.table = nullptr;

            o.clear_table();
        }

        return *this;
    }

    ~KmerHashTable() { clear_table(); }

    void clear_table() {

        if (table != nullptr) {

            delete[] table;
            //free(table);
            table = nullptr;
        }

        size_ = 0;
        pop  = 0;
        num_empty = 0;
    }

    size_t size() const { return pop; }

    bool empty() const { return pop == 0; }

    void clear() {

        std::fill(table, table+size_, empty_val);

        pop = 0;
        num_empty = size_;
    }

    void init_table(size_t sz) {

        clear_table();

        size_ = rndup(sz);

        table = new value_type[size_];
        //table = (value_type*) malloc(size_ * sizeof(value_type));

        clear();
    }

    iterator find(const Kmer& key) {

        size_t h = hasher(key) & (size_-1);
        size_t end_h = (h == 0) ? (size_-1) : h-1;

        for (;; h =  (h+1!=size_ ? h+1 : 0)) {

            if (table[h].first == empty_val.first) return iterator(this); // empty slot, not in table
            else if (table[h].first == key) return iterator(this, h); // same key, found
            // if it is deleted, we still have to continue
            if (h==end_h) return iterator(this); // we've gone throught the table, quit
        }
    }

    const_iterator find(const Kmer& key) const {

        size_t h = hasher(key) & (size_-1);
        size_t end_h = (h == 0) ? (size_-1) : h-1;

        for (;; h =  (h+1!=size_ ? h+1 : 0)) {

            if (table[h].first == empty_val.first) return const_iterator(this); // empty slot, not in table
            else if (table[h].first == key) return const_iterator(this, h); // same key, found

            if (h==end_h) return const_iterator(this);
        }
    }

    iterator find(const size_t h) {

        if ((h < size_) && (table[h].first != empty_val.first) && (table[h].first != deleted.first))
            return iterator(this, h);

        return iterator(this);
    }

    const_iterator find(const size_t h) const {

        if ((h < size_) && (table[h].first != empty_val.first) && (table[h].first != deleted.first))
            return const_iterator(this, h);

        return const_iterator(this);
    }

    iterator erase(const_iterator pos) {

        if (pos == this->end()) return this->end();

        size_t h = pos.h;

        table[h] = deleted;
        --pop;

        return ++iterator(this, h); // return pointer to next element
    }

    iterator erase(const size_t h) {

        if (h >= size_) return this->end();

        table[h] = deleted;
        --pop;

        return ++iterator(this, h); // return pointer to next element
    }

    size_t erase(const Kmer& km) {

        const_iterator pos = find(km);
        size_t oldpop = pop;

        if (pos != this->end()) erase(pos);

        return oldpop-pop;
    }

    std::pair<iterator,bool> insert(const value_type& val) {

        if ((5*num_empty) < size_) reserve(2*size_); // if more than 80% full, resize

        bool is_deleted = false;

        for (size_t h = hasher(val.first) & (size_-1), h_tmp;; h = (h+1 != size_ ? h+1 : 0)) {

            if (table[h].first == empty_val.first) {

                if (!is_deleted) num_empty--;
                else h = h_tmp;

                table[h] = val;
                ++pop;

                return {iterator(this, h), true};
            }
            else if (table[h].first == val.first) return {iterator(this, h), false};
            else if (!is_deleted && (table[h].first == deleted.first)) {
                is_deleted = true;
                h_tmp = h;
            }
        }
    }

    void reserve(size_t sz) {

        if (sz <= size_) return;

        value_type *old_table = table;
        size_t old_size_ = size_;

        size_ = rndup(sz);
        pop = 0;
        num_empty = size_;

        table = new value_type[size_];

        std::fill(table, table+size_, empty_val);

        for (size_t i = 0; i < old_size_; i++) {

            if (old_table[i].first != empty_val.first && old_table[i].first != deleted.first) insert(old_table[i]);
        }

        delete[] old_table;
        old_table = nullptr;

        /*if (sz <= size_) return;

        const size_t prev_size_ = size_;

        size_ = rndup(sz);
        pop = 0;
        num_empty = size_;

        table = (value_type*) realloc(table, size_ * sizeof(value_type));

        std::sort(table, table + prev_size_, sortKmerHashTable(*this));

        value_type* table_empty = table;

        for (; table_empty < table + prev_size_; table_empty++) {

            if (((*table_empty).first == empty_val.first) || ((*table_empty).first == deleted.first)) break;
        }

        std::fill(table_empty, table + size_, empty_val);

        std::vector<value_type> v;

        for (int64_t i = table_empty - table - 1; i >= 0; i--) {

            size_t h = hasher(table[i].first) & (size_-1);

            if (h > i){

                for ( ; h < size_; h++) {

                    if (table[h].first == empty_val.first) {

                        num_empty--;
                        pop++;

                        std::swap(table[h], table[i]);

                        break;
                    }
                }

                if (h == size_){

                    v.push_back(table[i]);
                    table[i] = empty_val;
                }
            }
        }

        for (auto& vt : v) insert(vt);*/
    }

    size_t rndup(size_t v) {
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

    iterator begin() {

        iterator it(this);
        it.find_first();
        return it;
    }

    const_iterator begin() const {

        const_iterator it(this);
        it.find_first();
        return it;
    }

    iterator end() { return iterator(this); }

    const_iterator end() const { return const_iterator(this); }

    /*private:

        struct sortKmerHashTable {

            sortKmerHashTable(const KmerHashTable& kht_) : kht(kht_) {}

            bool operator() (const value_type& a, const value_type& b) const {

                const size_t h_a = (a.first == kht.empty_val.first) || (a.first == kht.deleted.first) ? 0xffffffffffffffff : kht.hasher(a.first) & (kht.size_-1);
                const size_t h_b = (b.first == kht.empty_val.first) || (b.first == kht.deleted.first) ? 0xffffffffffffffff : kht.hasher(b.first) & (kht.size_-1);

                return (h_a < h_b);
            }

            const KmerHashTable& kht;
        };*/
};

template<typename T, typename Hash = MinimizerHash>
struct MinimizerHashTable {

    using value_type = std::pair<Minimizer, T>;
    using key_type = Minimizer;
    using mapped_type = T;

    Hash hasher;

    size_t size_, pop, num_empty;

    value_type* table;

    value_type empty_val;
    value_type deleted;

// ---- iterator ----
    template<bool is_const_iterator = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, value_type> {

        public:

            typedef typename std::conditional<is_const_iterator, const MinimizerHashTable *, MinimizerHashTable *>::type DataStructurePointerType;
            typedef typename std::conditional<is_const_iterator, const value_type&, value_type&>::type ValueReferenceType;
            typedef typename std::conditional<is_const_iterator, const value_type *, value_type *>::type ValuePointerType;

            DataStructurePointerType ht;
            size_t h;

            iterator_() : ht(nullptr), h(0) {}
            iterator_(DataStructurePointerType ht_) : ht(ht_), h(ht_->size_) {}
            iterator_(DataStructurePointerType ht_, size_t h_) :  ht(ht_), h(h_) {}
            iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h) {}
            iterator_& operator=(const iterator_& o) {ht=o.ht; h=o.h; return *this;}

            ValueReferenceType operator*() const {return ht->table[h];}
            ValuePointerType operator->() const {return &(ht->table[h]);}

            size_t getHash() const { return h; }

            void find_first() {

                h = 0;

                if (ht->table != nullptr && ht->size_>0) {

                    Minimizer& minz = ht->table[h].first;

                    if (minz == ht->empty_val.first || minz == ht->deleted.first) operator++();
                }
            }

            iterator_ operator++(int) {

                const iterator_ old(*this);
                ++(*this);
                return old;
            }

            iterator_& operator++() {

                if (h == ht->size_) return *this;

                ++h;

                for (; h < ht->size_; ++h) {

                    Minimizer& minz = ht->table[h].first;

                    if (minz != ht->empty_val.first && minz != ht->deleted.first) break;
                }

                return *this;
            }

            bool operator==(const iterator_ &o) const {return (ht->table == o.ht->table) && (h == o.h);}
            bool operator!=(const iterator_ &o) const {return !(this->operator==(o));}
            friend class iterator_<true>;
        };

        typedef iterator_<true> const_iterator;
        typedef iterator_<false> iterator;

        // --- hash table
        MinimizerHashTable(const Hash& h = Hash() ) : hasher(h), table(nullptr), size_(0), pop(0), num_empty(0) {

            empty_val.first.set_empty();
            deleted.first.set_deleted();
            init_table(1024);
        }

        MinimizerHashTable(size_t sz, const Hash& h = Hash() ) : hasher(h), table(nullptr), size_(0), pop(0), num_empty(0) {

            empty_val.first.set_empty();
            deleted.first.set_deleted();
            init_table((size_t) (1.2*sz));
        }

        MinimizerHashTable(MinimizerHashTable&& o){

            hasher = o.hasher;
            size_ = o.size_;
            pop = o.pop;
            num_empty = o.num_empty;
            table = o.table;
            empty_val = o.empty_val;
            deleted = o.deleted;

            o.table = nullptr;

            o.clear_table();
        }

        MinimizerHashTable& operator=(MinimizerHashTable&& o){

            if (this != &o) {

                clear_table();

                hasher = o.hasher;
                size_ = o.size_;
                pop = o.pop;
                num_empty = o.num_empty;
                table = o.table;
                empty_val = o.empty_val;
                deleted = o.deleted;

                o.table = nullptr;

                o.clear_table();
            }

            return *this;
        }

        ~MinimizerHashTable() { clear_table(); }

        void clear_table() {

            if (table != nullptr) {

                delete[] table;
                //free(table);
                table = nullptr;
            }

            size_ = 0;
            pop  = 0;
            num_empty = 0;
        }

        size_t size() const { return pop; }

        bool empty() const { return pop == 0; }

        void clear() {

            std::fill(table, table+size_, empty_val);

            pop = 0;
            num_empty = size_;
        }

        void init_table(size_t sz) {

            clear_table();

            size_ = rndup(sz);

            table = new value_type[size_];
            //table = (value_type*) malloc(size_ * sizeof(value_type));

            clear();
        }

        iterator find(const Minimizer& key) {

            size_t h = hasher(key) & (size_-1);
            size_t end_h = (h == 0) ? (size_-1) : h-1;

            for (;; h =  (h+1!=size_ ? h+1 : 0)) {

                if (table[h].first == empty_val.first) return iterator(this); // empty slot, not in table
                else if (table[h].first == key) return iterator(this, h); // same key, found

                // if it is deleted, we still have to continue
                if (h==end_h) return iterator(this); // we've gone throught the table, quit
            }
        }

        const_iterator find(const Minimizer& key) const {

            size_t h = hasher(key) & (size_-1);
            size_t end_h = (h == 0) ? (size_-1) : h-1;

            for (;; h =  (h+1!=size_ ? h+1 : 0)) {

                if (table[h].first == empty_val.first) return const_iterator(this); // empty slot, not in table
                else if (table[h].first == key) return const_iterator(this, h); // same key, found

                if (h==end_h) return const_iterator(this);
            }
        }

        iterator find(const size_t h) {

            if ((h < size_) && (table[h].first != empty_val.first) && (table[h].first != deleted.first))
                return iterator(this, h);

            return iterator(this);
        }

        const_iterator find(const size_t h) const {

            if ((h < size_) && (table[h].first != empty_val.first) && (table[h].first != deleted.first))
                return const_iterator(this, h);

            return const_iterator(this);
        }

        iterator erase(const_iterator pos) {

            if (pos == this->end()) return this->end();

            table[pos.h] = deleted;
            --pop;

            return ++iterator(this, pos.h); // return pointer to next element
        }

        size_t erase(const Minimizer& minz) {

            const_iterator pos = find(minz);

            size_t oldpop = pop;
            if (pos != this->end()) erase(pos);

            return oldpop-pop;
        }

        std::pair<iterator,bool> insert(const value_type& val) {

            if ((5*num_empty) < size_) reserve(2*size_); // if more than 80% full, resize

            bool is_deleted = false;

            for (size_t h = hasher(val.first) & (size_-1), h_tmp;; h = (h+1 != size_ ? h+1 : 0)) {

                if (table[h].first == empty_val.first) {

                    if (!is_deleted) num_empty--;
                    else h = h_tmp;

                    table[h] = val;
                    ++pop;

                    return {iterator(this, h), true};
                }
                else if (table[h].first == val.first) return {iterator(this, h), false};
                else if (!is_deleted && (table[h].first == deleted.first)) {
                    is_deleted = true;
                    h_tmp = h;
                }
            }
        }

    void reserve(size_t sz) {

        if (sz <= size_) return;

        value_type *old_table = table;
        size_t old_size_ = size_;

        size_ = rndup(sz);
        pop = 0;
        num_empty = size_;

        table = new value_type[size_];

        std::fill(table, table+size_, empty_val);

        for (size_t i = 0; i < old_size_; i++) {

            if (old_table[i].first != empty_val.first && old_table[i].first != deleted.first) insert(old_table[i]);
        }

        free(old_table);
        old_table = nullptr;

        /*if (sz <= size_) return;

        const size_t prev_size_ = size_;

        size_ = rndup(sz);
        pop = 0;
        num_empty = size_;

        table = (value_type*) realloc(table, size_ * sizeof(value_type));

        std::sort(table, table + prev_size_, sortMinimizerHashTable(*this));

        value_type* table_empty = table;

        for (; table_empty < table + prev_size_; table_empty++) {

            if (((*table_empty).first == empty_val.first) || ((*table_empty).first == deleted.first)) break;
        }

        std::fill(table_empty, table + size_, empty_val);

        std::vector<value_type> v;

        for (int64_t i = table_empty - table - 1; i >= 0; i--) {

            size_t h = hasher(table[i].first) & (size_-1);

            if (h > i){

                for ( ; h < size_; h++) {

                    if (table[h].first == empty_val.first) {

                        num_empty--;
                        pop++;

                        std::swap(table[h], table[i]);

                        break;
                    }
                }

                if (h == size_){

                    v.push_back(table[i]);
                    table[i] = empty_val;
                }
            }
        }

        for (auto& vt : v) insert(vt);*/
    }

    size_t rndup(size_t v) {

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

    iterator begin() {

        iterator it(this);
        it.find_first();
        return it;
    }

    const_iterator begin() const {

        const_iterator it(this);
        it.find_first();
        return it;
    }

    iterator end() { return iterator(this); }

    const_iterator end() const { return const_iterator(this); }

    /*private:

        struct sortMinimizerHashTable {

            sortMinimizerHashTable(const MinimizerHashTable& mht_) : mht(mht_) {}

            bool operator() (const value_type& a, const value_type& b) const {

                const size_t h_a = (a.first == mht.empty_val.first) || (a.first == mht.deleted.first) ? 0xffffffffffffffff : mht.hasher(a.first) & (mht.size_-1);
                const size_t h_b = (b.first == mht.empty_val.first) || (b.first == mht.deleted.first) ? 0xffffffffffffffff : mht.hasher(b.first) & (mht.size_-1);

                return (h_a < h_b);
            }

            const MinimizerHashTable& mht;
        };*/
};

#endif // KALLISTO_KMERHASHTABLE_H
