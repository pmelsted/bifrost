#ifndef BIFROST_KMER_HASHTABLE_HPP
#define BIFROST_KMER_HASHTABLE_HPP

#include <utility>
#include <string>
#include <iterator>
#include <algorithm>

#include "Kmer.hpp"

/*template<typename T, typename Hash = KmerHash>
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
                operator++();
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

    #if defined(__AVX2__)
    __m256i empty_val256;
    #endif

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
                operator++();
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

            #if defined(__AVX2__)
            empty_val256 = _mm256_set1_epi64x(empty_val.first.longs[0]);
            #endif
        }

        MinimizerHashTable(size_t sz, const Hash& h = Hash() ) : hasher(h), table(nullptr), size_(0), pop(0), num_empty(0) {

            empty_val.first.set_empty();
            deleted.first.set_deleted();
            init_table((size_t) (1.2*sz));

            #if defined(__AVX2__)
            empty_val256 = _mm256_set1_epi64x(empty_val.first.longs[0]);
            #endif
        }

        MinimizerHashTable(MinimizerHashTable&& o){

            hasher = o.hasher;
            size_ = o.size_;
            pop = o.pop;
            num_empty = o.num_empty;
            table = o.table;
            empty_val = o.empty_val;
            deleted = o.deleted;

            #if defined(__AVX2__)
            empty_val256 = o.empty_val256;
            #endif

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

                #if defined(__AVX2__)
                empty_val256 = o.empty_val256;
                #endif

                o.table = nullptr;

                o.clear_table();
            }

            return *this;
        }

        ~MinimizerHashTable() { clear_table(); }

        void clear_table() {

            if (table != nullptr) {

                delete[] table;
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

            clear();
        }

        iterator find(const Minimizer& key) {

            const size_t end_table = size_-1;

            size_t h = hasher(key) & end_table;
            const size_t end_h = (h-1) & end_table;

            #if defined(__AVX2__)

            const __m256i key_val256 = _mm256_set1_epi64x(key.longs[0]);

            for (; h <= end_table - 4; h += 4){

                const __m256i table_val256 = _mm256_set_epi64x(table[h+3].first.longs[0], table[h+2].first.longs[0],
                                                               table[h+1].first.longs[0], table[h].first.longs[0]);

                const uint32_t pos = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi64(key_val256, table_val256)));

                if (pos != 0) return iterator(this, h + ((pos >> 15) & 0x1) + (((pos >> 23) & 0x1) << 1) + (pos >> 30));
                if (_mm256_movemask_epi8(_mm256_cmpeq_epi64(empty_val256, table_val256)) != 0) return iterator(this);
            }

            #endif

            for (; h != end_h; h = (h+1) & end_table) {

                if (table[h].first == empty_val.first) return iterator(this); // empty slot, not in table
                if (table[h].first == key) return iterator(this, h); // same key, found
            }

            return iterator(this);
        }

        const_iterator find(const Minimizer& key) const {

            const size_t end_table = size_ - 1;

            size_t h = hasher(key) & end_table;
            const size_t end_h = (h-1) & end_table;

            #if defined(__AVX2__)

            const __m256i key_val256 = _mm256_set1_epi64x(key.longs[0]);

            for (; h <= end_table - 4; h += 4){

                const __m256i table_val256 = _mm256_set_epi64x(table[h+3].first.longs[0], table[h+2].first.longs[0],
                                                               table[h+1].first.longs[0], table[h].first.longs[0]);

                const uint32_t pos = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi64(key_val256, table_val256)));

                if (pos != 0) return const_iterator(this, h + ((pos >> 15) & 0x1) + (((pos >> 23) & 0x1) << 1) + (pos >> 30));
                if (_mm256_movemask_epi8(_mm256_cmpeq_epi64(empty_val256, table_val256)) != 0) return const_iterator(this);
            }

            #endif

            for (; h != end_h; h = (h+1) & end_table) {

                if (table[h].first == empty_val.first) return const_iterator(this); // empty slot, not in table
                if (table[h].first == key) return const_iterator(this, h); // same key, found
            }

            return const_iterator(this);
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

            //for (size_t h = hasher(val.first) & (size_-1), h_tmp;; h = (h+1 != size_ ? h+1 : 0)) {
            for (size_t h = hasher(val.first) & (size_-1), h_tmp;; h = (h+1) & (size_-1)) {

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
};*/

template<typename T>
struct KmerHashTable {

    size_t size_, pop, num_empty;

    Kmer* table_keys;
    T* table_values;

    Kmer empty_key;
    Kmer deleted_key;

// ---- iterator ----
    template<bool is_const = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, T> {

        public:

            typedef typename std::conditional<is_const, const KmerHashTable *, KmerHashTable *>::type MHT_ptr_t;
            typedef typename std::conditional<is_const, const T&, T&>::type MHT_val_ref_t;
            typedef typename std::conditional<is_const, const T*, T*>::type MHT_val_ptr_t;

            MHT_ptr_t ht;
            size_t h;

            iterator_() : ht(nullptr), h(0) {}
            iterator_(MHT_ptr_t ht_) : ht(ht_), h(ht_->size_) {}
            iterator_(MHT_ptr_t ht_, size_t h_) :  ht(ht_), h(h_) {}
            iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h) {}
            iterator_& operator=(const iterator_& o) { ht=o.ht; h=o.h; return *this; }

            MHT_val_ref_t operator*() const { return ht->table_values[h]; }
            MHT_val_ptr_t operator->() const { return &(ht->table_values[h]); }

            const Kmer& getKey() const { return ht->table_keys[h]; }

            size_t getHash() const { return h; }

            void find_first() {

                h = 0;

                if ((ht != nullptr) && (ht->size_ > 0) &&
                    ((ht->table_keys[h] == ht->empty_key) || (ht->table_keys[h] == ht->deleted_key))) operator++();
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
        };

    typedef iterator_<true> const_iterator;
    typedef iterator_<false> iterator;

    // --- hash table
    KmerHashTable() : table_keys(nullptr), table_values(nullptr), size_(0), pop(0), num_empty(0) {

        empty_key.set_empty();
        deleted_key.set_deleted();

        init_tables(1024);
    }

    KmerHashTable(const size_t sz) : table_keys(nullptr), table_values(nullptr), size_(0), pop(0), num_empty(0) {

        empty_key.set_empty();
        deleted_key.set_deleted();

        init_tables(std::max(static_cast<size_t>(1.2 * sz), static_cast<size_t>(2)));
    }

    KmerHashTable(const KmerHashTable& o) : size_(o.size_), pop(o.pop), num_empty(o.num_empty),
                                            empty_key(o.empty_key), deleted_key(o.deleted_key) {

        table_keys = new Kmer[size_];
        table_values = new T[size_];

        std::copy(o.table_keys, o.table_keys + size_, table_keys);
        std::copy(o.table_values, o.table_values + size_, table_values);
    }

    KmerHashTable(KmerHashTable&& o){

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        empty_key = o.empty_key;
        deleted_key = o.deleted_key;

        table_keys = o.table_keys;
        table_values = o.table_values;

        o.table_keys = nullptr;
        o.table_values = nullptr;

        o.clear_tables();
    }

    KmerHashTable& operator=(const KmerHashTable& o) {

        clear_tables();

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        empty_key = o.empty_key;
        deleted_key = o.deleted_key;

        table_keys = new Kmer[size_];
        table_values = new T[size_];

        std::copy(o.table_keys, o.table_keys + size_, table_keys);
        std::copy(o.table_values, o.table_values + size_, table_values);

        return *this;
    }

    KmerHashTable& operator=(KmerHashTable&& o){

        if (this != &o) {

            clear_tables();

            size_ = o.size_;
            pop = o.pop;
            num_empty = o.num_empty;

            empty_key = o.empty_key;
            deleted_key = o.deleted_key;

            table_keys = o.table_keys;
            table_values = o.table_values;

            o.table_keys = nullptr;
            o.table_values = nullptr;

            o.clear_tables();
        }

        return *this;
    }

    ~KmerHashTable() { clear_tables(); }

    inline size_t size() const { return pop; }

    inline bool empty() const { return pop == 0; }

    void clear() {

        std::fill(table_keys, table_keys + size_, empty_key);

        pop = 0;
        num_empty = size_;
    }

    void clear_tables() {

        if (table_keys != nullptr) {

            delete[] table_keys;
            table_keys = nullptr;
        }

        if (table_values != nullptr) {

            delete[] table_values;
            table_values = nullptr;
        }

        size_ = 0;
        pop  = 0;
        num_empty = 0;
    }

    void init_tables(const size_t sz) {

        clear_tables();

        size_ = rndup(sz);

        table_keys = new Kmer[size_];
        table_values = new T[size_];

        clear();
    }

    void reserve(const size_t sz) {

        if (sz <= size_) return;

        Kmer* old_table_keys = table_keys;
        T* old_table_values = table_values;

        const size_t old_size_ = size_;

        size_ = rndup(sz);
        pop = 0;
        num_empty = size_;

        table_keys = new Kmer[size_];
        table_values = new T[size_];

        std::fill(table_keys, table_keys + size_, empty_key);

        for (size_t i = 0; i < old_size_; ++i) {

            if (old_table_keys[i] != empty_key && old_table_keys[i] != deleted_key){

                insert(std::move(old_table_keys[i]), std::move(old_table_values[i]));
            }
        }

        delete[] old_table_keys;
        delete[] old_table_values;
    }

    iterator find(const Kmer& key) {

        const size_t end_table = size_-1;
        size_t h = key.hash() & end_table;
        const size_t end_h = (h-1) & end_table;

        for (; h != end_h; h = (h+1) & end_table) {

            if ((table_keys[h] == empty_key) || (table_keys[h] == key)) break;
        }

        if ((h != end_h) && (table_keys[h] == key)) return iterator(this, h);

        return iterator(this);
    }

    const_iterator find(const Kmer& key) const {

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

    iterator erase(const_iterator pos) {

        if (pos == end()) return end();

        table_keys[pos.h] = deleted_key;
        --pop;

        return ++iterator(this, pos.h); // return pointer to next element
    }

    size_t erase(const Kmer& minz) {

        const_iterator pos = find(minz);

        size_t oldpop = pop;

        if (pos != end()) erase(pos);

        return oldpop - pop;
    }

    std::pair<iterator, bool> insert(const Kmer& key, const T& value) {

        if ((5 * num_empty) < size_) reserve(2 * size_); // if more than 80% full, resize

        bool is_deleted = false;

        const size_t end_table = size_-1;

        for (size_t h = key.hash() & end_table, h_tmp;; h = (h+1) & end_table) {

            if (table_keys[h] == empty_key) {

                is_deleted ? h = h_tmp : --num_empty;

                table_keys[h] = key;
                table_values[h] = value;

                ++pop;

                return {iterator(this, h), true};
            }
            else if (table_keys[h] == key) return {iterator(this, h), false};
            else if (!is_deleted && (table_keys[h] == deleted_key)) {

                is_deleted = true;
                h_tmp = h;
            }
        }
    }

    std::pair<iterator, bool> insert(Kmer&& key, T&& value) {

        if ((5 * num_empty) < size_) reserve(2 * size_); // if more than 80% full, resize

        bool is_deleted = false;

        const size_t end_table = size_-1;

        for (size_t h = key.hash() & (size_-1), h_tmp;; h = (h+1) & end_table) {

            if (table_keys[h] == empty_key) {

                is_deleted ? h = h_tmp : --num_empty;

                table_keys[h] = key;
                table_values[h] = value;

                ++pop;

                return {iterator(this, h), true};
            }
            else if (table_keys[h] == key) return {iterator(this, h), false};
            else if (!is_deleted && (table_keys[h] == deleted_key)) {

                is_deleted = true;
                h_tmp = h;
            }
        }
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
};

template<typename T>
struct MinimizerHashTable {

    size_t size_, pop, num_empty;

    Minimizer* table_keys;
    T* table_values;

    Minimizer empty_key;
    Minimizer deleted_key;

// ---- iterator ----
    template<bool is_const = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, T> {

        public:

            typedef typename std::conditional<is_const, const MinimizerHashTable *, MinimizerHashTable *>::type MHT_ptr_t;
            typedef typename std::conditional<is_const, const T&, T&>::type MHT_val_ref_t;
            typedef typename std::conditional<is_const, const T*, T*>::type MHT_val_ptr_t;

            MHT_ptr_t ht;
            size_t h;

            iterator_() : ht(nullptr), h(0) {}
            iterator_(MHT_ptr_t ht_) : ht(ht_), h(ht_->size_) {}
            iterator_(MHT_ptr_t ht_, size_t h_) :  ht(ht_), h(h_) {}
            iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h) {}
            iterator_& operator=(const iterator_& o) { ht=o.ht; h=o.h; return *this; }

            MHT_val_ref_t operator*() const { return ht->table_values[h]; }
            MHT_val_ptr_t operator->() const { return &(ht->table_values[h]); }

            const Minimizer& getKey() const { return ht->table_keys[h]; }

            size_t getHash() const { return h; }

            void find_first() {

                h = 0;

                if ((ht != nullptr) && (ht->size_ > 0) &&
                    ((ht->table_keys[h] == ht->empty_key) || (ht->table_keys[h] == ht->deleted_key))) operator++();
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
        };

    typedef iterator_<true> const_iterator;
    typedef iterator_<false> iterator;

    // --- hash table
    MinimizerHashTable() : table_keys(nullptr), table_values(nullptr), size_(0), pop(0), num_empty(0) {

        empty_key.set_empty();
        deleted_key.set_deleted();

        init_tables(1024);
    }

    MinimizerHashTable(const size_t sz) : table_keys(nullptr), table_values(nullptr), size_(0), pop(0), num_empty(0) {

        empty_key.set_empty();
        deleted_key.set_deleted();

        init_tables(std::max(static_cast<size_t>(1.2 * sz), static_cast<size_t>(2)));
    }

    MinimizerHashTable(const MinimizerHashTable& o) :   size_(o.size_), pop(o.pop), num_empty(o.num_empty),
                                                        empty_key(o.empty_key), deleted_key(o.deleted_key) {

        table_keys = new Minimizer[size_];
        table_values = new T[size_];

        std::copy(o.table_keys, o.table_keys + size_, table_keys);
        std::copy(o.table_values, o.table_values + size_, table_values);
    }

    MinimizerHashTable(MinimizerHashTable&& o){

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        empty_key = o.empty_key;
        deleted_key = o.deleted_key;

        table_keys = o.table_keys;
        table_values = o.table_values;

        o.table_keys = nullptr;
        o.table_values = nullptr;

        o.clear_tables();
    }

    MinimizerHashTable& operator=(const MinimizerHashTable& o) {

        clear_tables();

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        empty_key = o.empty_key;
        deleted_key = o.deleted_key;

        table_keys = new Minimizer[size_];
        table_values = new T[size_];

        std::copy(o.table_keys, o.table_keys + size_, table_keys);
        std::copy(o.table_values, o.table_values + size_, table_values);

        return *this;
    }

    MinimizerHashTable& operator=(MinimizerHashTable&& o){

        if (this != &o) {

            clear_tables();

            size_ = o.size_;
            pop = o.pop;
            num_empty = o.num_empty;

            empty_key = o.empty_key;
            deleted_key = o.deleted_key;

            table_keys = o.table_keys;
            table_values = o.table_values;

            o.table_keys = nullptr;
            o.table_values = nullptr;

            o.clear_tables();
        }

        return *this;
    }

    ~MinimizerHashTable() { clear_tables(); }

    inline size_t size() const { return pop; }

    inline bool empty() const { return pop == 0; }

    void clear() {

        std::fill(table_keys, table_keys + size_, empty_key);

        pop = 0;
        num_empty = size_;
    }

    void clear_tables() {

        if (table_keys != nullptr) {

            delete[] table_keys;
            table_keys = nullptr;
        }

        if (table_values != nullptr) {

            delete[] table_values;
            table_values = nullptr;
        }

        size_ = 0;
        pop  = 0;
        num_empty = 0;
    }

    void init_tables(const size_t sz) {

        clear_tables();

        size_ = rndup(sz);

        table_keys = new Minimizer[size_];
        table_values = new T[size_];

        clear();
    }

    void reserve(const size_t sz) {

        if (sz <= size_) return;

        Minimizer* old_table_keys = table_keys;
        T* old_table_values = table_values;

        size_t old_size_ = size_;

        size_ = rndup(sz);
        pop = 0;
        num_empty = size_;

        table_keys = new Minimizer[size_];
        table_values = new T[size_];

        std::fill(table_keys, table_keys + size_, empty_key);

        for (size_t i = 0; i < old_size_; i++) {

            if (old_table_keys[i] != empty_key && old_table_keys[i] != deleted_key){

                insert(std::move(old_table_keys[i]), std::move(old_table_values[i]));
            }
        }

        delete[] old_table_keys;
        delete[] old_table_values;
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

    iterator erase(const_iterator pos) {

        if (pos == end()) return end();

        table_keys[pos.h] = deleted_key;
        --pop;

        return ++iterator(this, pos.h); // return pointer to next element
    }

    size_t erase(const Minimizer& minz) {

        const_iterator pos = find(minz);

        size_t oldpop = pop;

        if (pos != end()) erase(pos);

        return oldpop - pop;
    }

    std::pair<iterator, bool> insert(const Minimizer& key, const T& value) {

        if ((5 * num_empty) < size_) reserve(2 * size_); // if more than 80% full, resize

        bool is_deleted = false;

        const size_t end_table = size_-1;

        for (size_t h = key.hash() & (size_-1), h_tmp;; h = (h+1) & end_table) {

            if (table_keys[h] == empty_key) {

                is_deleted ? h = h_tmp : --num_empty;

                table_keys[h] = key;
                table_values[h] = value;

                ++pop;

                return {iterator(this, h), true};
            }
            else if (table_keys[h] == key) return {iterator(this, h), false};
            else if (!is_deleted && (table_keys[h] == deleted_key)) {

                is_deleted = true;
                h_tmp = h;
            }
        }
    }

    std::pair<iterator, bool> insert(Minimizer&& key, T&& value) {

        if ((5 * num_empty) < size_) reserve(2 * size_); // if more than 80% full, resize

        bool is_deleted = false;

        const size_t end_table = size_-1;

        for (size_t h = key.hash() & (size_-1), h_tmp;; h = (h+1) & end_table) {

            if (table_keys[h] == empty_key) {

                is_deleted ? h = h_tmp : --num_empty;

                table_keys[h] = key;
                table_values[h] = value;

                ++pop;

                return {iterator(this, h), true};
            }
            else if (table_keys[h] == key) return {iterator(this, h), false};
            else if (!is_deleted && (table_keys[h] == deleted_key)) {

                is_deleted = true;
                h_tmp = h;
            }
        }
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
};

#endif
