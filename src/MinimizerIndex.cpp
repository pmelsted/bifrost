#include "MinimizerIndex.hpp"

MinimizerIndex::MinimizerIndex() :  table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                    size_(0), pop(0), num_empty(0)  {

    init_tables(max(static_cast<size_t>(1024), lck_block_sz));
}

MinimizerIndex::MinimizerIndex(const size_t sz) :   table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                                    size_(0), pop(0), num_empty(0) {

    if (sz == 0) init_tables(lck_block_sz);
    else {

        const size_t sz_with_empty = static_cast<size_t>(1.2 * sz);

        size_t rdnup_sz = rndup(sz);

        while (rdnup_sz < sz_with_empty) rdnup_sz <<= 1;

        init_tables(max(rdnup_sz, lck_block_sz));
    }
}

MinimizerIndex::MinimizerIndex(const MinimizerIndex& o) :   size_(o.size_), pop(o.pop), num_empty(o.num_empty) {

    table_keys = new Minimizer[size_];
    table_tinyv = new packed_tiny_vector[size_];
    table_tinyv_sz = new uint8_t[size_];

    lck_min = vector<SpinLock>(o.lck_min.size());

    std::copy(o.table_keys, o.table_keys + size_, table_keys);

    for (size_t i = 0; i < size_; ++i){

        table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;
        table_tinyv[i].copy(table_tinyv_sz[i], o.table_tinyv[i], o.table_tinyv_sz[i]);
    }
}

MinimizerIndex::MinimizerIndex(MinimizerIndex&& o){

    size_ = o.size_;
    pop = o.pop;
    num_empty = o.num_empty;

    table_keys = o.table_keys;
    table_tinyv = o.table_tinyv;
    table_tinyv_sz = o.table_tinyv_sz;

    lck_min = vector<SpinLock>(o.lck_min.size());

    o.table_keys = nullptr;
    o.table_tinyv = nullptr;
    o.table_tinyv_sz = nullptr;

    o.clear();
}

MinimizerIndex& MinimizerIndex::operator=(const MinimizerIndex& o) {

    if (this != &o) {

        clear();

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        table_keys = new Minimizer[size_];
        table_tinyv = new packed_tiny_vector[size_];
        table_tinyv_sz = new uint8_t[size_];

        lck_min = vector<SpinLock>(o.lck_min.size());

        std::copy(o.table_keys, o.table_keys + size_, table_keys);

        for (size_t i = 0; i < size_; ++i){

            table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;
            table_tinyv[i].copy(table_tinyv_sz[i], o.table_tinyv[i], o.table_tinyv_sz[i]);
        }
    }

    return *this;
}

MinimizerIndex& MinimizerIndex::operator=(MinimizerIndex&& o){

    if (this != &o) {

        clear();

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        table_keys = o.table_keys;
        table_tinyv = o.table_tinyv;
        table_tinyv_sz = o.table_tinyv_sz;

        lck_min = vector<SpinLock>(o.lck_min.size());

        o.table_keys = nullptr;
        o.table_tinyv = nullptr;
        o.table_tinyv_sz = nullptr;

        o.clear();
    }

    return *this;
}

MinimizerIndex::~MinimizerIndex() {

    clear();
}

void MinimizerIndex::clear() {

    if (table_tinyv != nullptr){

        for (size_t i = 0; i < size_; ++i) table_tinyv[i].destruct(table_tinyv_sz[i]);
    }

    clear_tables();

    lck_min.clear();
    lck_edit_table.release_all();
}

MinimizerIndex::iterator MinimizerIndex::find(const Minimizer& key) {

    const size_t end_table = size_-1;

    size_t h = key.hash() & end_table;
    size_t i = 0;

    while (i != size_) {

        if (table_keys[h].isEmpty() || (table_keys[h] == key)) break;

        h = (h+1) & end_table;
        ++i;
    }

    if ((i != size_) && (table_keys[h] == key)) return iterator(this, h);

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const Minimizer& key) const {

    const size_t end_table = size_-1;

    size_t h = key.hash() & end_table;
    size_t i = 0;

    while (i != size_) {

        if (table_keys[h].isEmpty() || (table_keys[h] == key)) break;

        h = (h+1) & end_table;
        ++i;
    }

    if ((i != size_) && (table_keys[h] == key)) return const_iterator(this, h);

    return const_iterator(this);
}

MinimizerIndex::iterator MinimizerIndex::find(const size_t h) {

    if ((h < size_) && !table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return iterator(this, h);

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const size_t h) const {

    if ((h < size_) && !table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return const_iterator(this, h);

    return const_iterator(this);
}

MinimizerIndex::iterator MinimizerIndex::erase(const_iterator it) {

    if (it == end()) return end();

    table_keys[it.h].set_deleted();
    table_tinyv[it.h].destruct(table_tinyv_sz[it.h]);
    table_tinyv_sz[it.h] = packed_tiny_vector::FLAG_EMPTY;

    --pop;

    return ++iterator(this, it.h); // return pointer to next element
}

size_t MinimizerIndex::erase(const Minimizer& minz) {

    const size_t end_table = size_-1;
    const size_t oldpop = pop;

    size_t h = minz.hash() & end_table;
    size_t i = 0;

    while (i != size_) {

        if (table_keys[h].isEmpty() || (table_keys[h] == minz)) break;

        h = (h+1) & end_table;
        ++i;
    }

    if ((i != size_) && (table_keys[h] == minz)){

        table_keys[h].set_deleted();
        table_tinyv[h].destruct(table_tinyv_sz[h]);
        table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;

        --pop;
    }

    return oldpop - pop;
}

pair<MinimizerIndex::iterator, bool> MinimizerIndex::insert(const Minimizer& key, const packed_tiny_vector& ptv, const uint8_t& flag) {

    if ((5 * num_empty) < size_) reserve(2 * size_); // if more than 80% full, resize

    const size_t end_table = size_-1;
    
    size_t h = key.hash() & end_table, h_del;

    bool is_deleted = false;

    while (true) {

        if (table_keys[h].isEmpty()) {

            h = ((static_cast<size_t>(is_deleted) - 1) & h) + ((static_cast<size_t>(!is_deleted) - 1) & h_del);
            num_empty -= static_cast<size_t>(!is_deleted);

            table_keys[h] = key;
            table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;

            table_tinyv[h].copy(table_tinyv_sz[h], ptv, flag);

            ++pop;

            return {iterator(this, h), true};
        }
        else if (table_keys[h] == key) return {iterator(this, h), false};
        else if (table_keys[h].isDeleted()) {

            h_del = ((static_cast<size_t>(!is_deleted) - 1) & h_del) + ((static_cast<size_t>(is_deleted) - 1) & h);
            is_deleted = true;
        }

        h = (h+1) & end_table;
    }
}

void MinimizerIndex::init_threads() {

    lck_min = vector<SpinLock>((size_ + lck_block_sz - 1) / lck_block_sz);

    pop_p = pop;
    num_empty_p = num_empty;
}

void MinimizerIndex::release_threads() {

    pop = pop_p;
    num_empty = num_empty_p;

    lck_min.clear();
    lck_edit_table.release_all();
}

MinimizerIndex::iterator MinimizerIndex::find_p(const Minimizer& key) {

    lck_edit_table.acquire_reader();

    const size_t end_table = size_-1;

    size_t i = 0;
    size_t h = key.hash() & end_table;
    size_t id_block = h >> lck_block_div_shift;

    lck_min[id_block].acquire();

    while (i != size_) {

        if ((h >> lck_block_div_shift) != id_block){

            lck_min[id_block].release();
            id_block = h >> lck_block_div_shift;
            lck_min[id_block].acquire();
        }

        if (table_keys[h].isEmpty()){

            lck_min[id_block].release();
            lck_edit_table.release_reader();

            iterator(this);
        }
        else if (table_keys[h] == key) return iterator(this, h);

        h = (h+1) & end_table;
        ++i;
    }

    lck_min[id_block].release();
    lck_edit_table.release_reader();

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find_p(const Minimizer& key) const {

    lck_edit_table.acquire_reader();

    const size_t end_table = size_-1;

    size_t i = 0;
    size_t h = key.hash() & end_table;
    size_t id_block = h >> lck_block_div_shift;

    lck_min[id_block].acquire();

    while (i != size_) {

        if ((h >> lck_block_div_shift) != id_block){

            lck_min[id_block].release();
            id_block = h >> lck_block_div_shift;
            lck_min[id_block].acquire();
        }

        if (table_keys[h].isEmpty()){

            lck_min[id_block].release();
            lck_edit_table.release_reader();

            const_iterator(this);
        }
        else if (table_keys[h] == key) return const_iterator(this, h);

        h = (h+1) & end_table;
        ++i;
    }

    lck_min[id_block].release();
    lck_edit_table.release_reader();

    return const_iterator(this);
}

MinimizerIndex::iterator MinimizerIndex::find_p(const size_t h) {

    lck_edit_table.acquire_reader();

    if (h < size_){

        const size_t id_block = h >> lck_block_div_shift;

        lck_min[id_block].acquire();

        if (!table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return iterator(this, h);

        lck_min[id_block].release();
    }

    lck_edit_table.release_reader();

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find_p(const size_t h) const {

    lck_edit_table.acquire_reader();

    if (h < size_){

        const size_t id_block = h >> lck_block_div_shift;

        lck_min[id_block].acquire();

        if (!table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return const_iterator(this, h);

        lck_min[id_block].release();
    }

    lck_edit_table.release_reader();

    return const_iterator(this);
}

void MinimizerIndex::release_p(const_iterator it) const {

    if (it != end()){

        lck_min[it.h >> lck_block_div_shift].release();
        lck_edit_table.release_reader();
    }
}

void MinimizerIndex::release_p(iterator it) {

    if (it != end()){

        lck_min[it.h >> lck_block_div_shift].release();
        lck_edit_table.release_reader();
    }
}

size_t MinimizerIndex::erase_p(const Minimizer& minz) {

    lck_edit_table.acquire_reader();

    const size_t end_table = size_ - 1;

    size_t i = 0;
    size_t h = minz.hash() & end_table;
    size_t id_block = h >> lck_block_div_shift;
    size_t l_pop = pop;

    lck_min[id_block].acquire();

    while (i != size_) {

        if ((h >> lck_block_div_shift) != id_block){

            lck_min[id_block].release();
            id_block = h >> lck_block_div_shift;
            lck_min[id_block].acquire();
        }

        if (table_keys[h].isEmpty()){

            lck_min[id_block].release();
            lck_edit_table.release_reader();

            return 0;
        }
        else if (table_keys[h] == minz) break;

        h = (h+1) & end_table;
        ++i;
    }

    if ((i != size_) && (table_keys[h] == minz)){

        table_keys[h].set_deleted();
        table_tinyv[h].destruct(table_tinyv_sz[h]);
        table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;

        lck_min[id_block].release();

        --pop;
    }

    l_pop -= pop;

    lck_edit_table.release_reader();

    return l_pop;
}

pair<MinimizerIndex::iterator, bool> MinimizerIndex::insert_p(const Minimizer& key, const packed_tiny_vector& v, const uint8_t& flag) {

    bool is_deleted = false;

    lck_edit_table.acquire_reader();

    if ((5 * num_empty_p) < size_){

        lck_edit_table.release_reader();
        lck_edit_table.acquire_writer();

        reserve(2 * size_); // if more than 80% full, resize

        pop_p = pop;
        num_empty_p = num_empty;

        lck_edit_table.release_writer_acquire_reader();
    }

    const size_t end_table = size_-1;
    const size_t h = key.hash() & end_table;

    size_t id_block = h >> lck_block_div_shift;

    size_t h1 = h;
    size_t h2;

    lck_min[id_block].acquire();

    while (true) {

        if ((h1 >> lck_block_div_shift) != id_block){

            lck_min[id_block].release();
            id_block = h1 >> lck_block_div_shift;
            lck_min[id_block].acquire();
        }

        if (table_keys[h1].isEmpty()) {

            if (is_deleted){

                const size_t id_block2 = h2 >> lck_block_div_shift;

                if (id_block2 == id_block) h1 = h2;
                else {

                    lck_min[id_block2].acquire();

                    if (table_keys[h2].isDeleted()){

                        lck_min[id_block].release();
                        h1 = h2;
                    }
                    else {

                        lck_min[id_block2].release();
                        --num_empty_p;
                    }
                }
            }
            else --num_empty_p;

            table_keys[h1] = key;
            table_tinyv_sz[h1] = packed_tiny_vector::FLAG_EMPTY;

            table_tinyv[h1].copy(table_tinyv_sz[h1], v, flag);

            ++pop_p;

            return {iterator(this, h1), true};
        }
        else if (table_keys[h1] == key){

            return {iterator(this, h1), false};
        }
        else if (!is_deleted && table_keys[h1].isDeleted()) {

            is_deleted = true;
            h2 = h1;
        }

        h1 = (h1+1) & end_table;
    }

    lck_min[id_block].release(); // Just for safety
    lck_edit_table.release_reader(); // Just for safety
}

MinimizerIndex::iterator MinimizerIndex::begin() {

    iterator it(this, 0xffffffffffffffffULL);
    it.operator++();
    return it;
}

MinimizerIndex::const_iterator MinimizerIndex::begin() const {

    const_iterator it(this, 0xffffffffffffffffULL);
    it.operator++();
    return it;
}

MinimizerIndex::iterator MinimizerIndex::end() {

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::end() const {

    return const_iterator(this);
}

void MinimizerIndex::clear_tables() {

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

void MinimizerIndex::init_tables(const size_t sz) {

    clear_tables();

    Minimizer empty_key;

    pop = 0;
    size_ = rndup(sz);
    num_empty = size_;

    table_keys = new Minimizer[size_];
    table_tinyv = new packed_tiny_vector[size_];
    table_tinyv_sz = new uint8_t[size_];

    empty_key.set_empty();

    std::fill(table_keys, table_keys + size_, empty_key);

    memset(table_tinyv_sz, packed_tiny_vector::FLAG_EMPTY, size_ * sizeof(uint8_t));
}

void MinimizerIndex::reserve(const size_t sz) {

    if (sz <= size_) return;

    const size_t old_size_ = size_;

    Minimizer empty_key;

    Minimizer* old_table_keys = table_keys;
    packed_tiny_vector* old_table_tinyv = table_tinyv;
    uint8_t* old_table_tinyv_sz = table_tinyv_sz;

    size_ = rndup(sz);
    pop = 0;
    num_empty = size_;

    table_keys = new Minimizer[size_];
    table_tinyv = new packed_tiny_vector[size_];
    table_tinyv_sz = new uint8_t[size_];

    empty_key.set_empty();

    if (!lck_min.empty()) lck_min = vector<SpinLock>((size_ + lck_block_sz - 1) / lck_block_sz);

    std::fill(table_keys, table_keys + size_, empty_key);

    memset(table_tinyv_sz, packed_tiny_vector::FLAG_EMPTY, size_ * sizeof(uint8_t));

    for (size_t i = 0; i < old_size_; ++i) {

        if (!old_table_keys[i].isEmpty() && !old_table_keys[i].isDeleted()){

            insert(old_table_keys[i], old_table_tinyv[i], old_table_tinyv_sz[i]);

            old_table_tinyv[i].destruct(old_table_tinyv_sz[i]);
        }
    }

    delete[] old_table_keys;
    delete[] old_table_tinyv;
    delete[] old_table_tinyv_sz;
}

const size_t MinimizerIndex::lck_block_sz = 64;
const size_t MinimizerIndex::lck_block_div_shift = 6;

/*CompactedMinimizerIndex::CompactedMinimizerIndex() :    table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr), table_hash_bits(nullptr) {

    clear();
}

CompactedMinimizerIndex::CompactedMinimizerIndex(const MinimizerIndex& mi, const size_t nb_hashes, const size_t nb_threads) :
                                                table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr), table_hash_bits(nullptr) {

    clear();

    if (nb_hashes == 0) {

        cerr << "CompactedMinimizerIndex::CompactedMinimizerIndex(): Number of hash functions to use cannot be 0." << endl;
        sys.exit(1);
    }

    if (nb_hashes > 64) {

        cerr << "CompactedMinimizerIndex::CompactedMinimizerIndex(): Number of hash functions to use cannot exceed 64." << endl;
        sys.exit(1);
    }

    if (nb_threads == 0) {

        cerr << "CompactedMinimizerIndex::CompactedMinimizerIndex(): Number of threads to use cannot be 0." << endl;
        sys.exit(1);
    }

    size = mi.size();
    nb_h = nb_hashes;

    if (size != 0) {

        hbits_per_minz = rndup(nb_h);

        table_keys = new MinimizerIndex[size];
        table_tinyv = new packed_tiny_vector[size];
        table_tinyv_sz = new uint8_t[size];
        table_hash_bits = new uint8_t[(hbits_per_minz * size + 63) / 64];

        for (size_t i = 0; i < size; ++i) table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;

        if (nb_threads == 1) {

            MinimizerIndex::const_iterator its = mi.begin(); ite = mi.end();

            while (its != ite) {

                const Minimizer minz = its.getKey();

                size_t hseed = 0;
                uint64_t h_pos = 0, h_bit_pos = 0;
                bool is_used = true;

                while ((hseed != nb_hashes) && is_used) {

                    h_pos = minz.hash(hseed) % size;
                    h_bit_pos = h_pos * bits_per_hash + hseed;
                    is_used = static_cast<bool>(table_hash_bits[h_bit_pos >> 6] & (1ULL << (h_bit_pos & 0x3fULL)));
                    hseed += static_cast<uint64_t>(is_used);
                }

                if (hseed != nb_hashes) {

                    table_keys[h_pos] = minz;
                    table_hash_bits[h_bit_pos >> 6] |= (1ULL << (h_bit_pos & 0x3fULL));
                    table_tinyv[h_pos].copy(table_tinyv_sz[h_pos], its.getVector(), its.getVectorSize());
                }
                else mi_overflow.insert(minz, its.getVector(), its.getVectorSize());
            }
        }
        else {

            SpinLock* table_splk = new SpinLock[(size + 1023) / 1024];
            SpinLock splk_overflow;

            auto insert = [&](MinimizerIndex::const_iterator& its, MinimizerIndex::const_iterator& ite) {

                while (its != ite) {

                    const Minimizer minz = its.getKey();

                    size_t hseed = 0, pos_splk = 0;
                    uint64_t h_pos = 0, h_bit_pos = 0;
                    bool is_used = true;

                    while (hseed != nb_hashes) {

                        h_pos = minz.hash(hseed) % size;
                        h_bit_pos = h_pos * hbits_per_minz + hseed;

                        pos_splk = h_pos / 1024;

                        table_splk[pos_splk].acquire();

                        is_used = static_cast<bool>(table_hash_bits[h_bit_pos >> 6] & (1ULL << (h_bit_pos & 0x3fULL)));

                        if (!is_used) break;

                        table_splk[pos_splk].release();

                        ++hseed;
                    }

                    if (hseed != nb_hashes) {

                        table_keys[h_pos] = minz;
                        table_hash_bits[h_bit_pos >> 6] |= (1ULL << (h_bit_pos & 0x3fULL));
                        table_tinyv[h_pos].copy(table_tinyv_sz[h_pos], its.getVector(), its.getVectorSize());

                        table_splk[pos_splk].release();
                    }
                    else {

                        splk_overflow.acquire();

                        mi_overflow.insert(minz, its.getVector(), its.getVectorSize());

                        splk_overflow.release();
                    }
                }
            };

            {
                bool stop = false;

                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mutex_mi;

                MinimizerIndex::const_iterator its = mi.begin(); ite = mi.end();

                for (size_t t = 0; t < nb_threads; ++t) {

                    workers.emplace_back(

                        [&]{

                            MinimizerIndex::const_iterator lits;
                            MinimizerIndex::const_iterator lite;

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_mi);

                                    if (its == ite) return;

                                    lits = its;
                                    lite = its;

                                    for (size_t i = 0; (i < 10000) && (lite != ite); ++i) ++lite;

                                    its = lite;
                                }

                                insert(lits, lite);
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }
        }
    }
}

CompactedMinimizerIndex::CompactedMinimizerIndex(MinimizerIndex&& mi, const size_t nb_hashes, const size_t nb_threads) :
                                                CompactedMinimizerIndex(static_cast<const MinimizerIndex&>(mi), nb_hashes, nb_threads) {

    mi.clear();
}

CompactedMinimizerIndex::~CompactedMinimizerIndex() {

    clear();
}

CompactedMinimizerIndex& CompactedMinimizerIndex::operator=(const CompactedMinimizerIndex& cmi) {

    if (this != &cmi) {

        clear();

        size = cmi.size;
        nb_h = cmi.nb_h;
        hbits_per_minz = cmi.bits_per_hash;

        if (cmi.table_keys != nullptr) {

            table_keys = new MinimizerIndex[size];

            std::copy(cmi.table_keys, cmi.table_keys + size, table_keys);
        }

        if (cmi.table_tinyv != nullptr) && (cmi.table_tinyv_sz != nullptr) {

            table_tinyv = new packed_tiny_vector[size];
            table_tinyv_sz = new uint8_t[size];

            for (size_t i = 0; i < size; ++i){

                table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;
                table_tinyv[i].copy(table_tinyv_sz[i], cmi.table_tinyv[i], cmi.table_tinyv_sz[i]);
            }
        }

        if (cmi.table_hash_bits != nullptr) {

            const size_t table_hash_bits_sz = (size * hbits_per_minz + 63) / 64;

            table_hash_bits = new uint64_t[table_hash_bits_sz];

            std::copy(cmi.table_hash_bits, cmi.table_hash_bits + table_hash_bits_sz, table_hash_bits);
        }

        mi_overflow = cmi.mi_overflow;
    }

    return *this;
}

CompactedMinimizerIndex& CompactedMinimizerIndex::operator=(CompactedMinimizerIndex&& cmi){

    if (this != &cmi) {

        clear();

        size = cmi.size;
        nb_h = cmi.nb_h;
        hbits_per_minz = cmi.bits_per_hash;

        table_keys = cmi.table_keys;
        table_tinyv = cmi.table_tinyv;
        table_tinyv_sz = cmi.table_tinyv_sz;
        table_hash_bits = cmi.table_hash_bits;

        mi_overflow = move(cmi.mi_overflow);

        cmi.table_keys = nullptr;
        cmi.table_tinyv = nullptr;
        cmi.table_tinyv_sz = nullptr;
        cmi.table_hash_bits = nullptr;

        cmi.clear();
    }

    return *this;
}

void CompactedMinimizerIndex::clear() {

    if (table_keys != nullptr) {

        delete[] table_keys;
        table_keys = nullptr;
    }

    if (table_tinyv != nullptr){

        for (size_t i = 0; i < size_; ++i) table_tinyv[i].destruct(table_tinyv_sz[i]);

        delete[] table_tinyv;
        delete[] table_tinyv_sz;

        table_tinyv = nullptr;
        table_tinyv_sz = nullptr;
    }

    if (table_hash_bits != nullptr) {

        delete[] table_hash_bits;
        table_hash_bits = nullptr;
    }

    mi_overflow.clear();

    size = 0;
    nb_h = 0;
    hbits_per_minz = 0;
}

MinimizerIndex::iterator MinimizerIndex::find(const Minimizer& key) {

    for (size_t hseed = 0; hseed != nb_hashes; ++hseed) {

        const uint64_t h_pos = key.hash(hseed) % size;

        for (uint64_t h_bit_pos = h_pos * bits_per_hash; h_bit_pos < (h_pos * bits_per_hash + nb_hashes); ++h_bit_pos) {

            const bool is_used = static_cast<bool>(table_hash_bits[h_bit_pos >> 6] & (1ULL << (h_bit_pos & 0x3fULL)));

            if (is_used && (table_keys[h_pos] == key)) return iterator(this, h_pos);
        }
    }

    return iterator(this, size, mi_overflow.find(minz));
}

MinimizerIndex::const_iterator MinimizerIndex::find(const Minimizer& key) const {

    for (size_t hseed = 0; hseed != nb_hashes; ++hseed) {

        const uint64_t h_pos = key.hash(hseed) % size;

        for (uint64_t h_bit_pos = h_pos * bits_per_hash; h_bit_pos < (h_pos * bits_per_hash + nb_hashes); ++h_bit_pos) {

            const bool is_used = static_cast<bool>(table_hash_bits[h_bit_pos >> 6] & (1ULL << (h_bit_pos & 0x3fULL)));

            if (is_used && (table_keys[h_pos] == key)) return const_iterator(this, h_pos);
        }
    }

    return const_iterator(this, size, mi_overflow.find(minz));
}

MinimizerIndex::iterator MinimizerIndex::find(const size_t h) {

    if ((h < size_) && !table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return iterator(this, h);

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const size_t h) const {

    if ((h < size_) && !table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return const_iterator(this, h);

    return const_iterator(this);
}

CompactedMinimizerIndex::iterator CompactedMinimizerIndex::begin() {

    iterator it(this, 0xffffffffffffffffULL);
    it.operator++();
    return it;
}

CompactedMinimizerIndex::const_iterator CompactedMinimizerIndex::begin() const {

    const_iterator it(this, 0xffffffffffffffffULL);
    it.operator++();
    return it;
}

CompactedMinimizerIndex::iterator CompactedMinimizerIndex::end() {

    return iterator(this);
}

CompactedMinimizerIndex::const_iterator CompactedMinimizerIndex::end() const {

    return const_iterator(this);
}*/