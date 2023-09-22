#include "MinimizerIndex.hpp"

MinimizerIndex::MinimizerIndex() :  table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                    size_(0), pop(0), sum_psl(0), max_psl(1), max_ratio_occupancy(0.95),
                                    M_u64(0)  {

    init_tables(1024);
}

MinimizerIndex::MinimizerIndex(const size_t sz, const double ratio_occupancy) : table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                                                                size_(0), pop(0), sum_psl(0), max_psl(1), M_u64(0),
                                                                                max_ratio_occupancy(ratio_occupancy) {

    if (sz == 0) init_tables(1024);
    else {

        const size_t sz_with_empty = static_cast<size_t>((1.0 + (1.0 - ratio_occupancy)) * sz);

        init_tables(sz_with_empty);
    }
}

MinimizerIndex::MinimizerIndex(const MinimizerIndex& o) :   size_(o.size_), pop(o.pop), sum_psl(o.sum_psl), max_psl(o.max_psl),
                                                            M_u64(o.M_u64), max_ratio_occupancy(o.max_ratio_occupancy) {

    table_keys = new Minimizer[size_];
    table_tinyv = new packed_tiny_vector[size_];
    table_tinyv_sz = new uint8_t[size_];

    std::copy(o.table_keys, o.table_keys + size_, table_keys);

    for (size_t i = 0; i < size_; ++i){

        table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;
        table_tinyv[i].copy(table_tinyv_sz[i], o.table_tinyv[i], o.table_tinyv_sz[i]);
    }
}

MinimizerIndex::MinimizerIndex(MinimizerIndex&& o) :    table_keys(o.table_keys), table_tinyv(o.table_tinyv), table_tinyv_sz(o.table_tinyv_sz),
                                                        size_(o.size_), pop(o.pop), sum_psl(o.sum_psl), max_psl(o.max_psl),
                                                        M_u64(o.M_u64), max_ratio_occupancy(o.max_ratio_occupancy) {

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
        sum_psl = o.sum_psl;
        max_psl = o.max_psl;
        max_ratio_occupancy = o.max_ratio_occupancy;
        M_u64 = o.M_u64;

        table_keys = new Minimizer[size_];
        table_tinyv = new packed_tiny_vector[size_];
        table_tinyv_sz = new uint8_t[size_];

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
        sum_psl = o.sum_psl;
        max_psl = o.max_psl;
        max_ratio_occupancy = o.max_ratio_occupancy;
        M_u64 = o.M_u64;

        table_keys = o.table_keys;
        table_tinyv = o.table_tinyv;
        table_tinyv_sz = o.table_tinyv_sz;

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
}

MinimizerIndex::iterator MinimizerIndex::find(const Minimizer& key) {

    const size_t end_table = size_-1;
    const size_t mean_psl = get_mean_psl();

    size_t psl = 0;
    size_t h = fastmod::fastmod_u64(key.hash(), M_u64, size_);

    if (mean_psl <= 2) {

        while ((psl != max_psl) && !table_keys[h].isEmpty() && (table_keys[h] != key)) {

            h = (h+1) & (static_cast<size_t>(h == end_table) - 1);
            ++psl;
        }

        if ((psl != max_psl) && (table_keys[h] == key)) return iterator(this, h, psl);
    }
    else {

        size_t h_mean = fastmod::fastmod_u64(h + mean_psl, M_u64, size_);
        size_t h_inc = h_mean, h_dec = h_mean;

        bool has_empty_key = false;

        size_t i = 0;

        // Check all elements located at positions mean + i and mean -i. Stop if empty key encountered. Stop if minimum key is encountered.
        for (; !has_empty_key && (i <= mean_psl); ++i) {

            if (table_keys[h_dec] == key) return iterator(this, h_dec, mean_psl - i);

            has_empty_key = table_keys[h_dec].isEmpty() || table_keys[h_inc].isEmpty();

            if (!has_empty_key && (table_keys[h_inc] == key)) return iterator(this, h_inc, mean_psl + i);

            h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        // Only check remaining acending positions if neither the empty key nor the minimum key were encountered
        for (; !has_empty_key && (mean_psl + i <= max_psl); ++i) {

            if (table_keys[h_inc] == key) return iterator(this, h_inc, mean_psl + i);

            has_empty_key = table_keys[h_inc].isEmpty();
            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        if (has_empty_key) { // Only check remaining descending positions if empty key was previously encountered but minimum key was not encountered

            for (; (i <= mean_psl); ++i) {

                if (table_keys[h_dec] == key) return iterator(this, h_dec, mean_psl - i);

                h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            }
        }
    }

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const Minimizer& key) const {

    const size_t end_table = size_-1;
    const size_t mean_psl = get_mean_psl();

    size_t h = fastmod::fastmod_u64(key.hash(), M_u64, size_);

    if (mean_psl <= 2) {

        size_t psl = 0;

        while ((psl != max_psl) && !table_keys[h].isEmpty() && (table_keys[h] != key)) {

            h = (h+1) & (static_cast<size_t>(h==end_table)-1);
            ++psl;
        }

        if ((psl != max_psl) && (table_keys[h] == key)) return const_iterator(this, h, psl);
    }
    else {

        //size_t h_mean = (h + mean_psl) % size_;
        size_t h_mean = fastmod::fastmod_u64(h + mean_psl, M_u64, size_);
        size_t h_inc = h_mean, h_dec = h_mean;

        bool has_empty_key = false;

        size_t i = 0;

        // Check all elements located at positions mean + i and mean -i. Stop if empty key encountered. Stop if minimum key is encountered.
        for (; !has_empty_key && (i <= mean_psl); ++i) {

            if (table_keys[h_dec] == key) return const_iterator(this, h_dec, mean_psl - i);

            has_empty_key = table_keys[h_dec].isEmpty() || table_keys[h_inc].isEmpty();

            if (!has_empty_key && (table_keys[h_inc] == key)) return const_iterator(this, h_inc, mean_psl + i);

            h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        // Only check remaining acending positions if neither the empty key nor the minimum key were encountered
        for (; !has_empty_key && (mean_psl + i <= max_psl); ++i) {

            if (table_keys[h_inc] == key) return const_iterator(this, h_inc, mean_psl + i);

            has_empty_key = table_keys[h_inc].isEmpty();
            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        // Only check remaining descending positions if empty key was previously encountered but minimum key was not encountered
        if (has_empty_key) {

            for (; (i <= mean_psl); ++i) {

                if (table_keys[h_dec] == key) return const_iterator(this, h_dec, mean_psl - i);

                h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            }
        }
    }

    return const_iterator(this);
}

MinimizerIndex::iterator MinimizerIndex::find(const size_t h) {

    if ((h < size_) && !table_keys[h].isEmpty()) return iterator(this, h);

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const size_t h) const {

    if ((h < size_) && !table_keys[h].isEmpty()) return const_iterator(this, h);

    return const_iterator(this);
}

size_t MinimizerIndex::erase(const_iterator it) {

    if (it == end()) return 0;

    const size_t end_table = size_-1;

    table_keys[it.h].set_empty();
    table_tinyv[it.h].destruct(table_tinyv_sz[it.h]);
    table_tinyv_sz[it.h] = packed_tiny_vector::FLAG_EMPTY;

    --pop;
    
    if (it.psl != 0xffffffffffffffffULL) sum_psl -= it.psl;
    else {

        //const size_t h = table_keys[it.h].hash() % size_;
        const size_t h = fastmod::fastmod_u64(table_keys[it.h].hash(), M_u64, size_);

        sum_psl -= (size_ - h + it.h) & (static_cast<size_t>(it.h >= h) - 1);
        sum_psl -= (it.h - h) & (static_cast<size_t>(it.h < h) - 1);
    }

    // Robin-hood hashing
    // Push the tombstone further away if subsequent keys can be closer to where they were supposed to be
    size_t i = 0;
    size_t j1 = it.h;
    size_t j2 = (it.h + 1) & (static_cast<size_t>(it.h == end_table) - 1);

    while ((i != size_) && !table_keys[j2].isEmpty() && (fastmod::fastmod_u64(table_keys[j2].hash(), M_u64, size_) != j2)) {

        swap(j1, j2);

        j1 = j2;
        j2 = (j1 + 1) & (static_cast<size_t>(j1 == end_table) - 1);

        --sum_psl; 
        ++i;
    }

    return 1;
}

// Insert with Robin Hood hashing
pair<MinimizerIndex::iterator, bool> MinimizerIndex::insert(const Minimizer& key, const packed_tiny_vector& ptv, const uint8_t& flag) {

    if (pop >= static_cast<size_t>(size_ * max_ratio_occupancy)) {

        size_t resize = 1.2 * size_;

        while (pop >= static_cast<size_t>(resize * max_ratio_occupancy)) resize *= 1.2;

        reserve(resize);
    }

    const size_t end_table = size_-1;

    bool has_rich_psl = false, cascade_ins = false;

    size_t h = fastmod::fastmod_u64(key.hash(), M_u64, size_);
    size_t h_rich_psl_ins = 0;
    size_t psl_ins_key = 0, psl_rich_key = 0, psl_curr_key = 0;

    pair<MinimizerIndex::iterator, bool> it_ret;

    Minimizer l_key = key;

    packed_tiny_vector l_ptv(ptv, flag);

    uint8_t l_flag = flag;

    ++pop; // Pre-emptively increase population

    while (true) {

        if (table_keys[h].isEmpty() || (has_rich_psl && (cascade_ins || (psl_ins_key >= max_psl)))) {

            if (has_rich_psl) {

                packed_tiny_vector l_ptv_swap;
                uint8_t l_flag_swap = packed_tiny_vector::FLAG_EMPTY;

                h = h_rich_psl_ins;

                std::swap(table_keys[h], l_key);

                l_ptv_swap.move(l_flag_swap, move(table_tinyv[h]), move(table_tinyv_sz[h]));
                table_tinyv[h].move(table_tinyv_sz[h], move(l_ptv), move(l_flag));
                l_ptv.move(l_flag, move(l_ptv_swap), move(l_flag_swap));

                if (!cascade_ins) it_ret = {iterator(this, h, psl_rich_key), true};

                max_psl = max(max_psl, psl_rich_key + 1);
                sum_psl -= psl_curr_key;
                sum_psl += psl_rich_key;

                psl_ins_key = psl_curr_key;
                has_rich_psl = false;
                cascade_ins = true;
            }
            else {

                table_keys[h] = l_key;
                table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;
                table_tinyv[h].move(table_tinyv_sz[h], move(l_ptv), move(l_flag));

                max_psl = max(max_psl, psl_ins_key + 1);
                sum_psl += psl_ins_key;

                if (!cascade_ins) it_ret = {iterator(this, h, psl_ins_key), true};

                return it_ret;
            }
        }
        else if (table_keys[h] == l_key) {

            --pop; // Key already in there, pre-emptive population increase was not necessary

            return {iterator(this, h, psl_ins_key), false}; // Can only happen when inserting the input key
        }
        else if (!has_rich_psl) {

            const size_t h_curr = fastmod::fastmod_u64(table_keys[h].hash(), M_u64, size_);
            
            psl_curr_key = ((size_ - h_curr + h) & (static_cast<size_t>(h >= h_curr) - 1)) + ((h - h_curr) & (static_cast<size_t>(h < h_curr) - 1));

            if (psl_ins_key > psl_curr_key) {

                h_rich_psl_ins = h;
                psl_rich_key = psl_ins_key;
                has_rich_psl = true;
            }
        }

        h = (h + 1) & (static_cast<size_t>(h == end_table) - 1);
        ++psl_ins_key;
    }
}

void MinimizerIndex::recomputeMaxPSL(const size_t nb_threads) {

    max_psl = 1;

    if (pop != 0) {

        if (nb_threads <= 1){

            for (size_t i = 0; i != size_; ++i) {

                if (!table_keys[i].isEmpty()) {

                    const size_t h = fastmod::fastmod_u64(table_keys[i].hash(), M_u64, size_);
                    const size_t psl = ((size_ - h + i) & (static_cast<size_t>(i >= h) - 1)) + ((i - h) & (static_cast<size_t>(i < h) - 1));

                    max_psl = max(max_psl, psl + 1);
                }
            }
        }
        else {

            const size_t chunk_per_thread = (size_ + nb_threads - 1) / nb_threads;

            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mtx_max_psl;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        const size_t chunk_start = t * chunk_per_thread;
                        const size_t chunk_end = min(((t+1) * chunk_per_thread), size_);

                        size_t l_max_psl = 1;

                        for (size_t i = chunk_start; i < chunk_end; ++i) {

                            if (!table_keys[i].isEmpty()) {

                                const size_t h = fastmod::fastmod_u64(table_keys[i].hash(), M_u64, size_);
                                const size_t psl = ((size_ - h + i) & (static_cast<size_t>(i >= h) - 1)) + ((i - h) & (static_cast<size_t>(i < h) - 1));

                                l_max_psl = max(l_max_psl, psl + 1);
                            }
                        }

                        {
                            unique_lock<mutex> lock(mtx_max_psl);

                            max_psl = max(max_psl, l_max_psl);
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }
    }
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
    sum_psl = 0;
    max_psl = 1;
    M_u64 = 0;
}

void MinimizerIndex::init_tables(const size_t sz) {

    clear_tables();

    Minimizer empty_key;

    pop = 0;
    size_ = sz;
    M_u64 = fastmod::computeM_u64(size_);

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

    size_ = sz;
    pop = 0;
    sum_psl = 0;
    max_psl = 1;
    M_u64 = fastmod::computeM_u64(size_);

    table_keys = new Minimizer[size_];
    table_tinyv = new packed_tiny_vector[size_];
    table_tinyv_sz = new uint8_t[size_];

    empty_key.set_empty();

    std::fill(table_keys, table_keys + size_, empty_key);

    memset(table_tinyv_sz, packed_tiny_vector::FLAG_EMPTY, size_ * sizeof(uint8_t));

    for (size_t i = 0; i < old_size_; ++i) {

        if (!old_table_keys[i].isEmpty()){

            insert(old_table_keys[i], old_table_tinyv[i], old_table_tinyv_sz[i]);

            old_table_tinyv[i].destruct(old_table_tinyv_sz[i]);
        }
    }

    delete[] old_table_keys;
    delete[] old_table_tinyv;
    delete[] old_table_tinyv_sz;
}

void MinimizerIndex::swap(const size_t i, const size_t j) {

    uint8_t ptv_sz = packed_tiny_vector::FLAG_EMPTY;

    packed_tiny_vector ptv;

    ptv.move(ptv_sz, move(table_tinyv[i]), move(table_tinyv_sz[i]));
    table_tinyv[i].move(table_tinyv_sz[i], move(table_tinyv[j]), move(table_tinyv_sz[j]));
    table_tinyv[j].move(table_tinyv_sz[j], move(ptv), move(ptv_sz));

    std::swap(table_keys[i], table_keys[j]);
}


/*MinimizerIndex::MinimizerIndex() :  table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                    size_(0), pop(0), sum_psl(0), max_psl(1), max_ratio_occupancy(0.95)  {

    init_tables(1024);
}

MinimizerIndex::MinimizerIndex(const size_t sz, const double ratio_occupancy) : table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                                                                size_(0), pop(0), sum_psl(0), max_psl(1), max_ratio_occupancy(ratio_occupancy) {

    if (sz == 0) init_tables(1024);
    else {

        const size_t sz_with_empty = static_cast<size_t>((1.0 + (1.0 - ratio_occupancy)) * sz);

        init_tables(sz_with_empty);
    }
}

MinimizerIndex::MinimizerIndex(const MinimizerIndex& o) :   size_(o.size_), pop(o.pop), sum_psl(o.sum_psl), max_psl(o.max_psl), max_ratio_occupancy(o.max_ratio_occupancy) {

    table_keys = new Minimizer[size_];
    table_tinyv = new packed_tiny_vector[size_];
    table_tinyv_sz = new uint8_t[size_];

    std::copy(o.table_keys, o.table_keys + size_, table_keys);

    for (size_t i = 0; i < size_; ++i){

        table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;
        table_tinyv[i].copy(table_tinyv_sz[i], o.table_tinyv[i], o.table_tinyv_sz[i]);
    }
}

MinimizerIndex::MinimizerIndex(MinimizerIndex&& o){

    size_ = o.size_;
    pop = o.pop;
    sum_psl = o.sum_psl;
    max_psl = o.max_psl;
    max_ratio_occupancy = o.max_ratio_occupancy;

    table_keys = o.table_keys;
    table_tinyv = o.table_tinyv;
    table_tinyv_sz = o.table_tinyv_sz;

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
        sum_psl = o.sum_psl;
        max_psl = o.max_psl;
        max_ratio_occupancy = o.max_ratio_occupancy;

        table_keys = new Minimizer[size_];
        table_tinyv = new packed_tiny_vector[size_];
        table_tinyv_sz = new uint8_t[size_];

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
        sum_psl = o.sum_psl;
        max_psl = o.max_psl;
        max_ratio_occupancy = o.max_ratio_occupancy;

        table_keys = o.table_keys;
        table_tinyv = o.table_tinyv;
        table_tinyv_sz = o.table_tinyv_sz;

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
}

MinimizerIndex::iterator MinimizerIndex::find(const Minimizer& key) {

    const size_t end_table = size_-1;
    const size_t mean_psl = get_mean_psl();

    size_t h = key.hash() % size_;

    if (mean_psl <= 2) {

        size_t psl = 0;

        while ((psl != max_psl) && !table_keys[h].isEmpty()) {

            if (table_keys[h] == key) return iterator(this, h, psl);

            h = (h+1) & (static_cast<size_t>(h==end_table)-1);
            ++psl;
        }
    }
    else {

        size_t h_mean = (h + mean_psl) % size_;
        size_t h_inc = h_mean, h_dec = h_mean;

        bool has_empty_key = false;

        size_t i = 0;

        // Check all elements located at positions mean + i and mean -i. Stop if empty key encountered. Stop if minimum key is encountered.
        for (; !has_empty_key && (i <= mean_psl); ++i) {

            if (table_keys[h_dec] == key) return iterator(this, h_dec, mean_psl - i);

            if (table_keys[h_dec].isEmpty() || table_keys[h_inc].isEmpty()) has_empty_key = true;
            else if (table_keys[h_inc] == key) return iterator(this, h_inc, mean_psl + i);

            h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        // Only check remaining acending positions if neither the empty key nor the minimum key were encountered
        for (; !has_empty_key && (mean_psl + i <= max_psl); ++i) {

            if (table_keys[h_inc] == key) return iterator(this, h_inc, mean_psl + i);
            if (table_keys[h_inc].isEmpty()) has_empty_key = true;

            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        if (has_empty_key) { // Only check remaining descending positions if empty key was previously encountered but minimum key was not encountered

            for (; (i <= mean_psl); ++i) {

                if (table_keys[h_dec] == key) return iterator(this, h_dec, mean_psl - i);

                h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            }
        }
    }

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const Minimizer& key) const {

    const size_t end_table = size_-1;
    const size_t mean_psl = get_mean_psl();

    size_t h = key.hash() % size_;

    if (mean_psl <= 2) {

        size_t psl = 0;

        while ((psl != max_psl) && !table_keys[h].isEmpty()) {

            if (table_keys[h] == key) return const_iterator(this, h, psl);

            h = (h+1) & (static_cast<size_t>(h==end_table)-1);
            ++psl;
        }
    }
    else {

        size_t h_mean = (h + mean_psl) % size_;
        size_t h_inc = h_mean, h_dec = h_mean;

        bool has_empty_key = false;

        size_t i = 0;

        // Check all elements located at positions mean + i and mean -i. Stop if empty key encountered. Stop if minimum key is encountered.
        for (; !has_empty_key && (i <= mean_psl); ++i) {

            if (table_keys[h_dec] == key) return const_iterator(this, h_dec, mean_psl - i);

            if (table_keys[h_dec].isEmpty() || table_keys[h_inc].isEmpty()) has_empty_key = true;
            else if (table_keys[h_inc] == key) return const_iterator(this, h_inc, mean_psl + i);

            h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        // Only check remaining acending positions if neither the empty key nor the minimum key were encountered
        for (; !has_empty_key && (mean_psl + i <= max_psl); ++i) {

            if (table_keys[h_inc] == key) return const_iterator(this, h_inc, mean_psl + i);
            if (table_keys[h_inc].isEmpty()) has_empty_key = true;

            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        // Only check remaining descending positions if empty key was previously encountered but minimum key was not encountered
        if (has_empty_key) {

            for (; (i <= mean_psl); ++i) {

                if (table_keys[h_dec] == key) return const_iterator(this, h_dec, mean_psl - i);

                h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            }
        }
    }

    return const_iterator(this);
}

MinimizerIndex::iterator MinimizerIndex::find(const size_t h) {

    if ((h < size_) && !table_keys[h].isEmpty()) return iterator(this, h);

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const size_t h) const {

    if ((h < size_) && !table_keys[h].isEmpty()) return const_iterator(this, h);

    return const_iterator(this);
}

size_t MinimizerIndex::erase(const_iterator it) {

    if (it == end()) return 0;

    const size_t end_table = size_-1;

    table_keys[it.h].set_empty();
    table_tinyv[it.h].destruct(table_tinyv_sz[it.h]);
    table_tinyv_sz[it.h] = packed_tiny_vector::FLAG_EMPTY;

    --pop;
    
    if (it.psl != 0xffffffffffffffffULL) sum_psl -= it.psl;
    else {

        const size_t h = table_keys[it.h].hash() % size_;

        sum_psl -= (it.h < h) ? (size_ - h + it.h) : (it.h - h);
    }

    // Robin-hood hashing
    // Push the tombstone further away if subsequent keys can be closer to where they were supposed to be
    for (size_t i = 0, j1 = it.h; i != size_; ++i) {

        const size_t j2 = (j1 + 1) & (static_cast<size_t>(j1==end_table)-1);

        if (table_keys[j2].isEmpty() || ((table_keys[j2].hash() % size_) == j2)) break;

        swap(j1, j2);

        j1 = j2;
        --sum_psl;
    }

    return 1;
}

// Insert with Robin Hood hashing
pair<MinimizerIndex::iterator, bool> MinimizerIndex::insert(const Minimizer& key, const packed_tiny_vector& ptv, const uint8_t& flag) {

    if (pop >= static_cast<size_t>(size_ * max_ratio_occupancy)) {

        size_t resize = 1.2 * size_;

        while (pop >= static_cast<size_t>(resize * max_ratio_occupancy)) resize *= 1.2;

        reserve(resize);
    }

    const size_t end_table = size_-1;

    bool has_rich_psl = false, cascade_ins = false;

    size_t h = key.hash() % size_;
    size_t h_rich_psl_ins = 0;
    size_t psl_ins_key = 0, psl_rich_key = 0, psl_curr_key = 0;

    pair<MinimizerIndex::iterator, bool> it_ret;

    Minimizer l_key = key;

    packed_tiny_vector l_ptv(ptv, flag);

    uint8_t l_flag = flag;

    ++pop; // Pre-emptively increase population

    while (true) {

        if (table_keys[h].isEmpty() || (has_rich_psl && (cascade_ins || (psl_ins_key >= max_psl)))) {

            if (has_rich_psl) {

                packed_tiny_vector l_ptv_swap;
                uint8_t l_flag_swap = packed_tiny_vector::FLAG_EMPTY;

                h = h_rich_psl_ins;

                std::swap(table_keys[h], l_key);

                l_ptv_swap.move(l_flag_swap, move(table_tinyv[h]), move(table_tinyv_sz[h]));
                table_tinyv[h].move(table_tinyv_sz[h], move(l_ptv), move(l_flag));
                l_ptv.move(l_flag, move(l_ptv_swap), move(l_flag_swap));

                if (!cascade_ins) it_ret = {iterator(this, h, psl_rich_key), true};

                max_psl = max(max_psl, psl_rich_key + 1);
                sum_psl -= psl_curr_key;
                sum_psl += psl_rich_key;

                psl_ins_key = psl_curr_key;
                has_rich_psl = false;
                cascade_ins = true;
            }
            else {

                table_keys[h] = l_key;
                table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;
                table_tinyv[h].move(table_tinyv_sz[h], move(l_ptv), move(l_flag));

                max_psl = max(max_psl, psl_ins_key + 1);
                sum_psl += psl_ins_key;

                if (!cascade_ins) it_ret = {iterator(this, h, psl_ins_key), true};

                return it_ret;
            }
        }
        else if (table_keys[h] == l_key) {

            --pop; // Key already in there, pre-emptive population increase was not necessary

            return {iterator(this, h, psl_ins_key), false}; // Can only happen when inserting the input key
        }
        else if (!has_rich_psl) {

            const size_t h_curr = table_keys[h].hash() % size_;
            
            psl_curr_key = (h < h_curr) ? (size_ - h_curr + h) : (h - h_curr);

            if (psl_ins_key > psl_curr_key) {

                h_rich_psl_ins = h;
                psl_rich_key = psl_ins_key;
                has_rich_psl = true;
            }
        }

        h = (h + 1) & (static_cast<size_t>(h == end_table) - 1);
        ++psl_ins_key;
    }
}

void MinimizerIndex::recomputeMaxPSL(const size_t nb_threads) {

    max_psl = 1;

    if (pop != 0) {

        if (nb_threads <= 1){

            for (size_t i = 0; i != size_; ++i) {

                if (!table_keys[i].isEmpty()) {

                    const size_t h = table_keys[i].hash() % size_;
                    const size_t psl = (i < h) ? (size_ - h + i) : (i - h);

                    max_psl = max(max_psl, psl + 1);
                }
            }
        }
        else {

            const size_t chunk_per_thread = (size_ + nb_threads - 1) / nb_threads;

            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mtx_max_psl;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        const size_t chunk_start = t * chunk_per_thread;
                        const size_t chunk_end = min(((t+1) * chunk_per_thread), size_);

                        size_t l_max_psl = 1;

                        for (size_t i = chunk_start; i < chunk_end; ++i) {

                            if (!table_keys[i].isEmpty()) {

                                const size_t h = table_keys[i].hash() % size_;
                                const size_t psl = (i < h) ? (size_ - h + i) : (i - h);

                                l_max_psl = max(l_max_psl, psl + 1);
                            }
                        }

                        {
                            unique_lock<mutex> lock(mtx_max_psl);

                            max_psl = max(max_psl, l_max_psl);
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }
    }
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
    sum_psl = 0;
    max_psl = 1;
}

void MinimizerIndex::init_tables(const size_t sz) {

    clear_tables();

    Minimizer empty_key;

    pop = 0;
    size_ = sz;

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

    size_ = sz;
    pop = 0;
    sum_psl = 0;
    max_psl = 1;

    table_keys = new Minimizer[size_];
    table_tinyv = new packed_tiny_vector[size_];
    table_tinyv_sz = new uint8_t[size_];

    empty_key.set_empty();

    std::fill(table_keys, table_keys + size_, empty_key);

    memset(table_tinyv_sz, packed_tiny_vector::FLAG_EMPTY, size_ * sizeof(uint8_t));

    for (size_t i = 0; i < old_size_; ++i) {

        if (!old_table_keys[i].isEmpty()){

            insert(old_table_keys[i], old_table_tinyv[i], old_table_tinyv_sz[i]);

            old_table_tinyv[i].destruct(old_table_tinyv_sz[i]);
        }
    }

    delete[] old_table_keys;
    delete[] old_table_tinyv;
    delete[] old_table_tinyv_sz;
}

void MinimizerIndex::swap(const size_t i, const size_t j) {

    uint8_t ptv_sz = packed_tiny_vector::FLAG_EMPTY;
    packed_tiny_vector ptv;

    ptv.move(ptv_sz, move(table_tinyv[i]), move(table_tinyv_sz[i]));
    table_tinyv[i].move(table_tinyv_sz[i], move(table_tinyv[j]), move(table_tinyv_sz[j]));
    table_tinyv[j].move(table_tinyv_sz[j], move(ptv), move(ptv_sz));

    std::swap(table_keys[i], table_keys[j]);
}*/

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