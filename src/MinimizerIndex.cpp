#include "MinimizerIndex.hpp"

MinimizerIndex::MinimizerIndex() :  table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                    size_(0), pop(0), sum_psl(0), max_psl(1), max_ratio_occupancy(0.95)  {

    //init_tables(max(static_cast<size_t>(1024), lck_block_sz));
    init_tables(1024);
}

MinimizerIndex::MinimizerIndex(const size_t sz, const double ratio_occupancy) : table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                                                                size_(0), pop(0), sum_psl(0), max_psl(1), max_ratio_occupancy(ratio_occupancy) {

    if (sz == 0) /*init_tables(lck_block_sz);*/ init_tables(1024);
    else {

        const size_t sz_with_empty = static_cast<size_t>((1.0 + (1.0 - ratio_occupancy)) * sz);

        //size_t rdnup_sz = rndup(sz);

        //while (rdnup_sz < sz_with_empty) rdnup_sz <<= 1;

        //init_tables(max(rdnup_sz, lck_block_sz));

        init_tables(sz_with_empty);
    }
}

MinimizerIndex::MinimizerIndex(const MinimizerIndex& o) :   size_(o.size_), pop(o.pop), sum_psl(o.sum_psl), max_psl(o.max_psl), max_ratio_occupancy(o.max_ratio_occupancy) {

    table_keys = new Minimizer[size_];
    table_tinyv = new packed_tiny_vector[size_];
    table_tinyv_sz = new uint8_t[size_];

    //lck_min = vector<SpinLock>(o.lck_min.size());

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

    //lck_min = vector<SpinLock>(o.lck_min.size());

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

        //lck_min = vector<SpinLock>(o.lck_min.size());

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

        //lck_min = vector<SpinLock>(o.lck_min.size());

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

    //lck_min.clear();
    //lck_edit_table.release_all();
}

/*MinimizerIndex::iterator MinimizerIndex::find(const Minimizer& key) {

    const size_t end_table = size_-1;
    const size_t mean_psl = get_mean_psl();

    size_t h = key.hash() % size_;
    size_t psl = 0;

    if (mean_psl == 0) {

        while (psl != max_psl) {

            if (table_keys[h].isEmpty() || (table_keys[h] == key)) break;

            h = (h+1) & (static_cast<size_t>(h==end_table)-1);
            ++psl;
        }

        if ((psl != max_psl) && (table_keys[h] == key)) return iterator(this, h);
    }
    else {

        size_t h_mean = (h + mean_psl) % size_;
        size_t h_inc = h_mean, h_dec = h_mean;

        {
            if (table_keys[h_inc] == key) return iterator(this, h_inc);

            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
            h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & end_table);
        }

        for (size_t i = 1; mean_psl + i <= max_psl; ++i) {

            if (table_keys[h_inc] == key) return iterator(this, h_inc);
            if (table_keys[h_dec] == key) return iterator(this, h_dec);

            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
            h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & end_table);
        }
    }

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const Minimizer& key) const {

    const size_t end_table = size_-1;
    const size_t mean_psl = get_mean_psl();

    size_t h = key.hash() % size_;
    size_t psl = 0;

    if (mean_psl == 0) {

        while (psl != max_psl) {

            if (table_keys[h].isEmpty() || (table_keys[h] == key)) break;

            h = (h+1) & (static_cast<size_t>(h==end_table)-1);
            ++psl;
        }

        if ((psl != max_psl) && (table_keys[h] == key)) return const_iterator(this, h);
    }
    else {

        size_t h_mean = (h + mean_psl) % size_;
        size_t h_inc = h_mean, h_dec = h_mean;

        {
            if (table_keys[h_inc] == key) return const_iterator(this, h_inc);

            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
            h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & end_table);
        }

        for (size_t i = 1; mean_psl + i <= max_psl; ++i) {

            if (table_keys[h_inc] == key) return const_iterator(this, h_inc);
            if (table_keys[h_dec] == key) return const_iterator(this, h_dec);

            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
            h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & end_table);
        }
    }

    return const_iterator(this);
}*/

MinimizerIndex::iterator MinimizerIndex::find(const Minimizer& key) {

    const size_t end_table = size_-1;
    const size_t mean_psl = get_mean_psl();

    size_t h = key.hash() % size_;
    size_t psl = 0;

    if (mean_psl <= 2) {

        while (psl != max_psl) {

            if (table_keys[h].isEmpty() || (table_keys[h] == key)) break;

            h = (h+1) & (static_cast<size_t>(h==end_table)-1);
            ++psl;
        }

        if ((psl != max_psl) && (table_keys[h] == key)) return iterator(this, h);
    }
    else {

        size_t h_mean = (h + mean_psl) % size_;
        size_t h_inc = h_mean, h_dec = h_mean;

        bool has_empty_key = false;

        size_t i = 0;

        // Check all elements located at positions mean + i and mean -i. Stop if empty key encountered. Stop if minimum key is encountered.
        for (; !has_empty_key && (i <= mean_psl); ++i) {

            if (table_keys[h_dec] == key) return iterator(this, h_dec);

            if (table_keys[h_dec].isEmpty() || table_keys[h_inc].isEmpty()) has_empty_key = true;
            else if (table_keys[h_inc] == key) return iterator(this, h_inc);

            h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        // Only check remaining acending positions if neither the empty key nor the minimum key were encountered
        for (; !has_empty_key && (mean_psl + i <= max_psl); ++i) {

            if (table_keys[h_inc] == key) return iterator(this, h_inc);
            if (table_keys[h_inc].isEmpty()) has_empty_key = true;

            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        if (has_empty_key) { // Only check remaining descending positions if empty key was previously encountered but minimum key was not encountered

            for (; (i <= mean_psl); ++i) {

                if (table_keys[h_dec] == key) return iterator(this, h_dec);

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

        while (psl != max_psl) {

            if (table_keys[h].isEmpty() || (table_keys[h] == key)) break;

            h = (h+1) & (static_cast<size_t>(h==end_table)-1);
            ++psl;
        }

        if ((psl != max_psl) && (table_keys[h] == key)) return const_iterator(this, h);
    }
    else {

        size_t h_mean = (h + mean_psl) % size_;
        size_t h_inc = h_mean, h_dec = h_mean;

        bool has_empty_key = false;

        size_t i = 0;

        // Check all elements located at positions mean + i and mean -i. Stop if empty key encountered. Stop if minimum key is encountered.
        for (; !has_empty_key && (i <= mean_psl); ++i) {

            if (table_keys[h_dec] == key) return const_iterator(this, h_dec);

            if (table_keys[h_dec].isEmpty() || table_keys[h_inc].isEmpty()) has_empty_key = true;
            else if (table_keys[h_inc] == key) return const_iterator(this, h_inc);

            h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        // Only check remaining acending positions if neither the empty key nor the minimum key were encountered
        for (; !has_empty_key && (mean_psl + i <= max_psl); ++i) {

            if (table_keys[h_inc] == key) return const_iterator(this, h_inc);
            if (table_keys[h_inc].isEmpty()) has_empty_key = true;

            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        // Only check remaining descending positions if empty key was previously encountered but minimum key was not encountered
        if (has_empty_key) {

            for (; (i <= mean_psl); ++i) {

                if (table_keys[h_dec] == key) return const_iterator(this, h_dec);

                h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            }
        }
    }

    return const_iterator(this);
}

/*MinimizerIndex::const_iterator MinimizerIndex::find(const Minimizer& key) const {

    const size_t end_table = size_-1;
    const size_t mean_psl = get_mean_psl();

    const size_t h = key.hash() % size_;

    size_t h1 = h;
    size_t h2 = h;

    bool found_h1 = false;
    bool found_h2 = false;

    bool has_empty_key = false;
    size_t h_dec_empty = 0, h_inc_empty = 0;
    size_t h_mean = (h + mean_psl) % size_;
    size_t h_inc = h_mean, h_dec = h_mean;

    {
        size_t psl = 0;

        while (psl != max_psl) {

            if (table_keys[h1].isEmpty() || (table_keys[h1] == key)) break;

            h1 = (h1+1) & (static_cast<size_t>(h1==end_table)-1);
            ++psl;
        }

        if ((psl != max_psl) && (table_keys[h1] == key)) found_h1 = true;
    }

    {

        size_t i = 0;

        // Check all elements located at positions mean + i and mean -i. Stop if empty key encountered.
        for (; !has_empty_key && (mean_psl + i <= max_psl); ++i) {

            if (table_keys[h_dec] == key) {

                found_h2 = true;
                h2 = h_dec;

                break;
            }

            if (table_keys[h_dec].isEmpty() || table_keys[h_inc].isEmpty()) {

                has_empty_key = true;
                h_dec_empty = h_dec;
                h_inc_empty = h_inc;
            }
            else if (table_keys[h_inc] == key) {

                found_h2 = true;
                h2 = h_inc;

                break;
            }

            h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
        }

        // Only check remaining descending positions if empty key was previously encountered
        if (has_empty_key) {

            for (; (mean_psl + i <= max_psl); ++i) {

                if (table_keys[h_dec] == key) {

                    found_h2 = true;
                    h2 = h_dec;

                    break;
                }

                h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            }
        }
    }

    if (found_h1 != found_h2){

        cout << found_h1 << "\t" << found_h2 << endl;
        cout << h << "\t" << h1 << "\t" << h2 << endl;
        cout << h_mean << "\t" << h_dec << "\t" << h_inc << endl;
        cout << has_empty_key << "\t" << h_dec_empty << "\t" << h_inc_empty << endl;

        exit(1);
    }

    if (found_h1) return const_iterator(this, h1);

    return const_iterator(this);
}*/

MinimizerIndex::iterator MinimizerIndex::find(const size_t h) {

    if ((h < size_) && !table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return iterator(this, h);

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const size_t h) const {

    if ((h < size_) && !table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return const_iterator(this, h);

    return const_iterator(this);
}

void MinimizerIndex::erase(const_iterator it) {

    if (it == end()) return;

    const size_t end_table = size_-1;
    const size_t h = table_keys[it.h].hash() % size_;
    const size_t psl = (it.h < h) ? (size_ - h + it.h) : (it.h - h);

    table_keys[it.h].set_deleted();
    table_tinyv[it.h].destruct(table_tinyv_sz[it.h]);
    table_tinyv_sz[it.h] = packed_tiny_vector::FLAG_EMPTY;

    --pop;
    sum_psl -= psl;

    // Robin-hood hashing
    // Push the tombstone further away if subsequent keys can be closer to where they were supposed to be
    for (size_t i = 0, j1 = it.h; i != size_; ++i) {

        const size_t j2 = (j1 + 1) & (static_cast<size_t>(j1==end_table)-1);
        const bool is_del = table_keys[j2].isDeleted();

        if (table_keys[j2].isEmpty() || (!is_del && ((table_keys[j2].hash() % size_) == j2))) break;

        swap(j1, j2);

        j1 = j2;
        sum_psl -= static_cast<size_t>(!is_del);
    }
}

size_t MinimizerIndex::erase(const Minimizer& minz) {

    const size_t end_table = size_-1;
    const size_t oldpop = pop;

    size_t h = minz.hash() % size_;
    size_t psl = 0;

    while (psl != max_psl) {

        if (table_keys[h].isEmpty() || (table_keys[h] == minz)) break;

        h = (h+1) & (static_cast<size_t>(h==end_table)-1);
        ++psl;
    }

    if ((psl != max_psl) && (table_keys[h] == minz)){

        table_keys[h].set_deleted();
        table_tinyv[h].destruct(table_tinyv_sz[h]);
        table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;

        --pop;
        sum_psl -= psl;

        // Robin-hood hashing
        // Push the tombstone further away if subsequent keys can be closer to where they were supposed to be
        for (size_t i = 0, j1 = h; i != size_; ++i) {

            const size_t j2 = (j1 + 1) & (static_cast<size_t>(j1==end_table)-1);
            const bool is_del = table_keys[j2].isDeleted();

            if (table_keys[j2].isEmpty() || (!is_del && ((table_keys[j2].hash() % size_) == j2))) break;

            swap(j1, j2);

            j1 = j2;
            sum_psl -= static_cast<size_t>(!is_del);
        }
    }

    return oldpop - pop;
}

// Insert with Robin Hood hashing
pair<MinimizerIndex::iterator, bool> MinimizerIndex::insert(const Minimizer& key, const packed_tiny_vector& ptv, const uint8_t& flag) {

    if (pop >= static_cast<size_t>(size_ * max_ratio_occupancy)) {

        size_t resize = 1.2 * size_;

        while (pop >= static_cast<size_t>(resize * max_ratio_occupancy)) resize *= 1.2;

        reserve(resize);
    }

    const size_t end_table = size_-1;

    bool is_deleted = false, has_rich_psl = false;

    size_t h = key.hash() % size_;
    size_t h_del = 0, h_rich_psl_ins = 0, h_rich_psl0 = 0;
    size_t psl_ins_key = 0, psl_rich_key = 0, psl_del_key = 0, psl_curr_key = 0;

    pair<MinimizerIndex::iterator, bool> it_ret;

    Minimizer l_key = key;

    packed_tiny_vector l_ptv(ptv, flag);

    uint8_t l_flag = flag;

    ++pop; // Pre-emptively increase population

    while (true) {

        if (table_keys[h].isEmpty() || ((psl_ins_key > max_psl) && (is_deleted || has_rich_psl))) {

            if (has_rich_psl) {

                packed_tiny_vector l_ptv_swap;
                uint8_t l_flag_swap = packed_tiny_vector::FLAG_EMPTY;

                std::swap(table_keys[h_rich_psl_ins], l_key);

                l_ptv_swap.move(l_flag_swap, move(table_tinyv[h_rich_psl_ins]), move(table_tinyv_sz[h_rich_psl_ins]));
                table_tinyv[h_rich_psl_ins].move(table_tinyv_sz[h_rich_psl_ins], move(l_ptv), move(l_flag));
                l_ptv.move(l_flag, move(l_ptv_swap), move(l_flag_swap));

                if (table_keys[h_rich_psl_ins] == key) it_ret = {iterator(this, h_rich_psl_ins), true};

                max_psl = max(max_psl, psl_rich_key + 1);
                sum_psl -= psl_curr_key;
                sum_psl += psl_rich_key;

                h = h_rich_psl0;

                psl_ins_key = 0;
                is_deleted = false;
                has_rich_psl = false;

                --psl_ins_key;
                --h;
            }
            else {

                size_t psl = psl_ins_key;

                if (is_deleted) {

                    h = h_del;
                    psl = psl_del_key;
                }

                max_psl = max(max_psl, psl + 1);
                sum_psl += psl;

                table_keys[h] = l_key;
                table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;
                table_tinyv[h].move(table_tinyv_sz[h], move(l_ptv), move(l_flag));

                if (table_keys[h] == key) it_ret = {iterator(this, h), true};

                return it_ret;
            }
        }
        else if (table_keys[h] == l_key) {

            --pop; // Key already in there, pre-emptive population increase was not necessary

            return {iterator(this, h), false}; // Can only happen when inserting the input key
        }
        else if (!is_deleted && !has_rich_psl) {

            if (table_keys[h].isDeleted()) {

                if (l_key != key) {

                    table_keys[h] = l_key;
                    table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;
                    table_tinyv[h].move(table_tinyv_sz[h], move(l_ptv), move(l_flag));

                    max_psl = max(max_psl, psl_ins_key + 1);
                    sum_psl += psl_ins_key;

                    return it_ret;
                }

                h_del = h;
                psl_del_key = psl_ins_key;
                is_deleted = true;
            }
            else {

                const size_t h_curr = table_keys[h].hash() % size_;
                
                psl_curr_key = (h < h_curr) ? (size_ - h_curr + h) : (h - h_curr);

                if (psl_ins_key > psl_curr_key) {

                    h_rich_psl_ins = h;
                    h_rich_psl0 = h_curr;
                    psl_rich_key = psl_ins_key;
                    has_rich_psl = true;

                    if (l_key != key) psl_ins_key = max_psl;
                }
            }
        }

        h = (h+1) & (static_cast<size_t>(h==end_table)-1);
        ++psl_ins_key;
    }
}

void MinimizerIndex::recomputeMaxPSL(const size_t nb_threads) {

    max_psl = 1;

    if (pop != 0) {

        if (nb_threads <= 1){

            for (size_t i = 0; i != size_; ++i) {

                if (!table_keys[i].isDeleted() && !table_keys[i].isEmpty()) {

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

                            if (!table_keys[i].isDeleted() && !table_keys[i].isEmpty()) {

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


/*pair<MinimizerIndex::iterator, bool> MinimizerIndex::insert(const Minimizer& key, const packed_tiny_vector& ptv, const uint8_t& flag) {

    if (pop >= static_cast<size_t>(size_ * max_ratio_occupancy)) reserve(2 * size_); // if more than 80% full, resize

    const size_t end_table = size_-1;
    
    size_t h = key.hash() & end_table, h_del;

    bool is_deleted = false;

    while (true) {

        if (table_keys[h].isEmpty()) {

            h = ((static_cast<size_t>(is_deleted) - 1) & h) + ((static_cast<size_t>(!is_deleted) - 1) & h_del);

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
}*/

/*void MinimizerIndex::init_threads() {

    lck_min = vector<SpinLock>((size_ + lck_block_sz - 1) / lck_block_sz);

    pop_p = pop;
}

void MinimizerIndex::release_threads() {

    pop = pop_p;

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

    if (pop_p >= static_cast<size_t>(size_ * max_ratio_occupancy)) {

        lck_edit_table.release_reader();
        lck_edit_table.acquire_writer();

        reserve(2 * size_); // if more than 80% full, resize

        pop_p = pop;

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
                    else lck_min[id_block2].release();
                }
            }

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
}*/

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

    //if (!lck_min.empty()) lck_min = vector<SpinLock>((size_ + lck_block_sz - 1) / lck_block_sz);

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

void MinimizerIndex::swap(const size_t i, const size_t j) {

    uint8_t ptv_sz = packed_tiny_vector::FLAG_EMPTY;
    packed_tiny_vector ptv;

    ptv.move(ptv_sz, move(table_tinyv[i]), move(table_tinyv_sz[i]));
    table_tinyv[i].move(table_tinyv_sz[i], move(table_tinyv[j]), move(table_tinyv_sz[j]));
    table_tinyv[j].move(table_tinyv_sz[j], move(ptv), move(ptv_sz));

    std::swap(table_keys[i], table_keys[j]);
}

//const size_t MinimizerIndex::lck_block_sz = 64;
//const size_t MinimizerIndex::lck_block_div_shift = 6;

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