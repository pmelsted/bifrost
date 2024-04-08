#include "MinimizerIndex.hpp"

MinimizerIndex::MinimizerIndex() :  table_keys(nullptr), table_tinyv(nullptr),
                                    table_tinyv_sz(nullptr), table_outliers_psl(nullptr) {

    clear_tables();
}

MinimizerIndex::MinimizerIndex(const size_t sz, const double ratio_occupancy) : table_keys(nullptr), table_tinyv(nullptr),
                                                                                table_tinyv_sz(nullptr), table_outliers_psl(nullptr) {

    clear_tables();

    max_ratio_occupancy = ratio_occupancy;

    if (sz != 0) {

        const size_t sz_with_empty = static_cast<size_t>((1.0 + (1.0 - ratio_occupancy)) * sz);

        init_tables(max(sz_with_empty, static_cast<size_t>(BIFROST_MI_INIT_SZ)));
    }
}

MinimizerIndex::MinimizerIndex(const MinimizerIndex& o) :   table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                                            table_outliers_psl(nullptr), size_(o.size_), pop(o.pop),
                                                            sum_psl(o.sum_psl), max_psl(o.max_psl), std_psl(o.std_psl),
                                                            M_u64(o.M_u64), max_ratio_occupancy(o.max_ratio_occupancy) {

    if (size_ != 0) {

        table_keys = new Minimizer[size_];
        table_tinyv = new packed_tiny_vector[size_];
        table_tinyv_sz = new uint8_t[size_];

        if (o.table_outliers_psl != nullptr) {

            table_outliers_psl = new uint64_t[(size_ + 63) / 64];

            std::copy(o.table_outliers_psl, o.table_outliers_psl + size_, table_outliers_psl);
        }

        std::copy(o.table_keys, o.table_keys + size_, table_keys);

        for (size_t i = 0; i < size_; ++i){

            table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;
            table_tinyv[i].copy(table_tinyv_sz[i], o.table_tinyv[i], o.table_tinyv_sz[i]);
        }
    }
}

MinimizerIndex::MinimizerIndex(MinimizerIndex&& o) :    table_keys(o.table_keys), table_tinyv(o.table_tinyv), table_tinyv_sz(o.table_tinyv_sz),
                                                        size_(o.size_), pop(o.pop), table_outliers_psl(o.table_outliers_psl),
                                                        sum_psl(o.sum_psl), max_psl(o.max_psl), std_psl(o.std_psl),
                                                        M_u64(o.M_u64), max_ratio_occupancy(o.max_ratio_occupancy) {

    o.table_keys = nullptr;
    o.table_tinyv = nullptr;
    o.table_tinyv_sz = nullptr;
    o.table_outliers_psl = nullptr;

    o.clear();
}

MinimizerIndex& MinimizerIndex::operator=(const MinimizerIndex& o) {

    if (this != &o) {

        clear();

        size_ = o.size_;
        pop = o.pop;
        sum_psl = o.sum_psl;
        std_psl = o.std_psl;
        max_psl = o.max_psl;
        max_ratio_occupancy = o.max_ratio_occupancy;
        M_u64 = o.M_u64;

        if (size_ != 0) {

            table_keys = new Minimizer[size_];
            table_tinyv = new packed_tiny_vector[size_];
            table_tinyv_sz = new uint8_t[size_];

            if (o.table_outliers_psl != nullptr) {

                table_outliers_psl = new uint64_t[(size_ + 63) / 64];

                std::copy(o.table_outliers_psl, o.table_outliers_psl + size_, table_outliers_psl);
            }

            std::copy(o.table_keys, o.table_keys + size_, table_keys);

            for (size_t i = 0; i < size_; ++i){

                table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;
                table_tinyv[i].copy(table_tinyv_sz[i], o.table_tinyv[i], o.table_tinyv_sz[i]);
            }
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
        std_psl = o.std_psl;
        max_psl = o.max_psl;
        max_ratio_occupancy = o.max_ratio_occupancy;
        M_u64 = o.M_u64;

        table_keys = o.table_keys;
        table_tinyv = o.table_tinyv;
        table_tinyv_sz = o.table_tinyv_sz;
        table_outliers_psl = o.table_outliers_psl;

        o.table_keys = nullptr;
        o.table_tinyv = nullptr;
        o.table_tinyv_sz = nullptr;
        o.table_outliers_psl = nullptr;

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

    if (table_outliers_psl != nullptr) {

        delete[] table_outliers_psl;
        table_outliers_psl = nullptr;
    }

    size_ = 0;
    pop  = 0;
    sum_psl = 0;
    std_psl = 0;
    max_psl = 0;
    M_u64 = 0;

    max_ratio_occupancy = BIFROST_MI_MAX_OCCUPANCY;
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

MinimizerIndex::iterator MinimizerIndex::find(const Minimizer& key) {

    if ((pop != 0) && (size_ != 0)) {

        const size_t end_table = size_-1;
        const size_t mean_psl = get_mean_psl();

        size_t h = fastmod::fastmod_u64(key.hash(), M_u64, size_);
        size_t l_max_psl = max_psl;
        size_t l_mean_psl = mean_psl;

        // If std psl was computed, is not 0 and position in table is not annotated as outlier
        if ((table_outliers_psl != nullptr) && ((table_outliers_psl[h >> 6] & (0x1ULL << (h & 0x3fULL))) == 0)) {

            l_max_psl = mean_psl + std_psl;
            l_mean_psl = std_psl;
        }

        if (mean_psl <= 2) {

            size_t psl = 0;

            while ((psl <= l_max_psl) && !table_keys[h].isEmpty() && (table_keys[h] != key)) {

                h = (h + 1) & (static_cast<size_t>(h == end_table) - 1);
                ++psl;
            }

            if (table_keys[h] == key) return iterator(this, h, psl);
        }
        else {

            const size_t h_mean = fastmod::fastmod_u64(h + mean_psl, M_u64, size_);

            size_t h_inc = h_mean;
            size_t h_dec = h_mean;
            size_t i = 0;

            bool has_empty_key = false;

            // Check all elements located at positions mean+i and mean-i. Stop if empty key encountered. Stop if minimum key is encountered.
            for (; !has_empty_key && (i <= l_mean_psl); ++i) {

                if (table_keys[h_dec] == key) return iterator(this, h_dec, mean_psl - i);
                if (table_keys[h_inc] == key) return iterator(this, h_inc, mean_psl + i);

                has_empty_key = table_keys[h_dec].isEmpty() || table_keys[h_inc].isEmpty();

                h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
                h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
            }

            if (!has_empty_key) {

                // Only check remaining acending positions if neither the empty key nor the minimum key were encountered
                for (size_t j = i; (mean_psl + j <= l_max_psl) && !table_keys[h_inc].isEmpty(); ++j) {

                    if (table_keys[h_inc] == key) return iterator(this, h_inc, mean_psl + j);

                    h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
                }
            }

            for (; i <= mean_psl; ++i) {

                if (table_keys[h_dec] == key) return iterator(this, h_dec, mean_psl - i);

                h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
            }
        }
    }

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const Minimizer& key) const {

    if ((pop != 0) && (size_ != 0)) {

        const size_t end_table = size_-1;
        const size_t mean_psl = get_mean_psl();

        size_t h = fastmod::fastmod_u64(key.hash(), M_u64, size_);
        size_t l_max_psl = max_psl;
        size_t l_mean_psl = mean_psl;

        // If std psl was computed, is not 0 and position in table is not annotated as outlier
        if ((table_outliers_psl != nullptr) && ((table_outliers_psl[h >> 6] & (0x1ULL << (h & 0x3fULL))) == 0)) {

            l_max_psl = mean_psl + std_psl;
            l_mean_psl = std_psl;
        }

        if (mean_psl <= 2) {

            size_t psl = 0;

            while ((psl <= l_max_psl) && !table_keys[h].isEmpty() && (table_keys[h] != key)) {

                h = (h + 1) & (static_cast<size_t>(h == end_table) - 1);
                ++psl;
            }

            if (table_keys[h] == key) return const_iterator(this, h, psl);
        }
        else {

            const size_t h_mean = fastmod::fastmod_u64(h + mean_psl, M_u64, size_);

            size_t h_inc = h_mean;
            size_t h_dec = h_mean;
            size_t i = 0;

            bool has_empty_key = false;

            // Check all elements located at positions mean+i and mean-i. Stop if empty key encountered. Stop if minimum key is encountered.
            for (; !has_empty_key && (i <= l_mean_psl); ++i) {

                if (table_keys[h_dec] == key) return const_iterator(this, h_dec, mean_psl - i);
                if (table_keys[h_inc] == key) return const_iterator(this, h_inc, mean_psl + i);

                has_empty_key = table_keys[h_dec].isEmpty() || table_keys[h_inc].isEmpty();

                h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
                h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
            }

            if (!has_empty_key) {

                // Only check remaining acending positions if neither the empty key nor the minimum key were encountered
                for (size_t j = i; (mean_psl + j <= l_max_psl) && !table_keys[h_inc].isEmpty(); ++j) {

                    if (table_keys[h_inc] == key) return const_iterator(this, h_inc, mean_psl + j);

                    h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
                }
            }

            for (; i <= mean_psl; ++i) {

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

    if ((size_ == 0) || (it == end())) return 0;

    const size_t end_table = size_-1;

    // Remove <key, value>
    {
        if (it.psl != 0xffffffffffffffffULL) sum_psl -= it.psl;
        else {

            const size_t h = fastmod::fastmod_u64(table_keys[it.h].hash(), M_u64, size_);

            sum_psl -= (size_ - h + it.h) & (static_cast<size_t>(it.h >= h) - 1);
            sum_psl -= (it.h - h) & (static_cast<size_t>(it.h < h) - 1);
        }

        table_keys[it.h].set_empty();
        table_tinyv[it.h].destruct(table_tinyv_sz[it.h]);
        table_tinyv_sz[it.h] = packed_tiny_vector::FLAG_EMPTY;

        --pop;
    }

    // Computed psl standard deviation is only valid if no insertion/deletion has happened since it was computed
    // Current deletion invalidate the psl standard deviation
    if (table_outliers_psl != nullptr) {

        std_psl = 0;
        
        delete[] table_outliers_psl;
        table_outliers_psl = nullptr;
    }

    // Robin-hood hash: push the tombstone further away if subsequent keys can be closer to where they are supposed to be
    {
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
    }

    return 1;
}

// Insert with Robin Hood hashing
pair<MinimizerIndex::iterator, bool> MinimizerIndex::insert(const Minimizer& key, const packed_tiny_vector& ptv, const uint8_t& flag) {

    if (size_ == 0) init_tables(BIFROST_MI_INIT_SZ);
    else if (pop >= static_cast<size_t>(size_ * max_ratio_occupancy)) {

        size_t resize = max(1.2 * size_, static_cast<double>(1 + size_));

        while (pop >= static_cast<size_t>(resize * max_ratio_occupancy)) resize = max(1.2 * resize, static_cast<double>(1 + resize));

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

        if (table_keys[h].isEmpty() || (has_rich_psl && (cascade_ins || (psl_ins_key > max_psl)))) {

            if (has_rich_psl) {

                packed_tiny_vector l_ptv_swap;
                uint8_t l_flag_swap = packed_tiny_vector::FLAG_EMPTY;

                h = h_rich_psl_ins;

                std::swap(table_keys[h], l_key);

                l_ptv_swap.move(l_flag_swap, move(table_tinyv[h]), move(table_tinyv_sz[h]));
                table_tinyv[h].move(table_tinyv_sz[h], move(l_ptv), move(l_flag));
                l_ptv.move(l_flag, move(l_ptv_swap), move(l_flag_swap));

                if (!cascade_ins) it_ret = {iterator(this, h, psl_rich_key), true};

                max_psl = max(max_psl, psl_rich_key);
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

                max_psl = max(max_psl, psl_ins_key);
                sum_psl += psl_ins_key;

                if (!cascade_ins) it_ret = {iterator(this, h, psl_ins_key), true};

                // Computed psl standard deviation is only valid if no insertion/deletion has happened since it was computed
                // Current insertion invalidate the psl standard deviation
                if (table_outliers_psl != nullptr) {

                    std_psl = 0;
                    
                    delete[] table_outliers_psl;
                    table_outliers_psl = nullptr;
                }

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

    max_psl = 0;

    if (pop != 0) {

        if (nb_threads <= 1){

            for (size_t i = 0; i != size_; ++i) {

                if (!table_keys[i].isEmpty()) {

                    const size_t h = fastmod::fastmod_u64(table_keys[i].hash(), M_u64, size_);
                    const size_t psl = ((size_ - h + i) & (static_cast<size_t>(i >= h) - 1)) + ((i - h) & (static_cast<size_t>(i < h) - 1));

                    max_psl = max(max_psl, psl);
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

                        size_t l_max_psl = 0;

                        for (size_t i = chunk_start; i < chunk_end; ++i) {

                            if (!table_keys[i].isEmpty()) {

                                const size_t h = fastmod::fastmod_u64(table_keys[i].hash(), M_u64, size_);
                                const size_t psl = ((size_ - h + i) & (static_cast<size_t>(i >= h) - 1)) + ((i - h) & (static_cast<size_t>(i < h) - 1));

                                l_max_psl = max(l_max_psl, psl);
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

void MinimizerIndex::recomputeMaxStdPSL(const size_t nb_threads) {

    std_psl = 0;
    max_psl = 0;

    if (table_outliers_psl != nullptr) {
        
        delete[] table_outliers_psl;
        table_outliers_psl = nullptr;
    }

    if (pop != 0) {

        const size_t mean_psl = get_mean_psl();

        if (nb_threads <= 1) {

            for (size_t i = 0; i != size_; ++i) {

                if (!table_keys[i].isEmpty()) {

                    const size_t h = fastmod::fastmod_u64(table_keys[i].hash(), M_u64, size_);
                    const size_t psl = ((size_ - h + i) & (static_cast<size_t>(i >= h) - 1)) + ((i - h) & (static_cast<size_t>(i < h) - 1));

                    max_psl = max(max_psl, psl);
                    std_psl += pow(abs(static_cast<int64_t>(psl) - static_cast<int64_t>(mean_psl)), 2);
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

                        size_t l_max_psl = 0;
                        size_t l_std_psl = 0;

                        for (size_t i = chunk_start; i < chunk_end; ++i) {

                            if (!table_keys[i].isEmpty()) {

                                const size_t h = fastmod::fastmod_u64(table_keys[i].hash(), M_u64, size_);
                                const size_t psl = ((size_ - h + i) & (static_cast<size_t>(i >= h) - 1)) + ((i - h) & (static_cast<size_t>(i < h) - 1));

                                l_max_psl = max(l_max_psl, psl);
                                l_std_psl += pow(abs(static_cast<int64_t>(psl) - static_cast<int64_t>(mean_psl)), 2);
                            }
                        }

                        {
                            unique_lock<mutex> lock(mtx_max_psl);

                            max_psl = max(max_psl, l_max_psl);
                            std_psl += l_std_psl;
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }

        // Finish computing psl standard deviation
        std_psl = static_cast<size_t>(ceil(sqrt(static_cast<double>(std_psl)/static_cast<double>(pop))));

        // Update outlier psl list
        if (std_psl != 0) {

            table_outliers_psl = new uint64_t[(size_ + 63) / 64]();

            if (nb_threads <= 1) {

                for (size_t i = 0; i != size_; ++i) {

                    if (!table_keys[i].isEmpty()) {

                        const size_t h = fastmod::fastmod_u64(table_keys[i].hash(), M_u64, size_);
                        const size_t psl = ((size_ - h + i) & (static_cast<size_t>(i >= h) - 1)) + ((i - h) & (static_cast<size_t>(i < h) - 1));

                        if (psl > mean_psl + std_psl) table_outliers_psl[h >> 6] |= 0x1ULL << (h & 0x3fULL);
                    }
                }
            }
            else {

                const size_t chunk_per_thread = (size_ + nb_threads - 1) / nb_threads;

                vector<thread> workers; // need to keep track of threads so we can join them

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&, t]{

                            const size_t chunk_start = t * chunk_per_thread;
                            const size_t chunk_end = min(((t+1) * chunk_per_thread), size_);

                            for (size_t i = chunk_start; i < chunk_end; ++i) {

                                if (!table_keys[i].isEmpty()) {

                                    const size_t h = fastmod::fastmod_u64(table_keys[i].hash(), M_u64, size_);
                                    const size_t psl = ((size_ - h + i) & (static_cast<size_t>(i >= h) - 1)) + ((i - h) & (static_cast<size_t>(i < h) - 1));

                                    if (psl > mean_psl + std_psl) table_outliers_psl[h >> 6] |= 0x1ULL << (h & 0x3fULL);
                                }
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }
        }
    }
}

MinimizerIndex::iterator MinimizerIndex::begin() {

    iterator it(this);

    it.get_to_first();

    return it;
}

MinimizerIndex::const_iterator MinimizerIndex::begin() const {

    const_iterator it(this);

    it.get_to_first();

    return it;
}

MinimizerIndex::iterator MinimizerIndex::end() {

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::end() const {

    return const_iterator(this);
}

void MinimizerIndex::reserve(const size_t sz) {

    if (sz <= size_) return;

    if (size_ == 0) init_tables(sz);
    else {

        const size_t old_size_ = size_;

        Minimizer empty_key;

        Minimizer* old_table_keys = table_keys;
        packed_tiny_vector* old_table_tinyv = table_tinyv;
        uint8_t* old_table_tinyv_sz = table_tinyv_sz;

        size_ = sz;
        pop = 0;
        sum_psl = 0;
        std_psl = 0;
        max_psl = 0;

        M_u64 = fastmod::computeM_u64(size_);

         if (table_outliers_psl != nullptr){

            delete[] table_outliers_psl;

            table_outliers_psl = nullptr;
         }

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
}

void MinimizerIndex::swap(const size_t i, const size_t j) {

    uint8_t ptv_sz = packed_tiny_vector::FLAG_EMPTY;

    packed_tiny_vector ptv;

    ptv.move(ptv_sz, move(table_tinyv[i]), move(table_tinyv_sz[i]));
    table_tinyv[i].move(table_tinyv_sz[i], move(table_tinyv[j]), move(table_tinyv_sz[j]));
    table_tinyv[j].move(table_tinyv_sz[j], move(ptv), move(ptv_sz));

    std::swap(table_keys[i], table_keys[j]);
}