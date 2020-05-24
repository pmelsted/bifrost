#ifndef BIFROST_SEARCH_DBG_TCC
#define BIFROST_SEARCH_DBG_TCC


template<typename U, typename G>
vector<pair<size_t, UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence(   const string& seq, const bool exact, const bool insertion,
                                                                            const bool deletion, const bool substitution,
                                                                            const bool or_exclusive_match) {

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    vector<pair<size_t, UnitigMap<U, G>>> v_um;

    string seqs;

    Roaring r;

    auto worker_func = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const bool subst_or_ind = (subst || ins);
        const bool inexact = (subst_or_ind || del);

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t end = 1ULL << ((static_cast<size_t>(!subst_or_ind) - 1) & 0x2ULL);
        const size_t seq_len = seq.length();

        auto processUnitigMap = [&](const UnitigMap<U, G>& um, const size_t pos_seq){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_seq + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if (l_pos_seq + k_ - 1 < seq_len) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_seq + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if (l_pos_seq + k_ - 1 < seq_len) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
        };

        for (size_t i = 0; i != end; ++i){

            if (subst_or_ind){

                for (size_t j = shift; j < seqs.length(); j += k_) seqs[j] = alpha[i];
            }

            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(seqs.c_str(), seqs.length(), k_, g_, RepHash(), true);
            minHashResultIterator<RepHash> it_min = *mhi, it_min_end;
            minHashResult mhr = *it_min;

            Minimizer minz = Minimizer(seqs.c_str() + mhr.pos).rep();

            pair<size_t, bool> minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

            for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                mhr = *it_min;
                minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
            }

            size_t pos_seq = 0;
            size_t l_pos_seq = 0;
            size_t shift_pos_seq = 0;

            shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

            l_pos_seq -= (ins_mask & shift_pos_seq);
            l_pos_seq += (del_mask & shift_pos_seq);

            if (!inexact || !or_exclusive_match || !r.contains(l_pos_seq)){

                if (minz_pres.second){ // If at least one minimizer was present, search the kmer

                    const UnitigMap<U, G> um = findUnitig(seqs.c_str(), pos_seq, seqs.length(), mhi);

                    if (!um.isEmpty){

                        processUnitigMap(um, pos_seq);

                        pos_seq += um.len - 1;
                        mhi += pos_seq - mhi.getKmerPosition();
                    }
                }
            }

            ++pos_seq;
            ++mhi;

            while (pos_seq < seqs.length() - k_ + 1){

                shift_pos_seq = (pos_seq / k_) + (pos_seq % k_ > shift);
                l_pos_seq = pos_seq - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                if (!inexact || !or_exclusive_match || !r.contains(l_pos_seq)){

                    it_min = *mhi;
                    mhr = *it_min;

                    // If minimizers of new kmer are different from minimizers of previous kmer
                    // or if minimizers are the same but they were present, search them again
                    if ((mhr.pos != minz_pres.first) || minz_pres.second){

                        if (mhr.pos != minz_pres.first){

                            minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                            minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                            for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                                mhr = *it_min;
                                minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                                minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                            }
                        }

                        if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                            const UnitigMap<U, G> um = findUnitig(seqs.c_str(), pos_seq, seqs.length(), mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_seq);

                                pos_seq += um.len - 1;
                                mhi += pos_seq - mhi.getKmerPosition();
                            }
                        }
                    }
                }

                ++pos_seq;
                ++mhi;
            }
        }
    };

    if (exact){

        for (size_t i = 0; i < seq.length() - k_ + 1; ++i) {

            const UnitigMap<U, G> um = findUnitig(seq.c_str(), i, seq.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({i + j - um.dist, um.getKmerMapping(j)});
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({i + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                }

                i += um.len - 1;
            }
        }

        if (or_exclusive_match){

            for (const auto& p : v_um) r.add(p.first);
        }
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            seqs = seq;
            worker_func(true, false, false, i);
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << seq[j];

            for (size_t j = i, cpt = 0; j < seq.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << seq[j];
            }

            seqs = ss.str();

            worker_func(false, true, false, i);
        }
    }

    if (deletion && (seq.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << seq[j];

            for (size_t j = i, cpt = 0; j < seq.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << seq[j];
            }

            seqs = ss.str();

            worker_func(false, false, true, i);
        }
    }

    return v_um;
}

template<typename U, typename G>
vector<pair<size_t, UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence(   const string& seq, const bool exact, const bool insertion,
                                                                            const bool deletion, const bool substitution,
                                                                            const double ratio_kmers, const bool or_exclusive_match) {

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (ratio_kmers < 0.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is less than 0.0" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (ratio_kmers > 1.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is greater than 1.0" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    const size_t nb_km_min = static_cast<double>(seq.length() - k_ + 1) * ratio_kmers;

    vector<pair<size_t, UnitigMap<U, G>>> v_um;

    string seqs;

    Roaring r, r_pos;

    auto worker_func = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const bool subst_or_ind = (subst || ins);
        const bool inexact = (subst_or_ind || del);

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t end = 1ULL << ((static_cast<size_t>(!subst_or_ind) - 1) & 0x2ULL);
        const size_t seq_len = seq.length();

        auto processUnitigMap = [&](const UnitigMap<U, G>& um, const size_t pos_seq){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_seq + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if (l_pos_seq + k_ - 1 < seq_len){

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        r_pos.add(l_pos_seq);
                    }
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_seq + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if (l_pos_seq + k_ - 1 < seq_len){

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        r_pos.add(l_pos_seq);
                    }
                }
            }
        };

        for (size_t i = 0; i != end; ++i){

            if (subst_or_ind){

                for (size_t j = shift; j < seqs.length(); j += k_) seqs[j] = alpha[i];
            }

            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(seqs.c_str(), seqs.length(), k_, g_, RepHash(), true);
            minHashResultIterator<RepHash> it_min = *mhi, it_min_end;
            minHashResult mhr = *it_min;

            Minimizer minz = Minimizer(seqs.c_str() + mhr.pos).rep();

            pair<size_t, bool> minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

            for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                mhr = *it_min;
                minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
            }

            size_t pos_seq = 0;
            size_t l_pos_seq = 0;
            size_t shift_pos_seq = 0;

            shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

            l_pos_seq -= (ins_mask & shift_pos_seq);
            l_pos_seq += (del_mask & shift_pos_seq);

            if (!inexact || !or_exclusive_match || !r.contains(l_pos_seq)){

                if (minz_pres.second){ // If at least one minimizer was present, search the kmer

                    const UnitigMap<U, G> um = findUnitig(seqs.c_str(), pos_seq, seqs.length(), mhi);

                    if (!um.isEmpty){

                        processUnitigMap(um, pos_seq);

                        if (r_pos.cardinality() >= nb_km_min) return;

                        pos_seq += um.len - 1;
                        mhi += pos_seq - mhi.getKmerPosition();
                    }
                }
            }

            ++pos_seq;
            ++mhi;

            while (pos_seq < seqs.length() - k_ + 1){

                shift_pos_seq = (pos_seq / k_) + (pos_seq % k_ > shift);
                l_pos_seq = pos_seq - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                if (!inexact || !or_exclusive_match || !r.contains(l_pos_seq)){

                    it_min = *mhi;
                    mhr = *it_min;

                    // If minimizers of new kmer are different from minimizers of previous kmer
                    // or if minimizers are the same but they were present, search them again
                    if ((mhr.pos != minz_pres.first) || minz_pres.second){

                        if (mhr.pos != minz_pres.first){

                            minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                            minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                            for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                                mhr = *it_min;
                                minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                                minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                            }
                        }

                        if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                            const UnitigMap<U, G> um = findUnitig(seqs.c_str(), pos_seq, seqs.length(), mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_seq);

                                if (r_pos.cardinality() >= nb_km_min) return;

                                pos_seq += um.len - 1;
                                mhi += pos_seq - mhi.getKmerPosition();
                            }
                        }
                    }
                }

                ++pos_seq;
                ++mhi;
            }
        }
    };

    if (exact){

        for (size_t i = 0; i < seq.length() - k_ + 1; ++i) {

            const UnitigMap<U, G> um = findUnitig(seq.c_str(), i, seq.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j){

                        const size_t l_pos = i + j - um.dist;

                        v_um.push_back({l_pos, um.getKmerMapping(j)});
                        r_pos.add(l_pos);
                    }
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j){

                        const size_t l_pos = i + um.dist + um.len - j - 1;

                        v_um.push_back({l_pos, um.getKmerMapping(j)});
                        r_pos.add(l_pos);
                    }
                }

                if (r_pos.cardinality() >= nb_km_min) return v_um;

                i += um.len - 1;
            }
        }

        if (or_exclusive_match) r = r_pos;
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            seqs = seq;
            worker_func(true, false, false, i);

            if (r_pos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << seq[j];

            for (size_t j = i, cpt = 0; j < seq.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << seq[j];
            }

            seqs = ss.str();

            worker_func(false, true, false, i);

            if (r_pos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (deletion && (seq.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << seq[j];

            for (size_t j = i, cpt = 0; j < seq.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << seq[j];
            }

            seqs = ss.str();

            worker_func(false, false, true, i);

            if (r_pos.cardinality() >= nb_km_min) return v_um;
        }
    }

    return v_um;
}

template<typename U, typename G>
vector<pair<size_t, const_UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence(     const string& seq, const bool exact, const bool insertion,
                                                                                    const bool deletion, const bool substitution,
                                                                                    const bool or_exclusive_match) const {

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    vector<pair<size_t, const_UnitigMap<U, G>>> v_um;

    string seqs;

    Roaring r;

    auto worker_func = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const bool subst_or_ind = (subst || ins);
        const bool inexact = (subst_or_ind || del);

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t end = 1ULL << ((static_cast<size_t>(!subst_or_ind) - 1) & 0x2ULL);
        const size_t seq_len = seq.length();

        auto processUnitigMap = [&](const const_UnitigMap<U, G>& um, const size_t pos_seq){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_seq + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if (l_pos_seq + k_ - 1 < seq_len) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_seq + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if (l_pos_seq + k_ - 1 < seq_len) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
        };

        for (size_t i = 0; i != end; ++i){

            if (subst_or_ind){

                for (size_t j = shift; j < seqs.length(); j += k_) seqs[j] = alpha[i];
            }

            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(seqs.c_str(), seqs.length(), k_, g_, RepHash(), true);
            minHashResultIterator<RepHash> it_min = *mhi, it_min_end;
            minHashResult mhr = *it_min;

            Minimizer minz = Minimizer(seqs.c_str() + mhr.pos).rep();

            pair<size_t, bool> minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

            for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                mhr = *it_min;
                minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
            }

            size_t pos_seq = 0;
            size_t l_pos_seq = 0;
            size_t shift_pos_seq = 0;

            shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

            l_pos_seq -= (ins_mask & shift_pos_seq);
            l_pos_seq += (del_mask & shift_pos_seq);

            if (!inexact || !or_exclusive_match || !r.contains(l_pos_seq)){

                if (minz_pres.second){ // If at least one minimizer was present, search the kmer

                    const const_UnitigMap<U, G> um = findUnitig(seqs.c_str(), pos_seq, seqs.length(), mhi);

                    if (!um.isEmpty){

                        processUnitigMap(um, pos_seq);

                        pos_seq += um.len - 1;
                        mhi += pos_seq - mhi.getKmerPosition();
                    }
                }
            }

            ++pos_seq;
            ++mhi;

            while (pos_seq < seqs.length() - k_ + 1){

                shift_pos_seq = (pos_seq / k_) + (pos_seq % k_ > shift);
                l_pos_seq = pos_seq - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                if (!inexact || !or_exclusive_match || !r.contains(l_pos_seq)){

                    it_min = *mhi;
                    mhr = *it_min;

                    // If minimizers of new kmer are different from minimizers of previous kmer
                    // or if minimizers are the same but they were present, search them again
                    if ((mhr.pos != minz_pres.first) || minz_pres.second){

                        if (mhr.pos != minz_pres.first){

                            minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                            minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                            for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                                mhr = *it_min;
                                minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                                minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                            }
                        }

                        if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                            const const_UnitigMap<U, G> um = findUnitig(seqs.c_str(), pos_seq, seqs.length(), mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_seq);

                                pos_seq += um.len - 1;
                                mhi += pos_seq - mhi.getKmerPosition();
                            }
                        }
                    }
                }

                ++pos_seq;
                ++mhi;
            }
        }
    };

    if (exact){

        for (size_t i = 0; i < seq.length() - k_ + 1; ++i) {

            const const_UnitigMap<U, G> um = findUnitig(seq.c_str(), i, seq.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({i + j - um.dist, um.getKmerMapping(j)});
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({i + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                }

                i += um.len - 1;
            }
        }

        if (or_exclusive_match){

            for (const auto& p : v_um) r.add(p.first);
        }
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            seqs = seq;
            worker_func(true, false, false, i);
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << seq[j];

            for (size_t j = i, cpt = 0; j < seq.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << seq[j];
            }

            seqs = ss.str();

            worker_func(false, true, false, i);
        }
    }

    if (deletion && (seq.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << seq[j];

            for (size_t j = i, cpt = 0; j < seq.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << seq[j];
            }

            seqs = ss.str();

            worker_func(false, false, true, i);
        }
    }

    return v_um;
}

template<typename U, typename G>
vector<pair<size_t, const_UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence(     const string& seq, const bool exact, const bool insertion,
                                                                                    const bool deletion, const bool substitution,
                                                                                    const double ratio_kmers, const bool or_exclusive_match) const {

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (ratio_kmers < 0.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is less than 0.0" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (ratio_kmers > 1.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is greater than 1.0" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    const size_t nb_km_min = static_cast<double>(seq.length() - k_ + 1) * ratio_kmers;

    vector<pair<size_t, const_UnitigMap<U, G>>> v_um;

    string seqs;

    Roaring r, r_pos;

    auto worker_func = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const bool subst_or_ind = (subst || ins);
        const bool inexact = (subst_or_ind || del);

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t end = 1ULL << ((static_cast<size_t>(!subst_or_ind) - 1) & 0x2ULL);
        const size_t seq_len = seq.length();

        auto processUnitigMap = [&](const const_UnitigMap<U, G>& um, const size_t pos_seq){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_seq + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if (l_pos_seq + k_ - 1 < seq_len){

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        r_pos.add(l_pos_seq);
                    }
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_seq + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if (l_pos_seq + k_ - 1 < seq_len){

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        r_pos.add(l_pos_seq);
                    }
                }
            }
        };

        for (size_t i = 0; i != end; ++i){

            if (subst_or_ind){

                for (size_t j = shift; j < seqs.length(); j += k_) seqs[j] = alpha[i];
            }

            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(seqs.c_str(), seqs.length(), k_, g_, RepHash(), true);
            minHashResultIterator<RepHash> it_min = *mhi, it_min_end;
            minHashResult mhr = *it_min;

            Minimizer minz = Minimizer(seqs.c_str() + mhr.pos).rep();

            pair<size_t, bool> minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

            for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                mhr = *it_min;
                minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
            }

            size_t pos_seq = 0;
            size_t l_pos_seq = 0;
            size_t shift_pos_seq = 0;

            shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

            l_pos_seq -= (ins_mask & shift_pos_seq);
            l_pos_seq += (del_mask & shift_pos_seq);

            if (!inexact || !or_exclusive_match || !r.contains(l_pos_seq)){

                if (minz_pres.second){ // If at least one minimizer was present, search the kmer

                    const const_UnitigMap<U, G> um = findUnitig(seqs.c_str(), pos_seq, seqs.length(), mhi);

                    if (!um.isEmpty){

                        processUnitigMap(um, pos_seq);

                        if (r_pos.cardinality() >= nb_km_min) return;

                        pos_seq += um.len - 1;
                        mhi += pos_seq - mhi.getKmerPosition();
                    }
                }
            }

            ++pos_seq;
            ++mhi;

            while (pos_seq < seqs.length() - k_ + 1){

                shift_pos_seq = (pos_seq / k_) + (pos_seq % k_ > shift);
                l_pos_seq = pos_seq - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                if (!inexact || !or_exclusive_match || !r.contains(l_pos_seq)){

                    it_min = *mhi;
                    mhr = *it_min;

                    // If minimizers of new kmer are different from minimizers of previous kmer
                    // or if minimizers are the same but they were present, search them again
                    if ((mhr.pos != minz_pres.first) || minz_pres.second){

                        if (mhr.pos != minz_pres.first){

                            minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                            minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                            for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                                mhr = *it_min;
                                minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                                minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                            }
                        }

                        if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                            const const_UnitigMap<U, G> um = findUnitig(seqs.c_str(), pos_seq, seqs.length(), mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_seq);

                                if (r_pos.cardinality() >= nb_km_min) return;

                                pos_seq += um.len - 1;
                                mhi += pos_seq - mhi.getKmerPosition();
                            }
                        }
                    }
                }

                ++pos_seq;
                ++mhi;
            }
        }
    };

    if (exact){

        for (size_t i = 0; i < seq.length() - k_ + 1; ++i) {

            const const_UnitigMap<U, G> um = findUnitig(seq.c_str(), i, seq.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j){

                        const size_t l_pos = i + j - um.dist;

                        v_um.push_back({l_pos, um.getKmerMapping(j)});
                        r_pos.add(l_pos);
                    }
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j){

                        const size_t l_pos = i + um.dist + um.len - j - 1;

                        v_um.push_back({l_pos, um.getKmerMapping(j)});
                        r_pos.add(l_pos);
                    }
                }

                if (r_pos.cardinality() >= nb_km_min) return v_um;

                i += um.len - 1;
            }
        }

        if (or_exclusive_match) r = r_pos;
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            seqs = seq;
            worker_func(true, false, false, i);

            if (r_pos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << seq[j];

            for (size_t j = i, cpt = 0; j < seq.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << seq[j];
            }

            seqs = ss.str();

            worker_func(false, true, false, i);

            if (r_pos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (deletion && (seq.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << seq[j];

            for (size_t j = i, cpt = 0; j < seq.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << seq[j];
            }

            seqs = ss.str();

            worker_func(false, false, true, i);

            if (r_pos.cardinality() >= nb_km_min) return v_um;
        }
    }

    return v_um;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::search(const vector<string>& query_filenames, const string& out_filename_prefix,
                                const double ratio_kmers, const bool inexact_search, const size_t nb_threads,
                                const size_t verbose) const {

     if (invalid){

        cerr << "CompactedDBG::search(): Graph is invalid and cannot be searched" << endl;
        return false;
    }

    if (nb_threads > std::thread::hardware_concurrency()){

        cerr << "CompactedDBG::search(): Number of threads cannot be greater than or equal to " << std::thread::hardware_concurrency() << "." << endl;
        return false;
    }

    if (nb_threads <= 0){

        cerr << "CompactedDBG::search(): Number of threads cannot be less than or equal to 0." << endl;
        return false;
    }

    const string out_tmp = out_filename_prefix + ".tsv";

    FILE* fp_tmp = fopen(out_tmp.c_str(), "w");

    if (fp_tmp == NULL) {

        cerr << "CompactedDBG::search(): Could not open file " << out_tmp << " for writing." << endl;
        return false;
    }
    else {

        fclose(fp_tmp);

        if (std::remove(out_tmp.c_str()) != 0) cerr << "CompactedDBG::search(): Could not remove temporary file " << out_tmp << endl;
    }

    if (verbose) cout << "CompactedDBG::search(): Querying graph." << endl;

    const CompactedDBG<U, G>& dbg = *this;

    string s;

    size_t file_id = 0;

    const size_t max_len_seq = 1024;
    const size_t thread_seq_buf_sz = 64 * max_len_seq;

    FileParser fp(query_filenames);

    ofstream outfile;
    ostream out(0);

    outfile.open(out_tmp.c_str());
    out.rdbuf(outfile.rdbuf());
    //out.sync_with_stdio(false);

    const char query_pres[3] = {'\t', '1', '\n'};
    const char query_abs[3] = {'\t', '0', '\n'};

    const size_t l_query_res = 3;

    // Write header to TSV file
    out << "query_name\tpresence_query\n";

    if (nb_threads == 1){

        char* buffer_res = new char[thread_seq_buf_sz];

        size_t pos_buffer_out = 0;
        size_t nb_queries_found = 0;

        while (fp.read(s, file_id)){

            bool is_found = false;

            const size_t nb_km_min = static_cast<double>(s.length() - k_ + 1) * ratio_kmers;
            const char* query_name = fp.getNameString();
            const size_t l_query_name = strlen(query_name);

            for (auto& c : s) c &= 0xDF;

            const vector<pair<size_t, const_UnitigMap<U, G>>> v = dbg.searchSequence(   s, true, inexact_search, inexact_search,
                                                                                        inexact_search, ratio_kmers, true);

            if (inexact_search){

                Roaring r;

                for (const auto& p : v) r.add(p.first);

                is_found = (r.cardinality() >= nb_km_min);
            }
            else is_found = (v.size() >= nb_km_min);

            if (pos_buffer_out + l_query_name + l_query_res >= thread_seq_buf_sz){ // If next result cannot fit in the buffer

                out.write(buffer_res, pos_buffer_out); // Write result buffer
                pos_buffer_out = 0; // Reset position to 0;
            }

            // Copy new result to buffer
            std::memcpy(buffer_res + pos_buffer_out, query_name, l_query_name * sizeof(char));

            if (is_found){

                std::memcpy(buffer_res + pos_buffer_out + l_query_name, query_pres, l_query_res * sizeof(char));

                ++nb_queries_found;
            }
            else std::memcpy(buffer_res + pos_buffer_out + l_query_name, query_abs, l_query_res * sizeof(char));

            pos_buffer_out += l_query_name + l_query_res;
        }

        // Flush unresult written to final output
        if (pos_buffer_out > 0) out.write(buffer_res, pos_buffer_out);

        delete[] buffer_res;

        if (verbose) cout << "CompactedDBG::search(): Found " << nb_queries_found << " queries. " << endl;
    }
    else {

        {
            bool stop = false;

            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mutex_files_in, mutex_file_out;

            std::atomic<size_t> nb_queries_found;

            nb_queries_found = 0;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        char* buffer_res = new char[nb_threads];

                        vector<string> buffers_seq;
                        vector<string> buffers_name;

                        while (true) {

                            {
                                if (stop) {

                                    delete[] buffer_res;

                                    return;
                                }

                                size_t buffer_sz = 0;

                                unique_lock<mutex> lock(mutex_files_in);

                                stop = !fp.read(s, file_id);

                                while (!stop){

                                    buffer_sz += s.length();

                                    buffers_seq.push_back(std::move(s));
                                    buffers_name.push_back(string(fp.getNameString()));

                                    if (buffer_sz >= thread_seq_buf_sz) break;
                                    else stop = !fp.read(s, file_id);
                                }
                            }

                            size_t pos_buffer_out = 0;

                            const size_t buffers_seq_sz = buffers_seq.size();

                            for (size_t i = 0; i < buffers_seq_sz; ++i){

                                bool is_found = false;

                                const size_t nb_km_min = static_cast<double>(buffers_seq[i].length() - k_ + 1) * ratio_kmers;
                                const size_t l_name = buffers_name[i].length();

                                for (auto& c : buffers_seq[i]) c &= 0xDF;

                                const vector<pair<size_t, const_UnitigMap<U, G>>> v = dbg.searchSequence(   buffers_seq[i], true, inexact_search, inexact_search,
                                                                                                            inexact_search, ratio_kmers, true);

                                if (inexact_search){

                                    Roaring r;

                                    for (const auto& p : v) r.add(p.first);

                                    is_found = (r.cardinality() >= nb_km_min);
                                }
                                else is_found = (v.size() >= nb_km_min);

                                if (pos_buffer_out + l_name + l_query_res >= thread_seq_buf_sz){ // If next result cannot fit in the buffer

                                    unique_lock<mutex> lock(mutex_file_out); // Get the output lock

                                    out.write(buffer_res, pos_buffer_out); // Write result buffer

                                    pos_buffer_out = 0; // Reset position to 0;
                                }

                                // Copy new result to buffer
                                std::memcpy(buffer_res + pos_buffer_out, buffers_name[i].c_str(), l_name * sizeof(char));

                                if (is_found){

                                    std::memcpy(buffer_res + pos_buffer_out + l_name, query_pres, l_query_res * sizeof(char));

                                    ++nb_queries_found;
                                }
                                else std::memcpy(buffer_res + pos_buffer_out + l_name, query_abs, l_query_res * sizeof(char));

                                pos_buffer_out += l_name + l_query_res;
                            }

                            if (pos_buffer_out > 0){ // Flush unresult written to final output

                                unique_lock<mutex> lock(mutex_file_out);

                                out.write(buffer_res, pos_buffer_out);
                            }

                            // Clear buffers for next round
                            buffers_seq.clear();
                            buffers_name.clear();
                        }

                        delete[] buffer_res;
                    }
                );
            }

            for (auto& t : workers) t.join();

            if (verbose) cout << "CompactedDBG::search(): Found " << nb_queries_found << " queries. " << endl;
        }
    }

    outfile.close();
    fp.close();

    return true;
}

#endif
