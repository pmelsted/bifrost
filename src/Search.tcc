#ifndef BIFROST_SEARCH_CDBG_TCC
#define BIFROST_SEARCH_CDBG_TCC

template<typename U, typename G>
vector<pair<size_t, UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence(   const string& s, const bool exact, const bool insertion,
                                                                            const bool deletion, const bool substitution,
                                                                            const bool or_exclusive_match) {

    struct hash_pair {

        size_t operator()(const pair<size_t, Kmer>& p) const {

            return wyhash(&(p.first), sizeof(size_t), 0, _wyp) ^ p.second.hash();
        }
    };

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (s.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    Roaring rpos;

    vector<pair<size_t, UnitigMap<U, G>>> v_um;

    string s_inexact;

    unordered_set<pair<size_t, Kmer>, hash_pair> us_pos_km;

    auto comp_pair = [](const pair<size_t, UnitigMap<U, G>>& p1, const pair<size_t, UnitigMap<U, G>>& p2) {

        return (p1.first < p2.first);
    };

    auto worker_inexact = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t s_len = s.length();
        const char* s_str = s.c_str();

        const size_t s_inexact_len = s_inexact.length();
        const char* s_inexact_str = s_inexact.c_str();

        const size_t k_1 = k_-1;

        auto processUnitigMap = [&](const UnitigMap<U, G>& um, const size_t pos_s){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
        };

        for (size_t i = 0; i != ((subst || ins) ? 4 : 1); ++i){

            if (ins) {

                for (size_t j = shift; j < s_inexact_len; j += k_) s_inexact[j] = alpha[i];
            }
            else if (subst) {

                for (size_t j = shift; j < s_inexact_len; j += k_) {

                    if (!isDNA(s[j]) || (alpha[i] == s[j])) s_inexact[j] = 'N';
                    else s_inexact[j] = alpha[i];
                }
            } 

            KmerIterator ki_s(s_inexact_str), ki_e;
            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(s_inexact_str, s_inexact_len, k_, g_, RepHash(), true);

            minHashResultIterator<RepHash> it_min, it_min_end;
            minHashResult mhr;

            Minimizer minz;

            pair<size_t, bool> minz_pres = {0xffffffffffffffffULL, true};

            while (ki_s != ki_e) {

                const size_t pos_s = ki_s->second;

                mhi += (pos_s - mhi.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                it_min = *mhi;
                mhr = *it_min;

                // If minimizers of new kmer are different from minimizers of previous kmer
                // or if minimizers are the same but they were present, search them again
                if (minz_pres.second || (mhr.pos != minz_pres.first)){

                    if (mhr.pos != minz_pres.first){

                        minz = Minimizer(s_inexact_str + mhr.pos).rep();
                        minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                        for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                            mhr = *it_min;
                            minz = Minimizer(s_inexact_str + mhr.pos).rep();
                            minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                        }
                    }

                    if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                        const size_t shift_pos_seq = (pos_s / k_) + (pos_s % k_ > shift);
                        const size_t l_pos_s = pos_s - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                        if ((l_pos_s + k_1 < s_len) && isDNA(s_str[l_pos_s]) && isDNA(s_str[l_pos_s + k_1]) && (!or_exclusive_match || (!rpos.contains(l_pos_s) && (us_pos_km.find({l_pos_s, ki_s->first}) == us_pos_km.end())))) {

                            const UnitigMap<U, G> um = findUnitig(s_inexact_str, pos_s, s_inexact_len, mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_s);

                                ki_s += um.len - 1;
                            }
                        }
                    }
                }

                ++ki_s;
            }
        }
    };

    if (exact){

        /*for (KmerIterator ki_s(s.c_str()), ki_e; ki_s != ki_e; ++ki_s) {

            const size_t pos_s = ki_s->second;
            const UnitigMap<U, G> um = findUnitig(s.c_str(), pos_s, s.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + j - um.dist, um.getKmerMapping(j)});
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                }

                ki_s += um.len - 1;
            }
        }*/

        {
            const size_t s_len = s.length();
            const char* s_str = s.c_str();

            KmerIterator ki_s(s_str), ki_e;
            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(s_str, s_len, k_, g_, RepHash(), true);

            minHashResultIterator<RepHash> it_min, it_min_end;
            minHashResult mhr;

            Minimizer minz;

            pair<size_t, bool> minz_pres = {0xffffffffffffffffULL, true};

            while (ki_s != ki_e) {

                const size_t pos_s = ki_s->second;

                mhi += (pos_s - mhi.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                it_min = *mhi;
                mhr = *it_min;

                // If minimizers of new kmer are different from minimizers of previous kmer
                // or if minimizers are the same but they were present, search them again
                if (minz_pres.second || (mhr.pos != minz_pres.first)){

                    if (mhr.pos != minz_pres.first){

                        minz = Minimizer(s_str + mhr.pos).rep();
                        minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                        for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                            mhr = *it_min;
                            minz = Minimizer(s_str + mhr.pos).rep();
                            minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                        }
                    }

                    if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                        const UnitigMap<U, G> um = findUnitig(s_str, pos_s, s_len, mhi);

                        if (!um.isEmpty) { // Read maps to a Unitig

                            if (um.strand){

                                for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + j - um.dist, um.getKmerMapping(j)});
                            }
                            else {

                                for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                            }

                            ki_s += um.len - 1;
                        }
                    }
                }

                ++ki_s;
            }
        }

        if (or_exclusive_match && (insertion || deletion || substitution)){

            for (const auto& pum : v_um) {

                us_pos_km.insert({pum.first, pum.second.getMappedKmer(pum.second.dist)});
                rpos.add(pum.first);
            }
        }
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            s_inexact = s;

            worker_inexact(true, false, false, i);
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << s[j];
            }

            s_inexact = ss.str();

            worker_inexact(false, true, false, i);
        }
    }

    if (deletion && (s.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << s[j];
            }

            s_inexact = ss.str();

            worker_inexact(false, false, true, i);
        }
    }

    return v_um;
}

template<typename U, typename G>
vector<pair<size_t, UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence(   const string& s, const bool exact, const bool insertion,
                                                                            const bool deletion, const bool substitution,
                                                                            const double ratio_kmers, const bool or_exclusive_match) {

    struct hash_pair {

        size_t operator()(const pair<size_t, Kmer>& p) const {

            return wyhash(&(p.first), sizeof(size_t), 0, _wyp) ^ p.second.hash();
        }
    };

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (ratio_kmers <= 0.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is less than or equal to 0.0" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (ratio_kmers > 1.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is greater than 1.0" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (s.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    const size_t nb_km_min = max(static_cast<size_t>(1), static_cast<size_t>(round(static_cast<double>(s.length() - k_ + 1) * ratio_kmers)));

    Roaring rpos;

    vector<pair<size_t, UnitigMap<U, G>>> v_um;

    string s_inexact;

    unordered_set<pair<size_t, Kmer>, hash_pair> us_pos_km;

    auto comp_pair = [](const pair<size_t, UnitigMap<U, G>>& p1, const pair<size_t, UnitigMap<U, G>>& p2) {

        return (p1.first < p2.first);
    };

    auto worker_func = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t s_len = s.length();
        const char* s_str = s.c_str();

        const size_t s_inexact_len = s_inexact.length();
        const char* s_inexact_str = s_inexact.c_str();

        const size_t k_1 = k_-1;

        auto processUnitigMap = [&](const UnitigMap<U, G>& um, const size_t pos_s){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) {

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        rpos.add(l_pos_seq);
                    }
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) {

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        rpos.add(l_pos_seq);
                    }
                }
            }
        };

        for (size_t i = 0; i != ((subst || ins) ? 4 : 1); ++i){

            if (ins) {

                for (size_t j = shift; j < s_inexact_len; j += k_) s_inexact[j] = alpha[i];
            }
            else if (subst) {

                for (size_t j = shift; j < s_inexact_len; j += k_) {

                    if (!isDNA(s[j]) || (alpha[i] == s[j])) s_inexact[j] = 'N';
                    else s_inexact[j] = alpha[i];
                }
            } 

            KmerIterator ki_s(s_inexact_str), ki_e;
            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(s_inexact_str, s_inexact_len, k_, g_, RepHash(), true);

            minHashResultIterator<RepHash> it_min, it_min_end;
            minHashResult mhr;

            Minimizer minz;

            pair<size_t, bool> minz_pres = {0xffffffffffffffffULL, true};

            while (ki_s != ki_e) {

                const size_t pos_s = ki_s->second;

                mhi += (pos_s - mhi.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                it_min = *mhi;
                mhr = *it_min;

                // If minimizers of new kmer are different from minimizers of previous kmer
                // or if minimizers are the same but they were present, search them again
                if (minz_pres.second || (mhr.pos != minz_pres.first)){

                    if (mhr.pos != minz_pres.first){

                        minz = Minimizer(s_inexact_str + mhr.pos).rep();
                        minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                        for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                            mhr = *it_min;
                            minz = Minimizer(s_inexact_str + mhr.pos).rep();
                            minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                        }
                    }

                    if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                        const size_t shift_pos_seq = (pos_s / k_) + (pos_s % k_ > shift);
                        const size_t l_pos_s = pos_s - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                        if ((l_pos_s + k_1 < s_len) && isDNA(s_str[l_pos_s]) && isDNA(s_str[l_pos_s + k_1]) && (!rpos.contains(l_pos_s) && (us_pos_km.find({l_pos_s, ki_s->first}) == us_pos_km.end()))) {

                            const UnitigMap<U, G> um = findUnitig(s_inexact_str, pos_s, s_inexact_len, mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_s);

                                if (rpos.cardinality() >= nb_km_min) return;

                                ki_s += um.len - 1;
                            }
                        }
                    }
                }

                ++ki_s;
            }
        }
    };

    if (exact){

        for (KmerIterator ki_s(s.c_str()), ki_e; ki_s != ki_e; ++ki_s) {

            const size_t pos_s = ki_s->second;
            const UnitigMap<U, G> um = findUnitig(s.c_str(), pos_s, s.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + j - um.dist, um.getKmerMapping(j)});
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                }

                if (v_um.size() >= nb_km_min) return v_um;

                ki_s += um.len - 1;
            }
        }

        for (const auto& pum : v_um) {

            us_pos_km.insert({pum.first, pum.second.getMappedKmer(pum.second.dist)});
            rpos.add(pum.first);
        }
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            s_inexact = s;

            worker_func(true, false, false, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, true, false, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (deletion && (s.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, false, true, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    return v_um;
}

template<typename U, typename G>
vector<pair<size_t, const_UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence( const string& s, const bool exact, const bool insertion,
                                                                                const bool deletion, const bool substitution,
                                                                                const bool or_exclusive_match) const {

    struct hash_pair {

        size_t operator()(const pair<size_t, Kmer>& p) const {

            return wyhash(&(p.first), sizeof(size_t), 0, _wyp) ^ p.second.hash();
        }
    };

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (s.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    Roaring rpos;

    vector<pair<size_t, const_UnitigMap<U, G>>> v_um;

    string s_inexact;

    unordered_set<pair<size_t, Kmer>, hash_pair> us_pos_km;

    auto comp_pair = [](const pair<size_t, const_UnitigMap<U, G>>& p1, const pair<size_t, const_UnitigMap<U, G>>& p2) {

        return (p1.first < p2.first);
    };

    auto worker_inexact = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t s_len = s.length();
        const char* s_str = s.c_str();

        const size_t s_inexact_len = s_inexact.length();
        const char* s_inexact_str = s_inexact.c_str();

        const size_t k_1 = k_-1;

        auto processUnitigMap = [&](const const_UnitigMap<U, G>& um, const size_t pos_s){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
        };

        for (size_t i = 0; i != ((subst || ins) ? 4 : 1); ++i){

            if (ins) {

                for (size_t j = shift; j < s_inexact_len; j += k_) s_inexact[j] = alpha[i];
            }
            else if (subst) {

                for (size_t j = shift; j < s_inexact_len; j += k_) {

                    if (!isDNA(s[j]) || (alpha[i] == s[j])) s_inexact[j] = 'N';
                    else s_inexact[j] = alpha[i];
                }
            } 

            KmerIterator ki_s(s_inexact_str), ki_e;
            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(s_inexact_str, s_inexact_len, k_, g_, RepHash(), true);

            minHashResultIterator<RepHash> it_min, it_min_end;
            minHashResult mhr;

            Minimizer minz;

            pair<size_t, bool> minz_pres = {0xffffffffffffffffULL, true};

            while (ki_s != ki_e) {

                const size_t pos_s = ki_s->second;

                mhi += (pos_s - mhi.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                it_min = *mhi;
                mhr = *it_min;

                // If minimizers of new kmer are different from minimizers of previous kmer
                // or if minimizers are the same but they were present, search them again
                if (minz_pres.second || (mhr.pos != minz_pres.first)){

                    if (mhr.pos != minz_pres.first){

                        minz = Minimizer(s_inexact_str + mhr.pos).rep();
                        minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                        for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                            mhr = *it_min;
                            minz = Minimizer(s_inexact_str + mhr.pos).rep();
                            minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                        }
                    }

                    if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                        const size_t shift_pos_seq = (pos_s / k_) + (pos_s % k_ > shift);
                        const size_t l_pos_s = pos_s - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                        if ((l_pos_s + k_1 < s_len) && isDNA(s_str[l_pos_s]) && isDNA(s_str[l_pos_s + k_1]) && (!or_exclusive_match || (!rpos.contains(l_pos_s) && (us_pos_km.find({l_pos_s, ki_s->first}) == us_pos_km.end())))) {

                            const const_UnitigMap<U, G> um = findUnitig(s_inexact_str, pos_s, s_inexact_len, mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_s);

                                ki_s += um.len - 1;
                            }
                        }
                    }
                }

                ++ki_s;
            }
        }
    };

    if (exact){

        /*for (KmerIterator ki_s(s.c_str()), ki_e; ki_s != ki_e; ++ki_s) {

            const size_t pos_s = ki_s->second;
            const const_UnitigMap<U, G> um = findUnitig(s.c_str(), pos_s, s.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + j - um.dist, um.getKmerMapping(j)});
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                }

                ki_s += um.len - 1;
            }
        }*/

        {
            const size_t s_len = s.length();
            const char* s_str = s.c_str();

            KmerIterator ki_s(s_str), ki_e;
            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(s_str, s_len, k_, g_, RepHash(), true);

            minHashResultIterator<RepHash> it_min, it_min_end;
            minHashResult mhr;

            Minimizer minz;

            pair<size_t, bool> minz_pres = {0xffffffffffffffffULL, true};

            while (ki_s != ki_e) {

                const size_t pos_s = ki_s->second;

                mhi += (pos_s - mhi.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                it_min = *mhi;
                mhr = *it_min;

                // If minimizers of new kmer are different from minimizers of previous kmer
                // or if minimizers are the same but they were present, search them again
                if (minz_pres.second || (mhr.pos != minz_pres.first)){

                    if (mhr.pos != minz_pres.first){

                        minz = Minimizer(s_str + mhr.pos).rep();
                        minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                        for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                            mhr = *it_min;
                            minz = Minimizer(s_str + mhr.pos).rep();
                            minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                        }
                    }

                    if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                        const const_UnitigMap<U, G> um = findUnitig(s_str, pos_s, s_len, mhi);

                        if (!um.isEmpty) { // Read maps to a Unitig

                            if (um.strand){

                                for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + j - um.dist, um.getKmerMapping(j)});
                            }
                            else {

                                for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                            }

                            ki_s += um.len - 1;
                        }
                    }
                }

                ++ki_s;
            }
        }

        if (or_exclusive_match && (insertion || deletion || substitution)){

            for (const auto& pum : v_um) {

                us_pos_km.insert({pum.first, pum.second.getMappedKmer(pum.second.dist)});
                rpos.add(pum.first);
            }
        }
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            s_inexact = s;

            worker_inexact(true, false, false, i);
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << s[j];
            }

            s_inexact = ss.str();

            worker_inexact(false, true, false, i);
        }
    }

    if (deletion && (s.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << s[j];
            }

            s_inexact = ss.str();

            worker_inexact(false, false, true, i);
        }
    }

    return v_um;
}

template<typename U, typename G>
vector<pair<size_t, const_UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence( const string& s, const bool exact, const bool insertion,
                                                                                const bool deletion, const bool substitution,
                                                                                const double ratio_kmers, const bool or_exclusive_match) const {

    struct hash_pair {

        size_t operator()(const pair<size_t, Kmer>& p) const {

            return wyhash(&(p.first), sizeof(size_t), 0, _wyp) ^ p.second.hash();
        }
    };

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (ratio_kmers <= 0.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is less than or equal to 0.0" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (ratio_kmers > 1.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is greater than 1.0" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (s.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    const size_t nb_km_min = max(static_cast<size_t>(1), static_cast<size_t>(round(static_cast<double>(s.length() - k_ + 1) * ratio_kmers)));

    Roaring rpos;

    vector<pair<size_t, const_UnitigMap<U, G>>> v_um;

    string s_inexact;

    unordered_set<pair<size_t, Kmer>, hash_pair> us_pos_km;

    auto comp_pair = [](const pair<size_t, const_UnitigMap<U, G>>& p1, const pair<size_t, const_UnitigMap<U, G>>& p2) {

        return (p1.first < p2.first);
    };

    auto worker_func = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t s_len = s.length();
        const char* s_str = s.c_str();

        const size_t s_inexact_len = s_inexact.length();
        const char* s_inexact_str = s_inexact.c_str();

        const size_t k_1 = k_-1;

        auto processUnitigMap = [&](const const_UnitigMap<U, G>& um, const size_t pos_s){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) {

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        rpos.add(l_pos_seq);
                    }
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) {

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        rpos.add(l_pos_seq);
                    }
                }
            }
        };

        for (size_t i = 0; i != ((subst || ins) ? 4 : 1); ++i){

            if (ins) {

                for (size_t j = shift; j < s_inexact_len; j += k_) s_inexact[j] = alpha[i];
            }
            else if (subst) {

                for (size_t j = shift; j < s_inexact_len; j += k_) {

                    if (!isDNA(s[j]) || (alpha[i] == s[j])) s_inexact[j] = 'N';
                    else s_inexact[j] = alpha[i];
                }
            } 

            KmerIterator ki_s(s_inexact_str), ki_e;
            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(s_inexact_str, s_inexact_len, k_, g_, RepHash(), true);

            minHashResultIterator<RepHash> it_min, it_min_end;
            minHashResult mhr;

            Minimizer minz;

            pair<size_t, bool> minz_pres = {0xffffffffffffffffULL, true};

            while (ki_s != ki_e) {

                const size_t pos_s = ki_s->second;

                mhi += (pos_s - mhi.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                it_min = *mhi;
                mhr = *it_min;

                // If minimizers of new kmer are different from minimizers of previous kmer
                // or if minimizers are the same but they were present, search them again
                if (minz_pres.second || (mhr.pos != minz_pres.first)){

                    if (mhr.pos != minz_pres.first){

                        minz = Minimizer(s_inexact_str + mhr.pos).rep();
                        minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                        for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                            mhr = *it_min;
                            minz = Minimizer(s_inexact_str + mhr.pos).rep();
                            minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                        }
                    }

                    if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                        const size_t shift_pos_seq = (pos_s / k_) + (pos_s % k_ > shift);
                        const size_t l_pos_s = pos_s - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                        if ((l_pos_s + k_1 < s_len) && isDNA(s_str[l_pos_s]) && isDNA(s_str[l_pos_s + k_1]) && (!rpos.contains(l_pos_s) && (us_pos_km.find({l_pos_s, ki_s->first}) == us_pos_km.end()))) {

                            const const_UnitigMap<U, G> um = findUnitig(s_inexact_str, pos_s, s_inexact_len, mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_s);

                                if (rpos.cardinality() >= nb_km_min) return;

                                ki_s += um.len - 1;
                            }
                        }
                    }
                }

                ++ki_s;
            }
        }
    };

    if (exact){

        for (KmerIterator ki_s(s.c_str()), ki_e; ki_s != ki_e; ++ki_s) {

            const size_t pos_s = ki_s->second;
            const const_UnitigMap<U, G> um = findUnitig(s.c_str(), pos_s, s.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + j - um.dist, um.getKmerMapping(j)});
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                }

                if (v_um.size() >= nb_km_min) return v_um;

                ki_s += um.len - 1;
            }
        }

        for (const auto& pum : v_um) {

            us_pos_km.insert({pum.first, pum.second.getMappedKmer(pum.second.dist)});
            rpos.add(pum.first);
        }
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            s_inexact = s;

            worker_func(true, false, false, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, true, false, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (deletion && (s.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, false, true, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    return v_um;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::searchMinRatioKmer(const vector<string>& query_filenames, const string& out_filename_prefix,
                                            const double min_ratio_kmers,
                                            const bool inexact_search, const bool files_as_queries,
                                            const size_t nb_threads, const size_t verbose) const {

    const string out_tmp = out_filename_prefix + ".tsv";

    {
        FILE* fp_tmp = fopen(out_tmp.c_str(), "w");

        if (fp_tmp == NULL) {

            cerr << "CompactedDBG::searchMinRatioKmer(): Could not open file " << out_tmp << " for writing." << endl;
            return false;
        }
        else {

            fclose(fp_tmp);

            if (std::remove(out_tmp.c_str()) != 0) {

                cerr << "CompactedDBG::searchMinRatioKmer(): Could not remove temporary file " << out_tmp << endl;
            }
        }
    }

    ofstream outfile;
    ostream out(0);

    outfile.open(out_tmp.c_str());
    out.rdbuf(outfile.rdbuf());

    const bool ret = this->searchMinRatioKmer(  query_filenames, out, min_ratio_kmers,
                                                inexact_search, files_as_queries, nb_threads, verbose);

    outfile.close();

    return ret;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::search(const vector<string>& query_filenames, const string& out_filename_prefix,
                                const bool found_km_ratio_out, const bool inexact_search,
                                const bool files_as_queries, const size_t nb_threads, const bool verbose) const {

    const string out_tmp = out_filename_prefix + ".tsv";

    {
        FILE* fp_tmp = fopen(out_tmp.c_str(), "w");

        if (fp_tmp == NULL) {

            cerr << "CompactedDBG::search(): Could not open file " << out_tmp << " for writing." << endl;
            return false;
        }
        else {

            fclose(fp_tmp);

            if (std::remove(out_tmp.c_str()) != 0) {

                cerr << "CompactedDBG::search(): Could not remove temporary file " << out_tmp << endl;
            }
        }
    }

    ofstream outfile;
    ostream out(0);

    outfile.open(out_tmp.c_str());
    out.rdbuf(outfile.rdbuf());

    const bool ret = this->search(  query_filenames, out, found_km_ratio_out,
                                    inexact_search, files_as_queries, nb_threads, verbose);

    outfile.close();

    return ret;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::searchMinRatioKmer(const vector<string>& query_filenames, ostream& out, const double min_ratio_kmers,
                                            const bool inexact_search, const bool files_as_queries,
                                            const size_t nb_threads, const size_t verbose) const {

     if (invalid){

        cerr << "CompactedDBG::searchMinRatioKmer(): Graph is invalid and cannot be searched" << endl;
        return false;
    }

    if (nb_threads > std::thread::hardware_concurrency()){

        cerr << "CompactedDBG::search(): Number of threads cannot be greater than or equal to " << std::thread::hardware_concurrency() << "." << endl;
        return false;
    }

    if (nb_threads <= 0){

        cerr << "CompactedDBG::searchMinRatioKmer(): Number of threads cannot be less than or equal to 0." << endl;
        return false;
    }

    if (min_ratio_kmers <= 0.0){

        cerr << "CompactedDBG::searchMinRatioKmer(): Ratio of k-mers is less than or equal to 0.0." << endl;
        return false;
    }

    if (min_ratio_kmers > 1.0){

        cerr << "CompactedDBG::searchMinRatioKmer(): Ratio of k-mers is greater than 1.0." << endl;
        return false;
    }

    if (out.fail()) {

        cerr << "CompactedDBG::searchMinRatioKmer(): Output stream is in a failed state and cannot be written to." << endl;
        return false;
    }

    if (verbose) cout << "CompactedDBG::searchMinRatioKmer(): Querying graph." << endl;

    const CompactedDBG<U, G>& dbg = *this;

    string s;

    bool write_success = true;
    bool query_success = true;

    size_t file_id = 0;
    size_t prev_file_id = 0xffffffffffffffffULL; // Please don't input 2^64-1 files :D

    const size_t thread_seq_buf_sz = BUFFER_SIZE;

    const double ratio = files_as_queries ? 1.0 : min_ratio_kmers;

    FileParser fp(query_filenames);

    const char query_pres[3] = {'\t', '1', '\n'};
    const char query_abs[3] = {'\t', '0', '\n'};

    const size_t l_query_res = 3;

    if (write_success) {

        out << "query_name\tpresence_query\n"; // Write header to TSV file

        write_success = (write_success && !out.fail());
    }

    if (write_success) {

        if (nb_threads == 1){

            const char* query_name = nullptr;

            char* buffer_res = new char[thread_seq_buf_sz];

            size_t pos_buffer_out = 0;

            size_t nb_queries_found = 0;
            size_t nb_queries_processed = 0;

            size_t nb_km_found = 0;
            size_t nb_km_query = 0;

            auto writeBinaryOutput = [&]() {

                const size_t nb_km_min = max(static_cast<size_t>(1), static_cast<size_t>(round(static_cast<double>(nb_km_query) * min_ratio_kmers)));
                const size_t len_query_name = strlen(query_name);

                const bool is_found = (nb_km_found >= nb_km_min);

                if (pos_buffer_out + len_query_name + l_query_res >= thread_seq_buf_sz){ // If next result cannot fit in the buffer

                    out.write(buffer_res, pos_buffer_out); // Write result buffer

                    pos_buffer_out = 0; // Reset position to 0;
                    write_success = (write_success && !out.fail());
                }

                // Copy new result to buffer
                std::memcpy(buffer_res + pos_buffer_out, query_name, len_query_name * sizeof(char));

                if (is_found){

                    std::memcpy(buffer_res + pos_buffer_out + len_query_name, query_pres, l_query_res * sizeof(char));

                    ++nb_queries_found;
                }
                else std::memcpy(buffer_res + pos_buffer_out + len_query_name, query_abs, l_query_res * sizeof(char));

                pos_buffer_out += len_query_name + l_query_res;
            };

            while (write_success && query_success && fp.read(s, file_id)){

                if (files_as_queries) {

                    if (file_id != prev_file_id) {

                        if (prev_file_id != 0xffffffffffffffffULL) { // Push results to buffer, write buffer if overflow

                            writeBinaryOutput();

                            ++nb_queries_processed;
                        }

                        query_name = query_filenames[file_id].c_str(); // Query name is the filename

                        nb_km_found = 0;
                        nb_km_query = 0;
                    }

                    nb_km_query += s.length() - k_ + 1;
                }
                else {

                    // Push results to buffer, write buffer if overflow
                    if (prev_file_id != 0xffffffffffffffffULL) {

                        writeBinaryOutput();

                        ++nb_queries_processed;
                    }

                    query_name = fp.getNameString(); // Query name is the record name

                    nb_km_query = s.length() - k_ + 1;
                    nb_km_found = 0;
                }

                for (auto& c : s) c &= 0xDF; // Set all characters in uppercase

                const vector<pair<size_t, const_UnitigMap<U, G>>> v = dbg.searchSequence(   s, true, inexact_search, inexact_search,
                                                                                            inexact_search, ratio, true);

                if (inexact_search){

                    Roaring r;

                    for (const auto& p : v) r.add(p.first);

                    nb_km_found += r.cardinality();
                }
                else nb_km_found += v.size();

                prev_file_id = file_id;
            }

            // Flush rest of buffer result to final output
            if (write_success && (prev_file_id != 0xffffffffffffffffULL)) {

                writeBinaryOutput();

                ++nb_queries_processed;

                if (write_success && (pos_buffer_out > 0)) {

                    out.write(buffer_res, pos_buffer_out);

                    write_success = (write_success && !out.fail());
                }
            }

            delete[] buffer_res;

            if (write_success && verbose) {

                cout << "CompactedDBG::searchMinRatioKmer(): Processed " << nb_queries_processed << " queries. " << endl;
                cout << "CompactedDBG::searchMinRatioKmer(): Found " << nb_queries_found << " queries. " << endl;
            }
        }
        else {

            struct ResultFileQuery {

                size_t nb_km_found;
                size_t nb_km_queries;
                size_t nb_queries;

                bool is_read;

                ResultFileQuery() : nb_km_found(0), nb_km_queries(0), nb_queries(0), is_read(false) {}
            };

            bool stop = false;

            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mtx_files_in, mtx_file_out;

            std::atomic<size_t> nb_queries_found;
            std::atomic<size_t> nb_queries_processed;

            unordered_map<size_t, ResultFileQuery> um_file_id;

            nb_queries_found = 0;
            nb_queries_processed = 0;

            auto writeBinaryOutput = [&](   const string& query_name,
                                            const size_t nb_km_found, const size_t nb_km_query,
                                            size_t& pos_buffer_out, char* buffer_res) {

                const size_t nb_km_min = max(static_cast<size_t>(1), static_cast<size_t>(round(static_cast<double>(nb_km_query) * min_ratio_kmers)));
                const size_t len_query_name = query_name.length();

                const bool is_found = (nb_km_found >= nb_km_min);

                if (pos_buffer_out + len_query_name + l_query_res >= thread_seq_buf_sz){ // If next result cannot fit in the buffer

                    unique_lock<mutex> lock(mtx_file_out); // Get the output lock

                    out.write(buffer_res, pos_buffer_out); // Write result buffer

                    pos_buffer_out = 0; // Reset position to 0;
                    write_success = (write_success && !out.fail());
                }

                // Copy new result to buffer
                std::memcpy(buffer_res + pos_buffer_out, query_name.c_str(), len_query_name * sizeof(char));

                if (is_found){

                    std::memcpy(buffer_res + pos_buffer_out + len_query_name, query_pres, l_query_res * sizeof(char));

                    ++nb_queries_found;
                }
                else std::memcpy(buffer_res + pos_buffer_out + len_query_name, query_abs, l_query_res * sizeof(char));

                pos_buffer_out += len_query_name + l_query_res;
            };

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        size_t pos_buffer_out = 0;

                        char* buffer_res = new char[thread_seq_buf_sz];

                        vector<string> buffer_seq;
                        vector<string> buffer_name;
                        vector<pair<size_t, pair<size_t, size_t>>> buffer_file_id;

                        vector<pair<size_t, ResultFileQuery>> v_res_to_write;

                        while (true) {

                            bool l_stop;

                            {
                                size_t buffer_sz = 0;

                                unique_lock<mutex> lock(mtx_files_in);

                                l_stop = stop;

                                if (files_as_queries) {

                                    // Process results from previous search for this thread
                                    // If all queries have completed for this file, push result to buffer 
                                    for (const auto& p : buffer_file_id) {

                                        typename unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(p.first);

                                        if (it_um_file_id == um_file_id.end()) {

                                            query_success = false;
                                            break;
                                        }
                                        else {

                                            it_um_file_id->second.nb_queries -= 1;

                                            it_um_file_id->second.nb_km_found += p.second.first;
                                            it_um_file_id->second.nb_km_queries += p.second.second;

                                            if (it_um_file_id->second.is_read && (it_um_file_id->second.nb_queries == 0)) { // All queries for this file have been processed

                                                v_res_to_write.push_back(*it_um_file_id);
                                                um_file_id.erase(it_um_file_id);

                                                ++nb_queries_processed;
                                            }
                                        }
                                    }
                                }

                                if (query_success) {

                                    // Clear buffers for next round
                                    buffer_seq.clear();
                                    buffer_name.clear();
                                    buffer_file_id.clear();

                                    while (buffer_sz < thread_seq_buf_sz){

                                        stop = !fp.read(s, file_id);

                                        if (!stop) {

                                            buffer_sz += s.length();

                                            buffer_seq.push_back(std::move(s));

                                            if (files_as_queries) buffer_file_id.push_back(pair<size_t, pair<size_t, size_t>>(file_id, pair<size_t, size_t>(0, 0)));
                                            else buffer_name.push_back(string(fp.getNameString()));
                                        }
                                        else break;
                                    }

                                    if (files_as_queries) {

                                        for (const auto p : buffer_file_id) {

                                            pair<typename unordered_map<size_t, ResultFileQuery>::iterator, bool> p_it_um_file_id = um_file_id.insert(pair<size_t, ResultFileQuery>(p.first, ResultFileQuery()));

                                            p_it_um_file_id.first->second.nb_queries += 1;

                                            if ((p.first != prev_file_id) && (prev_file_id != 0xffffffffffffffffULL)) {

                                                typename unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(prev_file_id);

                                                if (it_um_file_id == um_file_id.end()) {

                                                    query_success = false;
                                                    break;
                                                }
                                                else it_um_file_id->second.is_read = true;
                                            }

                                            prev_file_id = p.first;
                                        }

                                        // This thread is the last one reading from input file(s), make sure we annotate last query file as fully read
                                        if (query_success && stop && !l_stop && (prev_file_id != 0xffffffffffffffffULL)) {

                                            typename unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(prev_file_id);

                                            if (it_um_file_id == um_file_id.end()) query_success = false;
                                            else it_um_file_id->second.is_read = true;
                                        }
                                    }
                                }
                            }

                            if (!v_res_to_write.empty()) { // Write results to output if any result available in buffer

                                for (const auto& p : v_res_to_write) writeBinaryOutput(query_filenames[p.first], p.second.nb_km_found, p.second.nb_km_queries, pos_buffer_out, buffer_res);

                                v_res_to_write.clear();
                            }

                            if (l_stop) break;

                            for (size_t i = 0; i < buffer_seq.size(); ++i){

                                const size_t nb_km_query = buffer_seq[i].length() - k_ + 1;

                                size_t nb_km_found = 0;

                                for (auto& c : buffer_seq[i]) c &= 0xDF;

                                const vector<pair<size_t, const_UnitigMap<U, G>>> v = dbg.searchSequence(   buffer_seq[i], true, inexact_search, inexact_search,
                                                                                                            inexact_search, ratio, true);

                                if (inexact_search){

                                    Roaring r;

                                    for (const auto& p : v) r.add(p.first);

                                    nb_km_found = r.cardinality();
                                }
                                else nb_km_found = v.size();

                                if (files_as_queries) {

                                    buffer_file_id[i].second.first = nb_km_found;
                                    buffer_file_id[i].second.second = nb_km_query;
                                }
                                else {

                                    writeBinaryOutput(buffer_name[i], nb_km_found, nb_km_query, pos_buffer_out, buffer_res);

                                    ++nb_queries_processed;
                                }
                            }
                        }

                        if (write_success && (pos_buffer_out > 0)) { // Flush unresult written to final output

                            unique_lock<mutex> lock(mtx_file_out);

                            out.write(buffer_res, pos_buffer_out);

                            write_success = (write_success && !out.fail());
                        }

                        delete[] buffer_res;
                    }
                );
            }

            for (auto& t : workers) t.join();

            if (files_as_queries && !um_file_id.empty()) query_success = false;

            if (write_success && query_success && verbose) {

                cout << "CompactedDBG::searchMinRatioKmer(): Processed " << nb_queries_processed << " queries. " << endl;
                cout << "CompactedDBG::searchMinRatioKmer(): Found " << nb_queries_found << " queries. " << endl;
            }
        }
    }

    fp.close();

    if (!query_success) cerr << "CompactedDBG::searchMinRatioKmer(): Unexpected error encountered. Please file an issue. Operation aborted." << endl;
    if (!write_success) cerr << "CompactedDBG::searchMinRatioKmer(): Output stream is in a failed state and cannot be written to. Operation aborted." << endl;

    return query_success && write_success;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::search(const vector<string>& query_filenames, ostream& out,
                                const bool found_km_ratio_out, const bool inexact_search,
                                const bool files_as_queries, const size_t nb_threads, const bool verbose) const {

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

    if (out.fail()) {

        cerr << "CompactedDBG::search(): Output stream is in a failed state and cannot be written to." << endl;
        return false;
    }

    if (verbose) cout << "CompactedDBG::search(): Querying graph." << endl;

    const CompactedDBG<U, G>& dbg = *this;

    string s;

    bool write_success = true;
    bool query_success = true;

    size_t file_id = 0;
    size_t prev_file_id = 0xffffffffffffffffULL; // Please don't input 2^64-1 files :D

    const size_t thread_seq_buf_sz = BUFFER_SIZE;

    FileParser fp(query_filenames);

    const char query_pres[3] = {'\t', '1', '\n'};
    const char query_abs[3] = {'\t', '0', '\n'};

    const size_t l_query_res = 3;

    if (write_success) {

        // Write header to TSV file
        if (found_km_ratio_out) out << "query_name\tratio_found_kmers\n";
        else out << "query_name\tnb_found_kmers\n";

        write_success = (write_success && !out.fail());
    }

    if (write_success) {

        if (nb_threads == 1){

            const char* query_name = nullptr;

            char* buffer_res = new char[thread_seq_buf_sz];

            size_t pos_buffer_out = 0;

            size_t nb_queries_processed = 0;

            size_t nb_km_found = 0;
            size_t nb_km_query = 0;

            auto writeQuantifiedOutput = [&]() {

                const string nb_found_str = to_string(found_km_ratio_out ? (static_cast<double>(nb_km_found) / static_cast<double>(nb_km_query)) : nb_km_found);

                const size_t len_nb_found_str = nb_found_str.length();
                const size_t len_query_name = strlen(query_name);

                if (pos_buffer_out + len_query_name + len_nb_found_str + 2 >= thread_seq_buf_sz){ // If next result cannot fit in the buffer

                    out.write(buffer_res, pos_buffer_out); // Write result buffer

                    pos_buffer_out = 0; // Reset position to 0;
                    write_success = (write_success && !out.fail());
                }

                // Add query name and tabulation to buffer
                {
                    std::memcpy(buffer_res + pos_buffer_out, query_name, len_query_name * sizeof(char));

                    buffer_res[pos_buffer_out + len_query_name] = '\t';

                    pos_buffer_out += len_query_name + 1;
                }

                // Add number of found km and end line character to buffer
                {
                    std::memcpy(buffer_res + pos_buffer_out, nb_found_str.c_str(), len_nb_found_str * sizeof(char));

                    buffer_res[pos_buffer_out + len_nb_found_str] = '\n';

                    pos_buffer_out += len_nb_found_str + 1;
                }
            };

            while (write_success && query_success && fp.read(s, file_id)){

                if (files_as_queries) {

                    if (file_id != prev_file_id) {

                        if (prev_file_id != 0xffffffffffffffffULL) { // Push results to buffer, write buffer if overflow

                            writeQuantifiedOutput();

                            ++nb_queries_processed;
                        }

                        query_name = query_filenames[file_id].c_str(); // Query name is the filename

                        nb_km_found = 0;
                        nb_km_query = 0;
                    }

                    nb_km_query += s.length() - k_ + 1;
                }
                else {

                    // Push results to buffer, write buffer if overflow
                    if (prev_file_id != 0xffffffffffffffffULL) {

                        writeQuantifiedOutput();

                        ++nb_queries_processed;
                    }

                    query_name = fp.getNameString(); // Query name is the record name

                    nb_km_query = s.length() - k_ + 1;
                    nb_km_found = 0;
                }

                for (auto& c : s) c &= 0xDF; // Set all characters in uppercase

                const vector<pair<size_t, const_UnitigMap<U, G>>> v = dbg.searchSequence(   s, true, inexact_search, inexact_search,
                                                                                            inexact_search, 1.0, true);

                if (inexact_search){

                    Roaring r;

                    for (const auto& p : v) r.add(p.first);

                    nb_km_found += r.cardinality();
                }
                else nb_km_found += v.size();

                prev_file_id = file_id;
            }

            // Flush rest of buffer result to final output
            if (write_success && (prev_file_id != 0xffffffffffffffffULL)) {

                writeQuantifiedOutput();

                ++nb_queries_processed;

                if (write_success && (pos_buffer_out > 0)) {

                    out.write(buffer_res, pos_buffer_out);

                    write_success = (write_success && !out.fail());
                }
            }

            delete[] buffer_res;

            if (write_success && verbose) cout << "CompactedDBG::search(): Processed " << nb_queries_processed << " queries. " << endl;
        }
        else {

            struct ResultFileQuery {

                size_t nb_km_found;
                size_t nb_km_queries;
                size_t nb_queries;

                bool is_read;

                ResultFileQuery() : nb_km_found(0), nb_km_queries(0), nb_queries(0), is_read(false) {}
            };

            bool stop = false;

            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mtx_files_in, mtx_file_out;

            std::atomic<size_t> nb_queries_processed;

            unordered_map<size_t, ResultFileQuery> um_file_id;

            nb_queries_processed = 0;

            auto writeQuantifiedOutput = [&](   const string& query_name,
                                                const size_t nb_km_found, const size_t nb_km_query,
                                                size_t& pos_buffer_out, char* buffer_res) {

                const string nb_found_str = to_string(found_km_ratio_out ? (static_cast<double>(nb_km_found) / static_cast<double>(nb_km_query)) : nb_km_found);

                const size_t len_nb_found_str = nb_found_str.length();
                const size_t len_query_name = query_name.length();

                if (pos_buffer_out + len_query_name + len_nb_found_str + 2 >= thread_seq_buf_sz){ // If next result cannot fit in the buffer

                    unique_lock<mutex> lock(mtx_file_out); // Get the output lock

                    out.write(buffer_res, pos_buffer_out); // Write result buffer

                    pos_buffer_out = 0; // Reset position to 0;
                    write_success = (write_success && !out.fail());
                }

                // Add query name and tabulation to buffer
                {
                    std::memcpy(buffer_res + pos_buffer_out, query_name.c_str(), len_query_name * sizeof(char));

                    buffer_res[pos_buffer_out + len_query_name] = '\t';

                    pos_buffer_out += len_query_name + 1;
                }

                // Add number of found km and end line character to buffer
                {
                    std::memcpy(buffer_res + pos_buffer_out, nb_found_str.c_str(), len_nb_found_str * sizeof(char));

                    buffer_res[pos_buffer_out + len_nb_found_str] = '\n';

                    pos_buffer_out += len_nb_found_str + 1;
                }
            };

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        size_t pos_buffer_out = 0;

                        char* buffer_res = new char[thread_seq_buf_sz];

                        vector<string> buffer_seq;
                        vector<string> buffer_name;
                        vector<pair<size_t, pair<size_t, size_t>>> buffer_file_id;

                        vector<pair<size_t, ResultFileQuery>> v_res_to_write;

                        while (true) {

                            bool l_stop;

                            {
                                size_t buffer_sz = 0;

                                unique_lock<mutex> lock(mtx_files_in);

                                l_stop = stop;

                                if (files_as_queries) {

                                    // Process results from previous search for this thread
                                    // If all queries have completed for this file, push result to buffer 
                                    for (const auto& p : buffer_file_id) {

                                        typename unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(p.first);

                                        if (it_um_file_id == um_file_id.end()) {

                                            query_success = false;
                                            break;
                                        }
                                        else {

                                            it_um_file_id->second.nb_queries -= 1;

                                            it_um_file_id->second.nb_km_found += p.second.first;
                                            it_um_file_id->second.nb_km_queries += p.second.second;

                                            if (it_um_file_id->second.is_read && (it_um_file_id->second.nb_queries == 0)) { // All queries for this file have been processed

                                                v_res_to_write.push_back(*it_um_file_id);
                                                um_file_id.erase(it_um_file_id);

                                                ++nb_queries_processed;
                                            }
                                        }
                                    }
                                }

                                if (query_success) {

                                    // Clear buffers for next round
                                    buffer_seq.clear();
                                    buffer_name.clear();
                                    buffer_file_id.clear();

                                    while (buffer_sz < thread_seq_buf_sz){

                                        stop = !fp.read(s, file_id);

                                        if (!stop) {

                                            buffer_sz += s.length();

                                            buffer_seq.push_back(std::move(s));

                                            if (files_as_queries) buffer_file_id.push_back(pair<size_t, pair<size_t, size_t>>(file_id, pair<size_t, size_t>(0, 0)));
                                            else buffer_name.push_back(string(fp.getNameString()));
                                        }
                                        else break;
                                    }

                                    if (files_as_queries) {

                                        for (const auto p : buffer_file_id) {

                                            pair<typename unordered_map<size_t, ResultFileQuery>::iterator, bool> p_it_um_file_id = um_file_id.insert(pair<size_t, ResultFileQuery>(p.first, ResultFileQuery()));

                                            p_it_um_file_id.first->second.nb_queries += 1;

                                            if ((p.first != prev_file_id) && (prev_file_id != 0xffffffffffffffffULL)) {

                                                typename unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(prev_file_id);

                                                if (it_um_file_id == um_file_id.end()) {

                                                    query_success = false;
                                                    break;
                                                }
                                                else it_um_file_id->second.is_read = true;
                                            }

                                            prev_file_id = p.first;
                                        }

                                        // This thread is the last one reading from input file(s), make sure we annotate last query file as fully read
                                        if (query_success && stop && !l_stop && (prev_file_id != 0xffffffffffffffffULL)) {

                                            typename unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(prev_file_id);

                                            if (it_um_file_id == um_file_id.end()) query_success = false;
                                            else it_um_file_id->second.is_read = true;
                                        }
                                    }
                                }
                            }

                            if (!v_res_to_write.empty()) { // Write results to output if any result available in buffer

                                for (const auto& p : v_res_to_write) writeQuantifiedOutput(query_filenames[p.first], p.second.nb_km_found, p.second.nb_km_queries, pos_buffer_out, buffer_res);

                                v_res_to_write.clear();
                            }

                            if (l_stop) break;

                            for (size_t i = 0; i < buffer_seq.size(); ++i){

                                const size_t nb_km_query = buffer_seq[i].length() - k_ + 1;

                                size_t nb_km_found = 0;

                                for (auto& c : buffer_seq[i]) c &= 0xDF;

                                const vector<pair<size_t, const_UnitigMap<U, G>>> v = dbg.searchSequence(   buffer_seq[i], true, inexact_search, inexact_search,
                                                                                                            inexact_search, 1.0, true);

                                if (inexact_search){

                                    Roaring r;

                                    for (const auto& p : v) r.add(p.first);

                                    nb_km_found = r.cardinality();
                                }
                                else nb_km_found = v.size();

                                if (files_as_queries) {

                                    buffer_file_id[i].second.first = nb_km_found;
                                    buffer_file_id[i].second.second = nb_km_query;
                                }
                                else {

                                    writeQuantifiedOutput(buffer_name[i], nb_km_found, nb_km_query, pos_buffer_out, buffer_res);

                                    ++nb_queries_processed;
                                }
                            }
                        }

                        if (write_success && (pos_buffer_out > 0)) { // Flush unresult written to final output

                            unique_lock<mutex> lock(mtx_file_out);

                            out.write(buffer_res, pos_buffer_out);

                            write_success = (write_success && !out.fail());
                        }

                        delete[] buffer_res;
                    }
                );
            }

            for (auto& t : workers) t.join();

            if (files_as_queries && !um_file_id.empty()) query_success = false;

            if (write_success && query_success && verbose) cout << "CompactedDBG::search(): Processed " << nb_queries_processed << " queries. " << endl;
        }
    }

    fp.close();

    if (!query_success) cerr << "CompactedDBG::search(): Unexpected error encountered. Please file an issue. Operation aborted." << endl;
    if (!write_success) cerr << "CompactedDBG::search(): Output stream is in a failed state and cannot be written to. Operation aborted." << endl;

    return query_success && write_success;
}

#endif
