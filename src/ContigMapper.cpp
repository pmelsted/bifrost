#include "ContigMapper.hpp"
#include "CompressedSequence.hpp"
#include "KmerIterator.hpp"

#include <string>
#include <iterator>
#include <algorithm>
#include <fstream>

#include <tuple>
// for debugging
#include <iostream>

size_t stringMatch(const string& a, const string& b, size_t pos) {

    return distance(a.begin(), mismatch(a.begin(), a.end(), b.begin() + pos).first);
}

// use: delete cm
// pre:
// post: all memory allocated has been released
ContigMapper::~ContigMapper() {

    // we do not own bf pointer
    for (auto contig : v_contigs) delete contig;
}

// use:  ContigMapper(sz)
// pre:  sz >= 0
// post: new contigmapper object
ContigMapper::ContigMapper() : bf(NULL) {}

ContigMapper::ContigMapper(const char* fp_min_name) : bf(NULL) {

    FILE* fp = fopen(fp_min_name, "rb");
    assert(fp != NULL);

    size_t hmap_min_sz;

    if (fread(&hmap_min_sz, sizeof(size_t), 1, fp) != 1){

        cerr << "ContigMapper(): cannot read minimizers from disk" << endl;
        exit(1);
    }

    hmap_min_contigs.init_table(hmap_min_sz);

    Minimizer minz;

    for (size_t i = 0; i < hmap_min_sz; i++){

        if (minz.read(fp)) hmap_min_contigs.insert(make_pair(minz, tiny_vector<size_t,tiny_vector_sz>()));
        else{

            cerr << "ContigMapper(): cannot read minimizers from disk" << endl;
            exit(1);
        }
    }

    fclose(fp);
}

// user: i = cm.contigCount()
// pre:
// post: i is the number of contigs in the mapper
size_t ContigMapper::contigCount() const {

    return v_contigs.size() + v_kmers.size() + h_kmers_ccov.size();
}

void ContigMapper::printContigCount() const {

    cerr << "v_contigs.size(): " << v_contigs.size() << endl;
    cerr << "v_kmers.size(): " << v_kmers.size() << endl;
    cerr << "h_kmers_ccov.size(): " << h_kmers_ccov.size() << endl;
}

// use:  cm.mapBloomFilter(bf)
// pre:  bf != null
// post: uses the bloom filter bf to map reads
void ContigMapper::mapBloomFilter(const BlockedBloomFilter *bf) { this->bf = bf; }

// use:  cm.mapRead(km,pos,cc)
// pre:  cc is a reference to a current contig in cm, km maps to cc
// post: the coverage information in cc has been updated
void ContigMapper::mapRead(const ContigMap& cc) {

    if (cc.isEmpty) return; // nothing maps, move on

    if (cc.isShort) v_kmers[cc.pos_contig].second.cover(cc.dist, cc.dist + cc.len - 1);
    else if (cc.isAbundant) h_kmers_ccov.find(cc.pos_contig)->second.cover(cc.dist, cc.dist + cc.len - 1);
    else v_contigs[cc.pos_contig]->ccov.cover(cc.dist, cc.dist + cc.len - 1);
}

// use: b = cm.addContig(km,read)
// pre:
// post: either contig string containsin has been added and b == true
//       or it was present and the coverage information was updated, b == false
//       NOT Threadsafe!
//bool ContigMapper::addContigSequence(Kmer km, const string& read, size_t pos, const string& seq) {
bool ContigMapper::addContigSequence(Kmer km, const string& read, size_t pos, const string& seq, vector<Kmer>& l_ignored_km_tip) {

    string s;
    bool selfLoop = false;
    bool isIsolated = false;

    if (!seq.empty()) s = seq;
    else findContigSequence(km, s, selfLoop, isIsolated, l_ignored_km_tip);

    size_t k = Kmer::k;

    if (selfLoop) {

        bool foundAny = false;

        for (KmerIterator it(s.c_str()), it_end; it != it_end; ++it) {

            ContigMap cm = find(it->first);

            if (!cm.isEmpty) {

                mapRead(cm);
                foundAny = true;
            }
        }

        if (!foundAny) {

            addContig(s, s.length() == k ? v_kmers.size() : v_contigs.size());

            for (KmerIterator it(s.c_str()), it_end; it != it_end; ++it) mapRead(find(it->first));
        }

        return true;
    }

    ContigMap cm = findContig(km, read, pos);

    if (cm.isEmpty){

        addContig(s, s.length() == k ? v_kmers.size() : v_contigs.size());
        cm = findContig(km, read, pos);
    }

    mapRead(cm);

    return !cm.isEmpty;
}

// use:  cm.findContigSequence(km, s, selfLoop)
// pre:  km is in the bloom filter
// post: s is the contig containing the kmer km
//       and the first k-mer in s is smaller (wrt. < operator)
//       than the last kmer
//       selfLoop is true of the contig is a loop or hairpin
/*size_t ContigMapper::findContigSequence(Kmer km, string& s, bool& selfLoop, bool& isIsolated) {

    string fw_s;

    Kmer end = km;
    Kmer last = end;
    Kmer twin = km.twin();

    char c;

    size_t j = 0;

    bool has_no_neighbor;

    selfLoop = false;
    isIsolated = false;

    while (fwBfStep(end, end, c, has_no_neighbor)) {

        j++;

        if (end == km) {
            selfLoop = true;
            break;
        }
        else if (end == twin) break;
        else if (end == last.twin()) break;

        fw_s.push_back(c);
        last = end;
    }

    string bw_s;
    Kmer front = km;
    Kmer first = front;

    if (!selfLoop) {

        isIsolated = (j == 0) && has_no_neighbor;
        j = 0;

        while (bwBfStep(front, front, c, has_no_neighbor)) {

            j++;

            if (front == km) {
                selfLoop = true;
                break;
            }
            else if (front == twin) break;
            else if (front == first.twin()) break;

            bw_s.push_back(c);
            first = front;
        }

        if (isIsolated) isIsolated = (j == 0) && has_no_neighbor;

        reverse(bw_s.begin(), bw_s.end());
    }

    s.reserve(Kmer::k + fw_s.size() + bw_s.size());
    s.append(bw_s);

    char tmp[Kmer::MAX_K];
    km.toString(tmp);
    s.append(tmp);

    s.append(fw_s);

    return bw_s.size();
}*/

size_t ContigMapper::findContigSequence(Kmer km, string& s, bool& selfLoop, bool& isIsolated, vector<Kmer>& l_ignored_km_tip) {

    string fw_s;

    Kmer end = km;
    Kmer last = end;
    Kmer twin = km.twin();

    char c;

    size_t j = 0;

    bool has_no_neighbor;

    selfLoop = false;
    isIsolated = false;

    while (fwBfStep(end, end, c, has_no_neighbor, l_ignored_km_tip)) {

        j++;

        if (end == km) {
            selfLoop = true;
            break;
        }
        else if (end == twin) break;
        else if (end == last.twin()) break;

        fw_s.push_back(c);
        last = end;
    }

    string bw_s;
    Kmer front = km;
    Kmer first = front;

    if (!selfLoop) {

        isIsolated = (j == 0) && has_no_neighbor;
        j = 0;

        while (bwBfStep(front, front, c, has_no_neighbor, l_ignored_km_tip)) {

            j++;

            if (front == km) {
                selfLoop = true;
                break;
            }
            else if (front == twin) break;
            else if (front == first.twin()) break;

            bw_s.push_back(c);
            first = front;
        }

        if (isIsolated) isIsolated = (j == 0) && has_no_neighbor;

        reverse(bw_s.begin(), bw_s.end());
    }

    s.reserve(Kmer::k + fw_s.size() + bw_s.size());
    s.append(bw_s);

    char tmp[Kmer::MAX_K];
    km.toString(tmp);
    s.append(tmp);

    s.append(fw_s);

    return bw_s.size();
}

// use:  cc = cm.findContig(km,s,pos)
// pre:  s[pos,pos+k-1] is the kmer km
// post: cc contains either the reference to the contig position
//       or empty if none found
ContigMap ContigMapper::findContig(const Kmer& km, const string& s, size_t pos) const {

    //cerr << "findContig()" << endl;

    assert(bf != NULL);

    // need to check if we find it right away, need to treat this common case
    ContigMap cc = find(km);

    if (!cc.isEmpty && !cc.isShort && !cc.isAbundant){

        const CompressedSequence& seq = v_contigs[cc.pos_contig]->seq;
        size_t km_dist = cc.dist;
        size_t jlen = 0;
        size_t k = Kmer::k;

        if (cc.strand) jlen = seq.jump(s.c_str(), pos, cc.dist, false) - k + 1;
        else {

            jlen = seq.jump(s.c_str(), pos, cc.dist + k - 1, true) - k + 1; // match s_fw to comp(seq)_bw
            km_dist -= jlen - 1;
        }

        //cerr << " exit findContig()" << endl;

        return ContigMap(cc.pos_contig, cc.pos_min, km_dist, jlen, cc.size, false, false, cc.strand);
    }

    //cerr << " exit findContig()" << endl;

    return cc;
}

ContigMap ContigMapper::findContig(const Kmer& km, const string& s, size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h) const {

    //cerr << "findContig() 2" << endl;

    assert(bf != NULL);

    // need to check if we find it right away, need to treat this common case
    ContigMap cc = find(km, it_min_h);

    if (!cc.isEmpty && !cc.isShort && !cc.isAbundant){

        const CompressedSequence& seq = v_contigs[cc.pos_contig]->seq;
        size_t km_dist = cc.dist;
        size_t jlen = 0;
        size_t k = Kmer::k;

        if (cc.strand) jlen = seq.jump(s.c_str(), pos, cc.dist, false) - k + 1;
        else {

            jlen = seq.jump(s.c_str(), pos, cc.dist + k - 1, true) - k + 1; // match s_fw to comp(seq)_bw
            km_dist -= jlen - 1;
        }

        return ContigMap(cc.pos_contig, cc.pos_min, km_dist, jlen, cc.size, false, false, cc.strand);
    }

    return cc;
}

/*ContigMap ContigMapper::findContig(const Kmer& km, const string& s, size_t pos_km_in_s, size_t pos_min_in_s, size_t it_h) {

    assert(bf != NULL);

    // need to check if we find it right away, need to treat this common case
    ContigMap cc = this->find(km, pos_min_in_s - pos_km_in_s, it_h);

    if (!cc.isEmpty && !cc.isShort && !cc.isAbundant){

        const CompressedSequence& seq = v_contigs[cc.pos_contig]->seq;
        size_t km_dist = cc.dist;
        size_t jlen = 0;
        size_t k = Kmer::k;

        if (cc.strand) jlen = seq.jump(s.c_str(), pos_km_in_s, cc.dist, false) - k + 1;
        else {

            jlen = seq.jump(s.c_str(), pos_km_in_s, cc.dist + k - 1, true) - k + 1; // match s_fw to comp(seq)_bw
            km_dist -= jlen - 1;
        }

        return ContigMap(cc.pos_contig, cc.pos_min, km_dist, jlen, cc.size, false, false, cc.strand);
    }

    return cc;
}*/

ContigMap ContigMapper::find(const Kmer& km, bool extremities_only) const {

    //cerr << "find()" << endl;

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort, isAbundant = false;

    size_t contig_id, len, it_h = 0;

    int64_t pos_match;

    int k = Kmer::k, g = Minimizer::g;

    char km_tmp[k + 1];
    km.toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k, k, g, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz = Minimizer(&km_tmp[min_h_res.pos]).rep();
        hmap_min_contigs_t::const_iterator it = hmap_min_contigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_contigs.end()){ // If the minimizer is found

            it_h = it.getHash();

            const tiny_vector<size_t,tiny_vector_sz>& v = it->second;

            it = hmap_min_contigs.end();

            for (auto contig_id_pos : v){ // For each contig associated with the minimizer

                contig_id = contig_id_pos >> 32;

                if (contig_id == RESERVED_ID){

                    if ((contig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()) return ContigMap(it_km.getHash(), it_h, 0, 1, Kmer::k, false, true, km == km_rep);
                    }

                    if ((contig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(&km_tmp[mhr.pos]).rep();
                            it = hmap_min_contigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (contig_id_pos & MASK_CONTIG_TYPE) != 0;
                    contig_id_pos &= MASK_CONTIG_POS;

                    pos_match = contig_id_pos - min_h_res.pos;

                    if (isShort){

                        if (min_h_res.pos == contig_id_pos){

                            if (v_kmers[contig_id].first == km_rep) return ContigMap(contig_id, it_h, 0, 1, k, true, false, true);
                        }
                        else if ((min_h_res.pos == k - contig_id_pos - g) && (v_kmers[contig_id].first == km_rep)){

                            return ContigMap(contig_id, it_h, 0, 1, k, true, false, false);
                        }
                    }
                    else {

                        len = v_contigs[contig_id]->length();

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len - k)) && (v_contigs[contig_id]->seq.getKmer(pos_match) == km)){

                                return ContigMap(contig_id, it_h, pos_match, 1, len, false, false, true);
                            }
                        }
                        else if ((pos_match >= 0) && (pos_match <= len - k) && (v_contigs[contig_id]->seq.getKmer(pos_match) == km)){

                            return ContigMap(contig_id, it_h, pos_match, 1, len, false, false, true);
                        }

                        pos_match = contig_id_pos - k + g + min_h_res.pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len - k)) && (v_contigs[contig_id]->seq.getKmer(pos_match) == km_twin)){

                                return ContigMap(contig_id, it_h, pos_match, 1, len, false, false, false);
                            }
                        }
                        else if ((pos_match >= 0) && (pos_match <= len - k) && (v_contigs[contig_id]->seq.getKmer(pos_match) == km_twin)){

                            return ContigMap(contig_id, it_h, pos_match, 1, len, false, false, false);
                        }
                    }
                }
            }
        }

        it_it_min++;
    }

    return ContigMap(it_h);
}

ContigMap ContigMapper::find(const Kmer& km, const preAllocMinHashIterator<RepHash>& it_min_h) const {

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort, isAbundant = false;

    size_t contig_id, len, it_h = 0;

    int64_t pos_match;

    int k = Kmer::k, g = Minimizer::g;

    preAllocMinHashIterator<RepHash> it_min = preAllocMinHashIterator<RepHash>(it_min_h, k);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz = Minimizer(&it_min.s[min_h_res.pos]).rep();
        hmap_min_contigs_t::const_iterator it = hmap_min_contigs.find(minz); // Look for the minimizer in the hash table

        mhr = min_h_res;

        while (it != hmap_min_contigs.end()){ // If the minimizer is found

            const tiny_vector<size_t,tiny_vector_sz>& v = it->second;

            it_h = it.getHash();
            it = hmap_min_contigs.end();

            for (auto contig_id_pos : v){ // For each contig associated with the minimizer

                contig_id = contig_id_pos >> 32;

                if (contig_id == RESERVED_ID){

                    if ((contig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()){

                            return ContigMap(it_km.getHash(), it_h, 0, 1, Kmer::k, false, true, km == km_rep);
                        }
                    }

                    if ((contig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(&it_min.s[mhr.pos]).rep();
                            it = hmap_min_contigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (contig_id_pos & MASK_CONTIG_TYPE) != 0;
                    contig_id_pos &= MASK_CONTIG_POS;

                    if (isShort){

                        if (min_h_res.pos == contig_id_pos){

                            if (v_kmers[contig_id].first == km_rep){

                                return ContigMap(contig_id, it_h, 0, 1, k, true, false, true);
                            }
                        }
                        else if ((min_h_res.pos == k - contig_id_pos - g) && (v_kmers[contig_id].first == km_rep)){

                            return ContigMap(contig_id, it_h, 0, 1, k, true, false, false);
                        }
                    }
                    else {

                        pos_match = contig_id_pos - min_h_res.pos;
                        len = v_contigs[contig_id]->length();

                        if ((pos_match >= 0) && (pos_match <= len - k) && (v_contigs[contig_id]->seq.getKmer(pos_match) == km)){

                            return ContigMap(contig_id, it_h, pos_match, 1, len, false, false, true);
                        }

                        pos_match = contig_id_pos - k + g + min_h_res.pos;

                        if ((pos_match >= 0) && (pos_match <= len - k) && (v_contigs[contig_id]->seq.getKmer(pos_match) == km_twin)){

                            return ContigMap(contig_id, it_h, pos_match, 1, len, false, false, false);
                        }
                    }
                }
            }
        }

        it_it_min++;
    }

    return ContigMap(it_h);
}

/*ContigMap ContigMapper::find(const Kmer& km, const size_t pos_min, const size_t it_h, bool extremities_only) {

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort, isAbundant = false;

    size_t contig_id, len;

    int64_t pos_match;

    int k = Kmer::k, g = Minimizer::g;

    hmap_min_contigs_t::iterator it = hmap_min_contigs.find(it_h); // Look for the minimizer in the hash table

    if (it != hmap_min_contigs.end()){ // If the minimizer is found

        for (auto contig_id_pos : it->second){ // For each contig associated with the minimizer

            contig_id = contig_id_pos >> 32;

            if (contig_id == RESERVED_ID) isAbundant = true;
            else {

                isShort = (contig_id_pos & MASK_CONTIG_TYPE) != 0;
                contig_id_pos &= MASK_CONTIG_POS;

                pos_match = contig_id_pos - pos_min;

                if (isShort){

                    if (pos_min == contig_id_pos){

                        if (v_kmers[contig_id].first == km_rep) return ContigMap(contig_id, it_h, 0, 1, k, true, false, true);
                    }
                    else if ((pos_min == k - contig_id_pos - g) && (v_kmers[contig_id].first == km_rep)){

                        return ContigMap(contig_id, it_h, 0, 1, k, true, false, false);
                    }
                }
                else {

                    len = v_contigs[contig_id]->length();

                    if (extremities_only){

                        if (((pos_match == 0) || (pos_match == len - k)) && (v_contigs[contig_id]->seq.getKmer(pos_match) == km)){

                            return ContigMap(contig_id, it_h, pos_match, 1, len, false, false, true);
                        }
                    }
                    else if ((pos_match >= 0) && (pos_match <= len - k) && (v_contigs[contig_id]->seq.getKmer(pos_match) == km)){

                        return ContigMap(contig_id, it_h, pos_match, 1, len, false, false, true);
                    }

                    pos_match = contig_id_pos - k + g + pos_min;

                    if (extremities_only){

                        if (((pos_match == 0) || (pos_match == len - k)) && (v_contigs[contig_id]->seq.getKmer(pos_match) == km_twin)){

                            return ContigMap(contig_id, it_h, pos_match, 1, len, false, false, false);
                        }
                    }
                    else if ((pos_match >= 0) && (pos_match <= len - k) && (v_contigs[contig_id]->seq.getKmer(pos_match) == km_twin)){

                        return ContigMap(contig_id, it_h, pos_match, 1, len, false, false, false);
                    }
                }
            }
        }
    }

    if (isAbundant){

        h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

        if (it_km != h_kmers_ccov.end()) return ContigMap(it_km.getHash(), it_h, 0, 1, Kmer::k, false, true, km == km_rep);
    }

    return ContigMap(it_h);
}*/

// use:  b = cm.bwBfStep(km,front,c,deg)
// pre:  km is in the bloom filter
// post: b is true if km is inside a contig, in that
//       case end is the bw link and c is the nucleotide used for the link.
//       if b is false, front and c are not updated
//       if km is an isolated self link (e.g. 'AAA') i.e. end == km then returns false
//       deg is the backwards degree of the front
/*bool ContigMapper::bwBfStep(Kmer km, Kmer& front, char& c, bool& has_no_neighbor) const {

    size_t i, j = -1, k = Kmer::k, g = Minimizer::g;
    size_t nb_neigh = 0;

    char km_tmp[k + 1];

    front.backwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
    std::pair<uint64_t*, uint64_t*> block = bf->getBlock(it_min_h);

    RepHash rep_h(k - 1), rep_h_cpy;
    rep_h.init(km_tmp + 1);

    for (i = 0; (i < 4) && (nb_neigh < 2); ++i) {

        km_tmp[0] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendBW(alpha[i]);

        if (bf->contains(rep_h_cpy.hash(), block)){

            j = i;
            ++nb_neigh;
        }
    }

    if (nb_neigh != 1) {

        has_no_neighbor = (nb_neigh == 0);
        return false;
    }
    else has_no_neighbor = false;

    // only one k-mer in the bw link
    nb_neigh = 0;

    Kmer bw = front.backwardBase(alpha[j]);

    bw.forwardBase('A').toString(km_tmp);

    it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
    block = bf->getBlock(it_min_h);

    rep_h.init(km_tmp);

    for (i = 0; (i < 4) && (nb_neigh < 2); ++i) {

        km_tmp[k - 1] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendFW(alpha[i]);

        if (bf->contains(rep_h_cpy.hash(), block)) ++nb_neigh;
    }

    assert(nb_neigh >= 1);
    if (nb_neigh != 1) return false;

    if (bw != km) {
        // exactly one k-mer in bw, character used is c
        front = bw;
        c = alpha[j];
        return true;
    }

    return false;
}*/

bool ContigMapper::bwBfStep(Kmer km, Kmer& front, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand) const {

    size_t i, j = -1, j_tmp, k = Kmer::k, g = Minimizer::g;
    size_t nb_neigh = 0;

    char km_tmp[k + 1];

    front.backwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
    std::pair<uint64_t*, uint64_t*> block = bf->getBlock(it_min_h);

    RepHash rep_h(k - 1), rep_h_cpy;
    rep_h.init(km_tmp + 1);

    int found_fp_bw = 0;
    Kmer km_fp;

    bool pres_neigh_bw[4] = {false, false, false, false};

    for (i = 0; i < 4; ++i) {

        km_tmp[0] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendBW(alpha[i]);

        if (bf->contains(rep_h_cpy.hash(), block)){

            j = i;
            pres_neigh_bw[i] = true;
            nb_neigh++;

            if (!check_fp_cand && (nb_neigh >= 2)) break;
        }
    }

    if (check_fp_cand && (nb_neigh >= 2)){

        for (i = 0; i < 4; ++i) {

            if (pres_neigh_bw[i]){

                char dummy;
                bool has_no_neighbor_tmp = false;

                km_tmp[0] = alpha[i];

                km_fp = Kmer(km_tmp);

                bwBfStep(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                if (has_no_neighbor_tmp && fwBfStep(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false)) found_fp_bw++;
                else {

                    j_tmp = i;
                    pres_neigh_bw[i] = false;
                }
            }
        }

        if (found_fp_bw != 0){

            if ((nb_neigh - found_fp_bw) != 0){

                j = j_tmp;
                nb_neigh -= found_fp_bw;
            }
            else found_fp_bw = 0;
        }
    }

    if (nb_neigh != 1) {

        has_no_neighbor = (nb_neigh == 0);
        return false;
    }
    else has_no_neighbor = false;

    if (check_fp_cand){

        nb_neigh = 0;

        int found_fp_fw = 0;
        bool pres_neigh_fw[4] = {false, false, false, false};

        Kmer bw = front.backwardBase(alpha[j]);

        bw.forwardBase('A').toString(km_tmp);

        it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
        block = bf->getBlock(it_min_h);

        rep_h.init(km_tmp);

        for (i = 0; i < 4; ++i) {

            km_tmp[k - 1] = alpha[i];

            rep_h_cpy = rep_h;
            rep_h_cpy.extendFW(alpha[i]);

            if (bf->contains(rep_h_cpy.hash(), block)){

                nb_neigh++;
                pres_neigh_fw[i] = true;
            }
        }

        if (nb_neigh >= 2){

            for (i = 0; i < 4; ++i) {

                if (pres_neigh_fw[i]){

                    char dummy;
                    bool add = true, has_no_neighbor_tmp = false;

                    km_tmp[k - 1] = alpha[i];

                    km_fp = Kmer(km_tmp);

                    fwBfStep(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                    if (has_no_neighbor_tmp && bwBfStep(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false)){

                        if (km_fp != km) found_fp_fw++;
                        else {

                            found_fp_fw = 0;
                            break;
                        }
                    }
                    else pres_neigh_fw[i] = false;
                }
            }

            if (found_fp_fw != 0){

                if ((nb_neigh - found_fp_fw) != 0) nb_neigh -= found_fp_fw;
                else found_fp_fw = 0;
            }
        }

        if (nb_neigh != 1) return false;

        if (bw != km) {
            // exactly one k-mer in bw, character used is c

            for (i = 0; (i < 4) && (found_fp_fw != 0); ++i) {

                if (pres_neigh_fw[i]){

                    km_tmp[k - 1] = alpha[i];
                    km_fp = Kmer(km_tmp).rep();

                    l_ignored_km_tip.push_back(km_fp);

                    found_fp_fw--;
                }
            }

            front.backwardBase('A').toString(km_tmp);

            for (i = 0; (i < 4) && (found_fp_bw != 0); ++i) {

                if (pres_neigh_bw[i]){

                    km_tmp[0] = alpha[i];
                    km_fp = Kmer(km_tmp).rep();

                    l_ignored_km_tip.push_back(km_fp);

                    found_fp_bw--;
                }
            }

            front = bw;
            c = alpha[j];

            return true;
        }

        return false;
    }

    return true;
}

// use:  b = cm.fwBfStep(km,end,c,deg)
// pre:  km is in the bloom filter
// post: b is true if km is inside a contig, in that
//       case end is the fw link and c is the nucleotide used for the link.
//       if b is false, end and c are not updated
//       if km is an isolated self link (e.g. 'AAA') i.e. end == km then returns false
//       deg is the degree of the end
/*bool ContigMapper::fwBfStep(Kmer km, Kmer& end, char& c, bool& has_no_neighbor) const {

    size_t i, j = -1, k = Kmer::k, g = Minimizer::g;
    size_t nb_neigh = 0;

    char km_tmp[k + 1];

    end.forwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
    std::pair<uint64_t*, uint64_t*> block = bf->getBlock(it_min_h);

    RepHash rep_h(k - 1), rep_h_cpy;
    rep_h.init(km_tmp);

    for (i = 0; (i < 4) && (nb_neigh < 2); ++i) {

        km_tmp[k - 1] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendFW(alpha[i]);

        if (bf->contains(rep_h_cpy.hash(), block)){

            j = i;
            ++nb_neigh;
        }
    }

    if (nb_neigh != 1) {

        has_no_neighbor = (nb_neigh == 0);
        return false;
    }
    else has_no_neighbor = false;

    // only one k-mer in fw link
    nb_neigh = 0;

    Kmer fw = end.forwardBase(alpha[j]);

    // check bw from fw link
    fw.backwardBase('A').toString(km_tmp);

    it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
    block = bf->getBlock(it_min_h);

    rep_h.init(km_tmp + 1);

    for (i = 0; (i < 4) && (nb_neigh < 2); ++i) {

        km_tmp[0] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendBW(alpha[i]);

        if (bf->contains(rep_h_cpy.hash(), block)) ++nb_neigh;
    }

    assert(nb_neigh >= 1);

    if (nb_neigh != 1) return false;

    if (fw != km) {

        // exactly one k-mer fw, character used is c
        end = fw;
        c = alpha[j];
        return true;
    }

    return false;
}*/

bool ContigMapper::fwBfStep(Kmer km, Kmer& end, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand) const {

    size_t i, j = -1, j_tmp, k = Kmer::k, g = Minimizer::g;
    size_t nb_neigh = 0;

    char km_tmp[k + 1];

    end.forwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
    std::pair<uint64_t*, uint64_t*> block = bf->getBlock(it_min_h);

    RepHash rep_h(k - 1), rep_h_cpy;
    rep_h.init(km_tmp);

    int found_fp_fw = 0;
    Kmer km_fp;

    bool pres_neigh_fw[4] = {false, false, false, false};

    for (i = 0; i < 4; ++i) {

        km_tmp[k - 1] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendFW(alpha[i]);

        if (bf->contains(rep_h_cpy.hash(), block)){

            j = i;
            pres_neigh_fw[i] = true;
            nb_neigh++;

            if (!check_fp_cand && (nb_neigh >= 2)) break;
        }
    }

    if (check_fp_cand && (nb_neigh >= 2)){

        for (i = 0; i < 4; ++i) {

            if (pres_neigh_fw[i]){

                km_tmp[k - 1] = alpha[i];

                char dummy;
                km_fp = Kmer(km_tmp);

                bool has_no_neighbor_tmp = false;

                fwBfStep(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                if (has_no_neighbor_tmp && bwBfStep(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false)) found_fp_fw++;
                else {

                    j_tmp = i;
                    pres_neigh_fw[i] = false;
                }
            }
        }

        if (found_fp_fw != 0){

            if ((nb_neigh - found_fp_fw) != 0){

                j = j_tmp;
                nb_neigh -= found_fp_fw;
            }
            else found_fp_fw = 0;
        }
    }

    if (nb_neigh != 1) {

        has_no_neighbor = (nb_neigh == 0);
        return false;
    }
    else has_no_neighbor = false;

    if (check_fp_cand){

        nb_neigh = 0;

        int found_fp_bw = 0;
        bool pres_neigh_bw[4] = {false, false, false, false};

        Kmer fw = end.forwardBase(alpha[j]);

        // check bw from fw link
        fw.backwardBase('A').toString(km_tmp);

        it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
        block = bf->getBlock(it_min_h);

        rep_h.init(km_tmp + 1);

        for (i = 0; i < 4; ++i) {

            km_tmp[0] = alpha[i];

            rep_h_cpy = rep_h;
            rep_h_cpy.extendBW(alpha[i]);

            if (bf->contains(rep_h_cpy.hash(), block)){

                nb_neigh++;
                pres_neigh_bw[i] = true;
            }
        }

        if (nb_neigh >= 2){

            for (i = 0; i < 4; ++i) {

                if (pres_neigh_bw[i]){

                    char dummy;
                    bool has_no_neighbor_tmp = false;

                    km_tmp[0] = alpha[i];

                    km_fp = Kmer(km_tmp);

                    bwBfStep(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                    if (has_no_neighbor_tmp && fwBfStep(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false)){

                        if (km_fp != km) found_fp_bw++;
                        else {

                            found_fp_bw = 0;
                            break;
                        }
                    }
                    else pres_neigh_bw[i] = false;
                }
            }

            if (found_fp_bw != 0){

                if ((nb_neigh - found_fp_bw) != 0) nb_neigh -= found_fp_bw;
                else found_fp_bw = 0;
            }
        }

        if (nb_neigh != 1) return false;

        if (fw != km) {

            for (i = 0; (i < 4) && (found_fp_bw != 0); ++i) {

                if (pres_neigh_bw[i]){

                    km_tmp[0] = alpha[i];
                    km_fp = Kmer(km_tmp).rep();

                    l_ignored_km_tip.push_back(km_fp);

                    found_fp_bw--;
                }
            }

            end.forwardBase('A').toString(km_tmp);

            for (i = 0; (i < 4) && (found_fp_fw != 0); ++i) {

                if (pres_neigh_fw[i]){

                    km_tmp[k - 1] = alpha[i];
                    km_fp = Kmer(km_tmp).rep();

                    l_ignored_km_tip.push_back(km_fp);

                    found_fp_fw--;
                }
            }

            end = fw;
            c = alpha[j];
            return true;
        }

        return false;
    }

    return true;
}

void ContigMapper::check_fp_tips(KmerHashTable<bool>& ignored_km_tips){

    uint64_t nb_real_short_tips = 0;

    size_t k = Kmer::k, g = Minimizer::g;

    size_t nxt_pos_insert_v_contigs = v_contigs.size();
    size_t v_contigs_sz = v_contigs.size();
    size_t v_kmers_sz = v_kmers.size();

    char km_tmp[k + 1];

    vector<pair<int,int>> sp;

    cerr << "ignored_km_tips.size() = " << ignored_km_tips.size() << endl;

    for (KmerHashTable<bool>::iterator it = ignored_km_tips.begin(); it != ignored_km_tips.end(); it++) {

        Kmer km = it->first;

        ContigMap cm = find(km, true); // Check if the (short) tip actually exists

        if (!cm.isEmpty){ // IF the tip exists

            nb_real_short_tips++;

            bool not_found = true;

            km.backwardBase('A').toString(km_tmp);

            uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            std::pair<uint64_t*, uint64_t*> block = bf->getBlock(it_min_h);

            RepHash rep_h(k - 1), rep_h_cpy;
            rep_h.init(km_tmp + 1);

            for (size_t i = 0; (i < 4) && not_found; ++i) {

                km_tmp[0] = alpha[i];

                rep_h_cpy = rep_h;
                rep_h_cpy.extendBW(alpha[i]);

                if (bf->contains(rep_h_cpy.hash(), block)){ // Query for all of its predecessors

                    ContigMap cm_bw = find(Kmer(km_tmp));

                    if (!cm_bw.isEmpty && !cm_bw.isAbundant && !cm_bw.isShort){

                        if (cm_bw.strand) cm_bw.dist++;

                        if ((cm_bw.dist != 0) && (cm_bw.dist != cm_bw.size - k + 1)){

                            sp.push_back(make_pair(0, cm_bw.dist));
                            sp.push_back(make_pair(cm_bw.dist, cm_bw.size - k + 1));

                            splitContig(cm_bw.pos_contig, nxt_pos_insert_v_contigs, v_contigs_sz, v_kmers_sz, sp);

                            sp.clear();
                        }

                        not_found = false;
                    }
                }
            }

            km.forwardBase('A').toString(km_tmp);

            it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            block = bf->getBlock(it_min_h);

            rep_h.init(km_tmp);

            for (size_t i = 0; (i < 4) && not_found; ++i) {

                km_tmp[k - 1] = alpha[i];

                rep_h_cpy = rep_h;
                rep_h_cpy.extendFW(alpha[i]);

                if (bf->contains(rep_h_cpy.hash(), block)){

                    ContigMap cm_fw = find(Kmer(km_tmp));

                    if (!cm_fw.isEmpty && !cm_fw.isAbundant && !cm_fw.isShort){

                        if (!cm_fw.strand) cm_fw.dist++;

                        if ((cm_fw.dist != 0) && (cm_fw.dist != cm_fw.size - k + 1)){

                            sp.push_back(make_pair(0, cm_fw.dist));
                            sp.push_back(make_pair(cm_fw.dist, cm_fw.size - k + 1));

                            splitContig(cm_fw.pos_contig, nxt_pos_insert_v_contigs, v_contigs_sz, v_kmers_sz, sp);

                            sp.clear();
                        }

                        not_found = false;
                    }
                }
            }
        }
    }

    if (nxt_pos_insert_v_contigs < v_contigs.size()) v_contigs.resize(nxt_pos_insert_v_contigs);
    if (v_kmers_sz < v_kmers.size()) v_kmers.resize(v_kmers_sz);

    cerr << "Number of real short tips = " << nb_real_short_tips << endl;
}

bool ContigMapper::isMinPresent(const Minimizer& minz) const{

    return (hmap_min_contigs.find(minz) != hmap_min_contigs.end());
}

bool ContigMapper::addContig(const string& str_contig, const size_t id_contig){

    int k = Kmer::k, g = Minimizer::g;

    size_t len = str_contig.size();
    size_t pos_id_contig = id_contig << 32;

    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    bool isShort = false, isAbundant = false, isForbidden = false;

    char* c_str = const_cast<char*>(str_contig.c_str());

    char km_tmp[k + 1];

    Kmer km_rep;

    if (len == k){ // Contig to add is short, maybe abundant as well

        isShort = true;

        pos_id_contig |= MASK_CONTIG_TYPE;

        km_rep = Kmer(c_str).rep();
        km_rep.toString(km_tmp);
        c_str = km_tmp;
    }

    minHashIterator<RepHash> it_min(c_str, len, k, Minimizer::g, RepHash(), true), it_min_end;

    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){

        //If current minimizer was not seen before
        if ((last_pos_min < it_min.getPosition()) || isForbidden){

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&c_str[min_h_res.pos]).rep(); //Get the minimizer to insert
                std::pair<hmap_min_contigs_t::iterator, bool> p = hmap_min_contigs.insert(make_pair(minz_rep, tiny_vector<size_t,tiny_vector_sz>()));
                size_t v_sz = p.first->second.size();

                pos_id_contig = (pos_id_contig & mask) | ((size_t) min_h_res.pos);

                if (!isShort){

                    mhr = min_h_res;

                    while ((v_sz >= max_abundance_lim) || ((v_sz > 0) && ((p.first->second[v_sz-1] & mask) == mask))){

                        mhr_tmp = it_min.getNewMin(mhr);
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            if ((p.first->second[v_sz-1] & mask) != mask){ // Minimizer was never signaled before as overcrowded

                                // If minimizer bin already contains abundant k-mer, just set flag for unitig overcrowding
                                if ((p.first->second[v_sz-1] & MASK_CONTIG_ID) == MASK_CONTIG_ID) p.first->second[v_sz-1] |= MASK_CONTIG_TYPE;
                                else p.first->second.push_back(mask);
                            }

                            mhr = mhr_tmp;
                            minz_rep = Minimizer(&c_str[mhr.pos]).rep();
                            p = hmap_min_contigs.insert(make_pair(minz_rep, tiny_vector<size_t,tiny_vector_sz>()));
                            v_sz = p.first->second.size();
                        }
                        else break;
                    }
                }

                tiny_vector<size_t,tiny_vector_sz>& v = p.first->second;

                if (v_sz == 0) v.push_back(pos_id_contig); //Newly created vector, just push contig ID
                else if (isShort && (v_sz >= min_abundance_lim)){ //The minimizer (is/might be) too abundant

                    isShort = false;
                    isAbundant = true;

                    it_min = it_min_end;

                    break;
                }
                else if ((v[v_sz-1] & MASK_CONTIG_ID) == MASK_CONTIG_ID){

                    if ((v_sz == 1) || (v[v_sz-2] != pos_id_contig)){

                        size_t tmp = v[v_sz-1];
                        v[v_sz-1] = pos_id_contig;
                        v.push_back(tmp);
                    }
                }
                else if (v[v_sz-1] != pos_id_contig) v.push_back(pos_id_contig);

                last_pos_min = min_h_res.pos;
                it_it_min++;
            }
        }
    }

    if (isAbundant){

        if (id_contig == v_kmers.size()) v_kmers.push_back(make_pair(km_rep, CompressedCoverage(1)));
        else v_kmers[id_contig] = make_pair(km_rep, CompressedCoverage(1));

        deleteContig(true, false, id_contig);
        if (id_contig == v_kmers.size() - 1) v_kmers.resize(v_kmers.size() - 1);

        it_min = minHashIterator<RepHash>(c_str, len, k, Minimizer::g, RepHash(), true);

        for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){

            if (last_pos_min < it_min.getPosition()){ //If current minimizer was not seen before

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&c_str[min_h_res.pos]).rep(); //Get the minimizer to insert

                    std::pair<hmap_min_contigs_t::iterator, bool> p = hmap_min_contigs.insert(make_pair(minz_rep, tiny_vector<size_t,tiny_vector_sz>()));

                    tiny_vector<size_t,tiny_vector_sz>& v = p.first->second;
                    size_t v_sz = v.size();

                    if ((v_sz > 0) && ((v[v_sz-1] & MASK_CONTIG_ID) == MASK_CONTIG_ID)) v[v_sz-1]++;
                    else v.push_back(MASK_CONTIG_ID + 1);

                    last_pos_min = min_h_res.pos;
                    it_it_min++;
                }
            }
        }

        h_kmers_ccov.insert(make_pair(km_rep, CompressedCoverage(1))).first;
    }
    else if (isShort){

        if (id_contig == v_kmers.size()) v_kmers.push_back(make_pair(km_rep, CompressedCoverage(1)));
        else v_kmers[id_contig] = make_pair(km_rep, CompressedCoverage(1));
    }
    else if (id_contig == v_contigs.size()) v_contigs.push_back(new Contig(c_str)); //Push contig to list of contigs
    else v_contigs[id_contig] = new Contig(c_str);

    return isAbundant;
}

void ContigMapper::deleteContig(const bool isShort, const bool isAbundant, const size_t id_contig){

    int k = Kmer::k;

    if (isAbundant){

        char km_str[k + 1];

        Kmer km = h_kmers_ccov.find(id_contig)->first;

        km.toString(km_str);

        minHashIterator<RepHash> it_min(km_str, k, k, Minimizer::g, RepHash(), true), it_min_end;

        for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ // Iterate over minimizers of contig to delete

            if (last_pos_min < it_min.getPosition()){ // If a new minimizer hash is found in contig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in contig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(&km_str[min_h_res.pos]).rep(); // Get canonical minimizer

                    hmap_min_contigs_t::iterator it_h = hmap_min_contigs.find(minz_rep); // Look for the minimizer in the hash table

                    if (it_h != hmap_min_contigs.end()){ // If the minimizer is found

                        tiny_vector<size_t, tiny_vector_sz>& v = it_h->second;
                        size_t last_pos_v = v.size() - 1;

                        v[last_pos_v]--;

                        if (((v[last_pos_v] & RESERVED_ID) == 0) && ((v[last_pos_v] & MASK_CONTIG_TYPE) != MASK_CONTIG_TYPE)){

                            if (last_pos_v == 0) hmap_min_contigs.erase(minz_rep);
                            else v.remove(v.size() - 1);

                            //if (last_pos_v != 0) v.remove(v.size() - 1);
                        }
                    }

                    last_pos_min = min_h_res.pos;
                    it_it_min++;
                }
            }
        }

        h_kmers_ccov.erase(km);

        //cerr << "exit deleteContig isAbundant" << endl;

        return;
    }

    bool isForbidden = false;

    size_t pos_id_contig = id_contig << 32;
    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    string str;

    if (isShort){

        str = v_kmers[id_contig].first.toString();
        pos_id_contig |= MASK_CONTIG_TYPE;
    }
    else str = v_contigs[id_contig]->seq.toString();

    const char* s = str.c_str();

    size_t len = str.size();

    minHashIterator<RepHash> it_min(s, len, k, Minimizer::g, RepHash(), true), it_min_end;

    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ // Iterate over minimizers of contig to delete

        if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer hash is found in contig to delete

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in contig to delete

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep(); // Get canonical minimizer
                hmap_min_contigs_t::iterator it_h = hmap_min_contigs.find(minz_rep); // Look for the minimizer in the hash table

                mhr = min_h_res;

                while (it_h != hmap_min_contigs.end()){ // If the minimizer is found

                    tiny_vector<size_t, tiny_vector_sz>& v = it_h->second;
                    const int v_sz = v.size();
                    bool found = false;

                    for (size_t i = 0; (i < v_sz) && !found; i++){

                        if ((v[i] & mask) == pos_id_contig){

                            v.remove(i);
                            found = true;
                        }
                    }

                    it_h = hmap_min_contigs.end();

                    if (v.size() == 0) { hmap_min_contigs.erase(minz_rep); }
                    else if (!isShort && ((v[v_sz-1] & mask) == mask)){ //Minimizer bin is overcrowded

                        mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz_rep = Minimizer(&s[mhr.pos]).rep();
                            it_h = hmap_min_contigs.find(minz_rep);
                        }
                        else break;
                    }
                }

                last_pos_min = min_h_res.pos;
                it_it_min++;
            }
        }
    }

    // The contig is deleted but its space in the contig vector is not because:
    // 1 - It would change indices in the minimizer hash table
    if (isShort) v_kmers[id_contig].first.set_deleted();
    else {

        delete v_contigs[id_contig];
        v_contigs[id_contig] = NULL;
    }

     //cerr << "exit deleteContig" << endl;
}

void ContigMapper::swapContigs(const bool isShort, const size_t id_a, const size_t id_b){

    size_t shift_id_contig_a = id_a << 32;
    size_t shift_id_contig_b = id_b << 32;

    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    string str;

    // Swap the contig pointers in v_contigs
    if (isShort){

        std::swap(v_kmers[id_a], v_kmers[id_b]);

        shift_id_contig_a |= MASK_CONTIG_TYPE;
        shift_id_contig_b |= MASK_CONTIG_TYPE;

        str = v_kmers[id_a].first.toString();
    }
    else {

        std::swap(v_contigs[id_a], v_contigs[id_b]);

        str = v_contigs[id_a]->seq.toString();
    }

    bool isForbidden = false;

    // Swap the contig IDs in the minimizer hash table
    const char* s = str.c_str();

    size_t len = str.size();

    vector<Minimizer> v_min_a;

    minHashIterator<RepHash> it_min(s, len, Kmer::k, Minimizer::g, RepHash(), true), it_min_end;

    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ // Iterate over minimizers of contig

        if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer is found in contig

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in contig

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep();

                if (!isShort){

                    hmap_min_contigs_t::iterator it_h = hmap_min_contigs.find(minz_rep);

                    if (it_h != hmap_min_contigs.end()){

                        v_min_a.push_back(minz_rep); //Add minimizer to list of minimizers

                        size_t v_sz = it_h->second.size();

                        mhr = min_h_res;

                        while ((it_h->second[v_sz-1] & mask) == mask){

                            mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                            isForbidden = true;

                            if (mhr_tmp.hash != mhr.hash){

                                minz_rep = Minimizer(&s[mhr_tmp.pos]).rep();
                                it_h = hmap_min_contigs.find(minz_rep);

                                if (it_h == hmap_min_contigs.end()) break;

                                mhr = mhr_tmp;
                                v_sz = it_h->second.size();

                                v_min_a.push_back(minz_rep);
                            }
                            else break;
                        }
                    }
                }
                else v_min_a.push_back(minz_rep); //Add minimizer to list of minimizers

                last_pos_min = min_h_res.pos;
                it_it_min++;
            }
        }
    }

    sort(v_min_a.begin(), v_min_a.end()); //Sort the list of minimizers

    for (vector<Minimizer>::iterator it_v_min = v_min_a.begin(); it_v_min != v_min_a.end(); it_v_min++){ // Iterate over minimizers

        if ((it_v_min == v_min_a.begin()) || (*it_v_min != *(it_v_min-1))){ //If the minimizer is diff. from the previous one

            hmap_min_contigs_t::iterator it_h = hmap_min_contigs.find(*it_v_min); // Look for the minimizer in the hash table

            if (it_h != hmap_min_contigs.end()){ // If the minimizer is found

                tiny_vector<size_t,tiny_vector_sz>& v_id_contigs = it_h->second;
                tiny_vector<size_t,tiny_vector_sz>::iterator it_v_c = v_id_contigs.begin();

                // Iterate over contig IDs associated with this minimizer
                for (; it_v_c != v_id_contigs.end(); it_v_c++){

                     // Swap the contig ids but do not change positions;
                    if ((*it_v_c & mask) == shift_id_contig_b) *it_v_c = shift_id_contig_a | (*it_v_c & MASK_CONTIG_POS);
                    else if ((*it_v_c & mask) == shift_id_contig_a) *it_v_c = shift_id_contig_b | (*it_v_c & MASK_CONTIG_POS);
                }
            }
        }
    }

    vector<Minimizer> v_min_b;

    str = isShort ? v_kmers[id_b].first.toString() : v_contigs[id_b]->seq.toString();
    s = str.c_str();
    len = str.size();

    isForbidden = false;

    it_min = minHashIterator<RepHash>(s, len, Kmer::k, Minimizer::g, RepHash(), true);

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ // Iterate over minimizers of contig

        if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer is found in contig

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in contig

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep();

                if (!isShort){

                    hmap_min_contigs_t::iterator it_h = hmap_min_contigs.find(minz_rep);

                    if (it_h != hmap_min_contigs.end()){

                        v_min_b.push_back(minz_rep); //Add minimizer to list of minimizers

                        size_t v_sz = it_h->second.size();

                        mhr = min_h_res;

                        while ((it_h->second[v_sz-1] & mask) == mask){

                            mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                            isForbidden = true;

                            if (mhr_tmp.hash != mhr.hash){

                                minz_rep = Minimizer(&s[mhr_tmp.pos]).rep();
                                it_h = hmap_min_contigs.find(minz_rep);

                                if (it_h == hmap_min_contigs.end()) break;

                                mhr = mhr_tmp;
                                v_sz = it_h->second.size();

                                v_min_b.push_back(minz_rep);
                            }
                            else break;
                        }
                    }
                }
                else v_min_b.push_back(minz_rep); //Add minimizer to list of minimizers

                last_pos_min = min_h_res.pos;
                it_it_min++;
            }
        }
    }

    sort(v_min_b.begin(), v_min_b.end()); //Sort the list of minimizers

    auto it_a = std::begin(v_min_a);
    auto it_a_end = std::end(v_min_a);

    //Remove Minimizers that we have already changed
    auto iter = std::remove_if (std::begin(v_min_b), std::end(v_min_b),
                                [&it_a, &it_a_end](Minimizer& minz) -> bool {
                                    while  (it_a != it_a_end && *it_a < minz) ++it_a;
                                    return (it_a != it_a_end && *it_a == minz);
                                });

    for (vector<Minimizer>::iterator it_v_min = v_min_b.begin(); it_v_min != iter; it_v_min++){ // Iterate over minimizers

        if ((it_v_min == v_min_b.begin()) || (*it_v_min != *(it_v_min-1))){ //If the minimizer is diff. from the previous one

            hmap_min_contigs_t::iterator it_h = hmap_min_contigs.find(*it_v_min); // Look for the minimizer in the hash table

            if (it_h != hmap_min_contigs.end()){ // If the minimizer is found

                tiny_vector<size_t,tiny_vector_sz>& v_id_contigs = it_h->second;
                tiny_vector<size_t,tiny_vector_sz>::iterator it_v_c = v_id_contigs.begin();

                // Iterate over contig IDs associated with this minimizer
                for (; it_v_c != v_id_contigs.end(); it_v_c++){

                     // Swap the contig ids;
                    if ((*it_v_c & mask) == shift_id_contig_a) *it_v_c = shift_id_contig_b | (*it_v_c & MASK_CONTIG_POS);
                }
            }
        }
    }
}


// use:  split, deleted = mapper.splitAllContigs()
// post: All contigs with 1 coverage somewhere have been split where the coverage is 1
//       split is the number of contigs splitted
//       deleted is the number of contigs deleted
//       Now every contig in mapper has coverage >= 2 everywhere
pair<size_t, size_t> ContigMapper::splitAllContigs() {

    size_t i;
    size_t split = 0, deleted = 0;
    size_t v_kmers_sz = v_kmers.size();
    size_t v_contigs_sz = v_contigs.size();
    size_t nxt_pos_insert = v_contigs.size();

    for (h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); it++) {

        if (!it->second.isFull()){

            deleteContig(false, true, it.getHash());
            deleted++;
        }
    }

    for (i = 0; i < v_kmers_sz;) {

        if (!v_kmers[i].second.isFull()) {

            v_kmers_sz--;

            if (i != v_kmers_sz) swapContigs(true, i, v_kmers_sz);

            deleteContig(true, false, v_kmers_sz);

            deleted++;
        }
        else i++;
    }

    for (i = 0; i < v_contigs_sz;) { // Iterate over contigs created so far

        if (!v_contigs[i]->ccov.isFull()) { //Coverage not full, contig must be splitted

            vector<pair<int,int>> sp = v_contigs[i]->ccov.splittingVector();

            if (splitContig(i, nxt_pos_insert, v_contigs_sz, v_kmers_sz, sp)) deleted++;
            else {

                split++;
                sp.clear();
            }
        }
        else i++;
    }

    if (nxt_pos_insert < v_contigs.size()) v_contigs.resize(nxt_pos_insert);
    if (v_kmers_sz < v_kmers.size()) v_kmers.resize(v_kmers_sz);

    return make_pair(split, deleted);
}

bool ContigMapper::splitContig(size_t& pos_v_contigs, size_t& nxt_pos_insert_v_contigs,
                               size_t& v_contigs_sz, size_t& v_kmers_sz, const vector<pair<int,int>>& sp){

    Contig* contig = v_contigs[pos_v_contigs];

    bool first_long_contig = true;
    bool deleted = true;

    if (!sp.empty()){

        pair<size_t, size_t> lowpair = contig->ccov.lowCoverageInfo();

        size_t lowcount = lowpair.first;
        size_t lowsum = lowpair.second;
        size_t totalcoverage = contig->coveragesum - lowsum;
        size_t ccov_size = contig->ccov.size();

        const string str = contig->seq.toString();

        int k = Kmer::k;

        for (vector<pair<int,int>>::const_iterator sit = sp.begin(); sit != sp.end(); ++sit) { //Iterate over created split contigs

            size_t pos = sit->first;
            size_t len = sit->second - pos;

            string split_str = str.substr(pos, len + k - 1); // Split contig sequence
            uint64_t cov_tmp = (totalcoverage * len) / (ccov_size - lowcount); // Split contig coverage

            if (split_str.length() == k){

                if (addContig(split_str, v_kmers_sz)) h_kmers_ccov.find(Kmer(split_str.c_str()).rep())->second.setFull();
                else {

                    v_kmers[v_kmers_sz].second.setFull(); // We don't care about the coverage per k-mer anymore
                    v_kmers_sz++;
                }
            }
            else if (first_long_contig){

                // The contig is deleted but its space in the contig vector is not because:
                // 1 - It would change indices in the minimizer hash table
                // 2 - It is going to be reused for one split contig (more efficient than deleting)
                deleteContig(false, false, pos_v_contigs);

                addContig(split_str, pos_v_contigs);

                v_contigs[pos_v_contigs]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                v_contigs[pos_v_contigs]->coveragesum = cov_tmp;

                first_long_contig = false;
            }
            else {

                addContig(split_str, nxt_pos_insert_v_contigs);

                v_contigs[nxt_pos_insert_v_contigs]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                v_contigs[nxt_pos_insert_v_contigs]->coveragesum = cov_tmp;

                nxt_pos_insert_v_contigs++;
            }
        }

        deleted = false;
    }

    if (first_long_contig){

        nxt_pos_insert_v_contigs--; //Position of the last contig in the vector which is not NULL

        if (pos_v_contigs != nxt_pos_insert_v_contigs){ // Do not proceed to swap if swap positions are the same

            swapContigs(false, pos_v_contigs, nxt_pos_insert_v_contigs); // Swap contigs

            // If the swapped contig, previously in position nxt_pos_insert, was a split contig
            // created in this method, do not try to split it again
            if (nxt_pos_insert_v_contigs >= v_contigs_sz) pos_v_contigs++;
            else v_contigs_sz--;
        }
        else v_contigs_sz--;

        deleteContig(false, false, nxt_pos_insert_v_contigs);
    }
    else pos_v_contigs++;

    return deleted;
}

// use:  joined = mapper.joinAllContigs()
// pre:  no short contigs exist in sContigs.
// post: all contigs that could be connected have been connected
//       joined is the number of joined contigs
size_t ContigMapper::joinAllContigs(vector<Kmer>* v_joins) {

    size_t i, k = Kmer::k;
    size_t joined = 0;
    size_t v_contigs_size = v_contigs.size();
    size_t v_kmers_size = v_kmers.size();

    // a and b are candidates for joining
    vector<pair<Kmer, Kmer>> joins;

    if (v_joins == NULL){

        for (h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); it++) {

            Kmer& tail = it->first;
            Kmer head_twin = tail.twin();
            Kmer fw, bw;

            const ContigMap cm(it.getHash(), 0, 0, 1, k, false, true, true);

            if (checkJoin(tail, cm, fw)) joins.push_back(make_pair(tail, fw));
            if (checkJoin(head_twin, cm, bw)) joins.push_back(make_pair(head_twin, bw));
        }

        for (i = 0; i != v_kmers_size; i++) {

            Kmer& tail = v_kmers[i].first;
            Kmer head_twin = tail.twin();
            Kmer fw, bw;

            const ContigMap cm(i, 0, 0, 1, k, true, false, true);

            if (checkJoin(tail, cm, fw)) joins.push_back(make_pair(tail, fw));
            if (checkJoin(head_twin, cm, bw)) joins.push_back(make_pair(head_twin, bw));
        }

        for (size_t i = 0; i != v_contigs_size; i++) {

            const CompressedSequence& seq = v_contigs[i]->seq;

            Kmer head_twin = seq.getKmer(0).twin();
            Kmer tail = seq.getKmer(seq.size()-k);
            Kmer fw, bw;

            const ContigMap cm(i, 0, 0, 1, seq.size(), false, false, true);

            if (checkJoin(tail, cm, fw)) joins.push_back(make_pair(tail, fw));
            if (checkJoin(head_twin, cm, bw)) joins.push_back(make_pair(head_twin, bw));
        }
    }
    else {

        Kmer fw;

        for (auto km : *v_joins){

            ContigMap cm = find(km, true);

            if (!cm.isEmpty){

                if (!cm.isShort && !cm.isAbundant){

                    if (cm.dist == 0 && cm.strand) km = km.twin();
                    else if (cm.dist != 0 && !cm.strand) km = km.twin();

                    if (checkJoin(km, cm, fw)) joins.push_back(make_pair(km, fw));
                }
                else {

                    if (checkJoin(km, cm, fw)) joins.push_back(make_pair(km, fw));
                    km = km.twin();
                    if (checkJoin(km, cm, fw)) joins.push_back(make_pair(km, fw));
                }
            }
        }

        (*v_joins).clear();
    }

    for (vector<pair<Kmer, Kmer>>::iterator it = joins.begin(); it != joins.end(); ++it) {

        const Kmer& head = it->first;
        const Kmer& tail = it->second;

        ContigMap cmHead = find(head, true);
        ContigMap cmTail = find(tail, true);

        if (!cmHead.isEmpty && !cmTail.isEmpty) {

            Kmer cmHead_head, cmTail_head;

            if (cmHead.isShort) cmHead_head = v_kmers[cmHead.pos_contig].first;
            else if (cmHead.isAbundant) cmHead_head = h_kmers_ccov.find(cmHead.pos_contig)->first;
            else cmHead_head = Kmer(v_contigs[cmHead.pos_contig]->seq.getKmer(0));

            if (cmTail.isShort) cmTail_head = v_kmers[cmTail.pos_contig].first;
            else if (cmTail.isAbundant) cmTail_head = h_kmers_ccov.find(cmTail.pos_contig)->first;
            else cmTail_head = Kmer(v_contigs[cmTail.pos_contig]->seq.getKmer(0));

            if (cmHead_head != cmTail_head) { // can't join a sequence with itself, either hairPin, loop or mobius loop

                // both kmers are still end-kmers
                bool headDir;
                bool len_k_head = cmHead.isShort || cmHead.isAbundant;

                if (len_k_head && (head == cmHead_head)) headDir = true;
                else if (!len_k_head && (head == v_contigs[cmHead.pos_contig]->seq.getKmer(v_contigs[cmHead.pos_contig]->numKmers()-1))) headDir = true;
                else if (head.twin() == cmHead_head) headDir = false;
                else continue;

                bool tailDir;
                bool len_k_tail = cmTail.isShort || cmTail.isAbundant;

                if (tail == cmTail_head) tailDir = true;
                else if (len_k_tail){
                    if (tail.twin() == cmTail_head) tailDir = false;
                    else continue;
                }
                else if (tail.twin() == v_contigs[cmTail.pos_contig]->seq.getKmer(v_contigs[cmTail.pos_contig]->numKmers()-1)) tailDir = false;
                else continue;

                //Compute join sequence
                string joinSeq, tailSeq;

                if (headDir) joinSeq = len_k_head ? cmHead_head.toString() : v_contigs[cmHead.pos_contig]->seq.toString();
                else joinSeq = len_k_head ? cmHead_head.twin().toString() : v_contigs[cmHead.pos_contig]->seq.rev().toString();

                if (tailDir) tailSeq = len_k_tail ? cmTail_head.toString() : v_contigs[cmTail.pos_contig]->seq.toString();
                else tailSeq = len_k_tail ? cmTail_head.twin().toString() : v_contigs[cmTail.pos_contig]->seq.rev().toString();

                assert(joinSeq.substr(joinSeq.size()-k+1) == tailSeq.substr(0,k-1));

                joinSeq.append(tailSeq, k-1, string::npos);

                //Compute new coverage
                uint64_t covsum;

                if (len_k_head){

                    CompressedCoverage& ccov = cmHead.isShort ? v_kmers[cmHead.pos_contig].second : h_kmers_ccov.find(cmHead.pos_contig)->second;
                    covsum = (ccov.isFull() ? ccov.cov_full : ccov.covAt(0));
                }
                else covsum = v_contigs[cmHead.pos_contig]->coveragesum;

                if (len_k_tail){

                    CompressedCoverage& ccov = cmTail.isShort ? v_kmers[cmTail.pos_contig].second : h_kmers_ccov.find(cmTail.pos_contig)->second;
                    covsum += (ccov.isFull()? ccov.cov_full : ccov.covAt(0));
                }
                else covsum += v_contigs[cmTail.pos_contig]->coveragesum;

                Contig* contig;

                if (cmHead.isShort){ //If head is a short contig, swap and delete it

                    v_kmers_size--;

                    if (cmHead.pos_contig != v_kmers_size){

                        swapContigs(true, cmHead.pos_contig, v_kmers_size);

                        // If the last contig of the vector used for the swap was the tail
                        if (cmTail.isShort && (v_kmers_size == cmTail.pos_contig)) cmTail.pos_contig = cmHead.pos_contig;
                    }

                    deleteContig(true, false, v_kmers_size);
                }
                else if (cmHead.isAbundant) deleteContig(false, true, cmHead.pos_contig);

                if (cmTail.isShort){ //If they are short

                    v_kmers_size--;

                    if (cmTail.pos_contig != v_kmers_size){

                        swapContigs(true, cmTail.pos_contig, v_kmers_size);

                        if (cmHead.isShort && (v_kmers_size == cmHead.pos_contig)) cmHead.pos_contig = cmTail.pos_contig;
                    }

                    deleteContig(true, false, v_kmers_size);
                }
                else if (cmTail.isAbundant) deleteContig(false, true, cmTail.pos_contig);

                if (len_k_head && len_k_tail){

                    addContig(joinSeq, v_contigs_size);
                    contig = v_contigs[v_contigs_size];
                    v_contigs_size++;
                }
                else if (len_k_head){

                    deleteContig(false, false, cmTail.pos_contig);
                    addContig(joinSeq, cmTail.pos_contig);
                    contig = v_contigs[cmTail.pos_contig];
                }
                else {

                    if (!len_k_tail){

                        v_contigs_size--;

                        if (cmTail.pos_contig != v_contigs_size){

                            swapContigs(false, cmTail.pos_contig, v_contigs_size);

                            if (v_contigs_size == cmHead.pos_contig) cmHead.pos_contig = cmTail.pos_contig;
                        }

                        deleteContig(false, false, v_contigs_size);
                    }

                    deleteContig(false, false, cmHead.pos_contig);
                    addContig(joinSeq, cmHead.pos_contig);
                    contig = v_contigs[cmHead.pos_contig];
                }

                contig->coveragesum = covsum;
                if (covsum >= contig->ccov.cov_full * contig->numKmers()) contig->ccov.setFull();

                joined++;
            }
        }
    }

    if (v_contigs_size < v_contigs.size()) v_contigs.resize(v_contigs_size);
    if (v_kmers_size < v_kmers.size()) v_kmers.resize(v_kmers_size);

    return joined;
}

bool ContigMapper::checkJoin(const Kmer& a, const ContigMap& cm_a, Kmer& b/*, bool& dir*/) {

    size_t k = Kmer::k, g = Minimizer::g;
    size_t fw_count = 0, bw_count = 0;

    Kmer fw_cand;

    ContigMap cm_cand, cm_cand_tmp;

    char km_tmp[k + 1];

    a.forwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
    std::pair<uint64_t*, uint64_t*> block = bf->getBlock(it_min_h);

    RepHash rep_h(k - 1), rep_h_cpy;
    rep_h.init(km_tmp);

    for (size_t i = 0; i < 4; i++) {

        km_tmp[k - 1] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendFW(alpha[i]);

        if (bf->contains(rep_h_cpy.hash(), block)){

            Kmer fw = a.forwardBase(alpha[i]);
            cm_cand_tmp = find(fw, true);

            if (!cm_cand_tmp.isEmpty) {

                fw_count++;
                if (fw_count > 1) break;
                fw_cand = fw;
                cm_cand = cm_cand_tmp;
            }
        }
    }

    if (fw_count == 1) {

        Kmer cand_head, ac_head;

        if (cm_cand.isShort) cand_head = v_kmers[cm_cand.pos_contig].first;
        else if (cm_cand.isAbundant) cand_head = h_kmers_ccov.find(cm_cand.pos_contig)->first;
        else cand_head = v_contigs[cm_cand.pos_contig]->seq.getKmer(0);

        if (cm_a.isShort) ac_head = v_kmers[cm_a.pos_contig].first;
        else if (cm_a.isAbundant) ac_head = h_kmers_ccov.find(cm_a.pos_contig)->first;
        else ac_head = v_contigs[cm_a.pos_contig]->seq.getKmer(0);

        if (cand_head != ac_head) { // not a self loop or hair-pin

            Kmer fw_cpy = fw_cand.twin();

            fw_cpy.forwardBase('A').toString(km_tmp);

            it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            block = bf->getBlock(it_min_h);

            rep_h.init(km_tmp);

            for (size_t j = 0; j < 4; j++) {

                km_tmp[k - 1] = alpha[j];

                rep_h_cpy = rep_h;
                rep_h_cpy.extendFW(alpha[j]);

                if (bf->contains(rep_h_cpy.hash(), block)){

                    Kmer fw = fw_cpy.forwardBase(alpha[j]);
                    cm_cand_tmp = find(fw, true);

                    if (!cm_cand_tmp.isEmpty) {

                        bw_count++;
                        if (bw_count > 1) break;
                    }
                }
            }

            if (bw_count == 1) { // join up

                if (cand_head == fw_cand) {

                    b = fw_cand;
                    return true;
                }

                Kmer candLast;

                if (cm_cand.isShort || cm_cand.isAbundant) candLast = cand_head;
                else candLast = v_contigs[cm_cand.pos_contig]->seq.getKmer(v_contigs[cm_cand.pos_contig]->seq.size()-k);

                if (candLast.twin() == fw_cand) {

                    b = fw_cand;
                    return true;
                }

                return true;
            }
        }
    }

    return false;
}

void ContigMapper::writeGFA(string graphfilename) {

    const size_t k = Kmer::k, g = Minimizer::g;
    const size_t v_contigs_sz = v_contigs.size();
    const size_t v_kmers_sz = v_kmers.size();

    size_t i, h_min, labelA, labelB, id = v_contigs_sz + v_kmers_sz + 1;

    char km_tmp[k + 1];

    //bool found;

    int nb_min, pos_min;

    uint64_t it_min_h;

    std::pair<uint64_t*, uint64_t*> block;

    ContigMap cand;

    minHashKmer<RepHash> mhk;

    RepHash rep_h(k - 1), rep_h_cpy;

    Contig* contig = NULL;

    ofstream graphfile;
    ostream graph(0);

    graphfile.open(graphfilename.c_str());
    graph.rdbuf(graphfile.rdbuf());
    assert(!graphfile.fail());

    KmerHashTable<size_t> idmap(h_kmers_ccov.size());

    // gfa header
    graph << "H\tVN:Z:1.0\n";

    for (labelA = 1; labelA <= v_contigs_sz; labelA++) {

        contig = v_contigs[labelA - 1];

        graph << "S\t" << labelA << "\t" << contig->seq.toString() << "\tLN:i:" <<
        contig->seq.size() << "\tXC:i:" << contig->coveragesum << "\n";
    }

    for (labelA = 1; labelA <= v_kmers_sz; labelA++) {

        const pair<Kmer, CompressedCoverage>& p = v_kmers[labelA - 1];

        graph << "S\t" << (labelA + v_contigs_sz) << "\t" << p.first.toString() << "\tLN:i:" <<
        k << "\tXC:i:" << (p.second.isFull() ? p.second.cov_full : p.second.covAt(0)) << "\n";
    }

    for (h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); it++) {

        id++;
        idmap.insert({it->first, id});

        graph << "S\t" << id << "\t" << it->first.toString() << "\tLN:i:" << k << "\tXC:i:" <<
        (it->second.isFull() ? it->second.cov_full : it->second.covAt(0)) << "\n";
    }

    // We need to deal with the tail of long contigs
    for (labelA = 1; labelA <= v_contigs_sz; labelA++) {

        Contig* contig = v_contigs[labelA - 1];

        Kmer head = contig->seq.getKmer(0);

        head.backwardBase('A').toString(km_tmp);

        mhk = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true);

        nb_min = mhk.getNbMin();
        pos_min = mhk.getPosition();
        it_min_h = mhk.getHash();

        block = bf->getBlock(it_min_h);

        //found = false;

        rep_h.init(km_tmp + 1);

        for (i = 0; i < 4; ++i) {

            km_tmp[0] = alpha[i];

            rep_h_cpy = rep_h;
            rep_h_cpy.extendBW(alpha[i]);

            if (bf->contains(rep_h_cpy.hash(), block)){

                Kmer b = head.backwardBase(alpha[i]);

                /*if (found && (nb_min == 1)) cand = find(b, pos_min, h_min, true);
                else {*/

                    cand = find(b, true);
                    /*h_min = cand.pos_min;
                    found = true;
                }*/

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_contig + 1 + (cand.isShort ? v_contigs_sz: 0);

                    graph << "L\t" << labelA << "\t-\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            }
        }

        Kmer tail = contig->seq.getKmer(contig->seq.size() - k);

        tail.forwardBase('A').toString(km_tmp);

        mhk = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true);

        nb_min = mhk.getNbMin();
        pos_min = mhk.getPosition();
        it_min_h = mhk.getHash();

        block = bf->getBlock(it_min_h);

        //found = false;

        rep_h.init(km_tmp);

        for (i = 0; i < 4; ++i) {

            km_tmp[k - 1] = alpha[i];

            rep_h_cpy = rep_h;
            rep_h_cpy.extendFW(alpha[i]);

            if (bf->contains(rep_h_cpy.hash(), block)){

                Kmer b = tail.forwardBase(alpha[i]);

                /*if (found && (nb_min == 1)) cand = find(b, pos_min, h_min, true);
                else {*/

                    cand = find(b, true);
                    /*h_min = cand.pos_min;
                    found = true;
                }*/

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_contig + 1 + (cand.isShort ? v_contigs_sz: 0);

                    graph << "L\t" << labelA << "\t+\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            }
        }
    }

    for (labelA = v_contigs_sz + 1; labelA <= v_kmers_sz + v_contigs_sz; labelA++) {

        const pair<Kmer, CompressedCoverage>& p = v_kmers[labelA - v_contigs_sz - 1];

        p.first.backwardBase('A').toString(km_tmp);

        mhk = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true);

        nb_min = mhk.getNbMin();
        pos_min = mhk.getPosition();
        it_min_h = mhk.getHash();

        block = bf->getBlock(it_min_h);

        //found = false;

        rep_h.init(km_tmp + 1);

        for (i = 0; i < 4; ++i) {

            km_tmp[0] = alpha[i];

            rep_h_cpy = rep_h;
            rep_h_cpy.extendBW(alpha[i]);

            if (bf->contains(rep_h_cpy.hash(), block)){

                Kmer b = p.first.backwardBase(alpha[i]);

                /*if (found && (nb_min == 1)) cand = find(b, pos_min, h_min, true);
                else {*/

                    cand = find(b, true);
                    /*h_min = cand.pos_min;
                    found = true;
                }*/

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_contig + 1 + (cand.isShort ? v_contigs_sz: 0);

                    graph << "L\t" << labelA << "\t-\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            }
        }

        p.first.forwardBase('A').toString(km_tmp);

        mhk = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true);

        nb_min = mhk.getNbMin();
        pos_min = mhk.getPosition();
        it_min_h = mhk.getHash();

        block = bf->getBlock(it_min_h);

        //found = false;

        rep_h.init(km_tmp);

        for (i = 0; i < 4; ++i) {

            km_tmp[k - 1] = alpha[i];

            rep_h_cpy = rep_h;
            rep_h_cpy.extendFW(alpha[i]);

            if (bf->contains(rep_h_cpy.hash(), block)){

                Kmer b = p.first.forwardBase(alpha[i]);

                /*if (found && (nb_min == 1)) cand = find(b, pos_min, h_min, true);
                else {*/

                    cand = find(b, true);
                    /*h_min = cand.pos_min;
                    found = true;
                }*/

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_contig + 1 + (cand.isShort ? v_contigs_sz: 0);

                    graph << "L\t" << labelA << "\t+\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            }
        }
    }

    for (KmerHashTable<size_t>::iterator it = idmap.begin(); it != idmap.end(); it++) {

        labelA = it->second;

        it->first.backwardBase('A').toString(km_tmp);

        mhk = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true);

        nb_min = mhk.getNbMin();
        pos_min = mhk.getPosition();
        it_min_h = mhk.getHash();

        block = bf->getBlock(it_min_h);

        //found = false;

        rep_h.init(km_tmp + 1);

        for (i = 0; i < 4; ++i) {

            km_tmp[0] = alpha[i];

            rep_h_cpy = rep_h;
            rep_h_cpy.extendBW(alpha[i]);

            if (bf->contains(rep_h_cpy.hash(), block)){

                Kmer b = it->first.backwardBase(alpha[i]);

                /*if (found && (nb_min == 1)) cand = find(b, pos_min, h_min, true);
                else {*/

                    cand = find(b, true);
                    /*h_min = cand.pos_min;
                    found = true;
                }*/

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_contig + 1 + (cand.isShort ? v_contigs_sz: 0);

                    graph << "L\t" << labelA << "\t-\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            }
        }

        it->first.forwardBase('A').toString(km_tmp);

        mhk = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true);

        nb_min = mhk.getNbMin();
        pos_min = mhk.getPosition();
        it_min_h = mhk.getHash();

        block = bf->getBlock(it_min_h);

        //found = false;

        rep_h.init(km_tmp);

        for (i = 0; i < 4; ++i) {

            km_tmp[k - 1] = alpha[i];

            rep_h_cpy = rep_h;
            rep_h_cpy.extendFW(alpha[i]);

            if (bf->contains(rep_h_cpy.hash(), block)){

                Kmer b = it->first.forwardBase(alpha[i]);

                /*if (found && (nb_min == 1)) cand = find(b, pos_min, h_min, true);
                else {*/

                    cand = find(b, true);
                    /*h_min = cand.pos_min;
                    found = true;
                }*/

                if (!cand.isEmpty) {

                    if (cand.isAbundant) labelB = idmap.find(b)->second;
                    else labelB = cand.pos_contig + 1 + (cand.isShort ? v_contigs_sz: 0);

                    graph << "L\t" << labelA << "\t+\t" << labelB << "\t" << (cand.strand ? "+" : "-") << "\t" << (k-1) << "M\n";
                }
            }
        }
    }

    graphfile.close();
}

size_t ContigMapper::removeUnitigs(bool rmIsolated, bool clipTips, vector<Kmer>& v){

    if (!rmIsolated && !clipTips) return 0;

    const bool rm_and_clip = rmIsolated && clipTips;

    const size_t k = Kmer::k, g = Minimizer::g;

    size_t v_contigs_sz = v_contigs.size();
    size_t v_kmers_sz = v_kmers.size();
    size_t i, h_min;
    size_t removed = 0;

    uint64_t it_min_h;

    int64_t j;

    char km_tmp[k + 1];

    const int lim = (clipTips ? 1 : 0);

    int nb_pred, nb_succ;

    std::pair<uint64_t*, uint64_t*> block;

    Kmer km;

    RepHash rep_h(k - 1), rep_h_cpy;

    Contig* contig = NULL;

    for (j = 0; j < v_contigs_sz; j++) {

        Contig* contig = v_contigs[j];

        if (contig->numKmers() < k){

            Kmer head = contig->seq.getKmer(0);

            head.backwardBase('A').toString(km_tmp);
            rep_h.init(km_tmp + 1);

            it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            block = bf->getBlock(it_min_h);
            nb_pred = 0;

            for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

                km_tmp[0] = alpha[i];

                rep_h_cpy = rep_h;
                rep_h_cpy.extendBW(alpha[i]);

                if (bf->contains(rep_h_cpy.hash(), block) && !find(head.backwardBase(alpha[i]), true).isEmpty){

                    nb_pred++;
                    if (clipTips) km = head.backwardBase(alpha[i]);
                }
            }

            if (nb_pred <= lim){

                Kmer tail = contig->seq.getKmer(contig->seq.size() - k);

                tail.forwardBase('A').toString(km_tmp);
                rep_h.init(km_tmp);

                it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
                block = bf->getBlock(it_min_h);
                nb_succ = 0;

                for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                    km_tmp[k - 1] = alpha[i];

                    rep_h_cpy = rep_h;
                    rep_h_cpy.extendFW(alpha[i]);

                    if (bf->contains(rep_h_cpy.hash(), block) && !find(tail.forwardBase(alpha[i]), true).isEmpty){

                        nb_succ++;
                        if (clipTips) km = tail.forwardBase(alpha[i]);
                    }
                }

                if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))) { //Unitig is isolated

                    removed++;
                    v_contigs_sz--;

                    if (j != v_contigs_sz){

                        swapContigs(false, j, v_contigs_sz),
                        j--;
                    }

                    if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
                }
            }
        }
    }

    for (j = 0; j < v_kmers_sz; j++) {

        const pair<Kmer, CompressedCoverage>& p = v_kmers[j];

        p.first.backwardBase('A').toString(km_tmp);

        rep_h.init(km_tmp + 1);

        it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
        block = bf->getBlock(it_min_h);
        nb_pred = 0;

        for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

            km_tmp[0] = alpha[i];

            rep_h_cpy = rep_h;
            rep_h_cpy.extendBW(alpha[i]);

            if (bf->contains(rep_h_cpy.hash(), block) && !find(p.first.backwardBase(alpha[i]), true).isEmpty){

                nb_pred++;
                if (clipTips) km = p.first.backwardBase(alpha[i]);
            }
        }

        if (nb_pred <= lim){

            p.first.forwardBase('A').toString(km_tmp);
            rep_h.init(km_tmp);

            it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            block = bf->getBlock(it_min_h);
            nb_succ = 0;

            for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                km_tmp[k - 1] = alpha[i];

                rep_h_cpy = rep_h;
                rep_h_cpy.extendFW(alpha[i]);

                if (bf->contains(rep_h_cpy.hash(), block) && !find(p.first.forwardBase(alpha[i]), true).isEmpty){

                    nb_succ++;
                    if (clipTips) km = p.first.forwardBase(alpha[i]);
                }
            }

            if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))) { //Unitig is isolated

                removed++;
                v_kmers_sz--;

                if (j != v_kmers_sz){

                    swapContigs(true, j, v_kmers_sz),
                    j--;
                }

                if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
            }
        }
    }

    for (h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); it++) {

        it->first.backwardBase('A').toString(km_tmp);

        rep_h.init(km_tmp + 1);

        it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
        block = bf->getBlock(it_min_h);
        nb_pred = 0;

        for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

            km_tmp[0] = alpha[i];

            rep_h_cpy = rep_h;
            rep_h_cpy.extendBW(alpha[i]);

            if (bf->contains(rep_h_cpy.hash(), block) && !find(it->first.backwardBase(alpha[i]), true).isEmpty){

                nb_pred++;
                if (clipTips) km = it->first.backwardBase(alpha[i]);
            }
        }

        if (nb_pred <= lim){

            it->first.forwardBase('A').toString(km_tmp);
            rep_h.init(km_tmp);

            it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
            block = bf->getBlock(it_min_h);
            nb_succ = 0;

            for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                km_tmp[k - 1] = alpha[i];

                rep_h_cpy = rep_h;
                rep_h_cpy.extendFW(alpha[i]);

                if (bf->contains(rep_h_cpy.hash(), block) && !find(it->first.forwardBase(alpha[i]), true).isEmpty){

                    nb_succ++;
                    if (clipTips) km = it->first.forwardBase(alpha[i]);
                }
            }

            if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))){

                removed++;

                it->second = CompressedCoverage();

                if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
            }
        }
    }

    for (j = v_contigs_sz; j < v_contigs.size(); j++) deleteContig(false, false, j);
    v_contigs.resize(v_contigs_sz);

    for (j = v_kmers_sz; j < v_kmers.size(); j++) deleteContig(true, false, j);
    v_kmers.resize(v_kmers_sz);

    for (h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); it++){

        if (it->second.size() == 0) deleteContig(false, true, it.getHash());
    }

    return removed;
}

void ContigMapper::checkIntegrity(){

    bool is_max;

    double avg_load_tiny_v = 0;

    size_t load_tiny_v[min_abundance_lim + 1] = {0};

    size_t max_sz_tiny_v = 0;
    size_t max_nb_contig_tiny_v = 0;
    size_t max_nb_kmer_tiny_v = 0;

    size_t v_sz;

    Minimizer minz_max_occ;

    for (hmap_min_contigs_t::iterator it = hmap_min_contigs.begin(); it != hmap_min_contigs.end(); it++) {

        is_max = false;

        tiny_vector<size_t,tiny_vector_sz>& v_id_contigs = it->second;

        v_sz = v_id_contigs.size();

        avg_load_tiny_v += v_sz;

        if ((v_sz > 0) && ((v_id_contigs[v_sz - 1] & MASK_CONTIG_ID) == MASK_CONTIG_ID)) load_tiny_v[min_abundance_lim]++;
        else load_tiny_v[v_sz <= min_abundance_lim ? v_sz - 1 : min_abundance_lim]++;

        Minimizer& minz = it->first;
        Minimizer minz_twin = minz.twin();

        if (v_sz > max_sz_tiny_v){

            max_sz_tiny_v = v_sz;
            minz_max_occ = minz < minz_twin ? minz : minz_twin;
            is_max = true;
            max_nb_contig_tiny_v = 0;
            max_nb_kmer_tiny_v = 0;
        }

        sort(v_id_contigs.begin(), v_id_contigs.end()); // O(N log N)

        if ((adjacent_find(v_id_contigs.begin(), v_id_contigs.end()) == v_id_contigs.end()) == false){
            cerr << "Non unique list of contig IDs for a minimizer" << endl;
            exit(EXIT_FAILURE);
        }

        for (auto id_pos: v_id_contigs){

            if (id_pos >> 32 != RESERVED_ID){

                string str;
                bool isShort = false;

                if ((id_pos & MASK_CONTIG_TYPE) == 0){

                    if (is_max) max_nb_contig_tiny_v++;

                    if ((id_pos >> 32) >= v_contigs.size()){
                        cerr << "(id = " << (id_pos >> 32) << ") >= (v_contigs.size = " << v_contigs.size() << ")" << endl;
                        cerr << "Minz = " << minz.toString() << "Minz_twin = " << minz_twin.toString() << endl;
                        exit(EXIT_FAILURE);
                    }
                    else if (v_contigs[id_pos >> 32] == NULL){
                        cerr << "v_contigs[id] == NULL" << endl;
                        exit(EXIT_FAILURE);
                    }
                    else if ((id_pos & MASK_CONTIG_POS) >= v_contigs[id_pos >> 32]->length()){
                        cerr << "(pos = " << (id_pos & MASK_CONTIG_POS) << ") >= (seq.length() = " << v_contigs[id_pos >> 32]->length() << ")" << endl;
                        exit(EXIT_FAILURE);
                    }

                    str = v_contigs[id_pos >> 32]->seq.toString().substr(id_pos & MASK_CONTIG_POS, Minimizer::g);
                }
                else {

                    isShort = true;

                    if (is_max) max_nb_kmer_tiny_v++;

                    if ((id_pos >> 32) >= v_kmers.size()){
                        cerr << "(id = " << (id_pos >> 32) << ") >= (v_kmers.size = " << v_kmers.size() << ")" << endl;
                        exit(EXIT_FAILURE);
                    }
                    else if ((id_pos & MASK_CONTIG_POS) >= Kmer::k){
                        cerr << "(pos = " << (id_pos & MASK_CONTIG_POS) << ") >= (kmer.length() = " << Kmer::k << ")" << endl;
                        exit(EXIT_FAILURE);
                    }

                    str = v_kmers[id_pos >> 32].first.toString().substr(id_pos & MASK_CONTIG_POS, Minimizer::g);
                }

                if ((str != minz.toString()) && (str != minz_twin.toString())){
                    cerr << "Couldn't find minimizer in contig at specified position, isShort = " << isShort << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    avg_load_tiny_v /= hmap_min_contigs.size();

    cerr << endl << "v_contigs.size() = " << v_contigs.size() << endl;
    cerr << "v_kmers.size() = " << v_kmers.size() << endl;
    cerr << "h_kmers_ccov.size() = " << h_kmers_ccov.size() << endl << endl;

    cerr << "Number of distinct minimizers is " << hmap_min_contigs.size() << endl << endl;

    cerr << "Average load tiny vectors is " << avg_load_tiny_v << endl;
    cerr << "Max load tiny vectors for minimizer " << minz_max_occ.toString() << " is " << max_sz_tiny_v <<
    ": " << max_nb_contig_tiny_v << " unitigs and " << max_nb_kmer_tiny_v << " kmers" << endl;

    for (size_t i = 0; i < min_abundance_lim; i++){

        cerr << "Number of tiny vectors with exactly " << (i+1) << " IDs is " << load_tiny_v[i] <<
        " (" << ((((double)load_tiny_v[i]) / ((double)hmap_min_contigs.size())) * 100) << " %)" << endl;
    }

    cerr << "Number of tiny vectors with more than " << min_abundance_lim << " ID is " << load_tiny_v[min_abundance_lim] <<
    " (" << ((((double)load_tiny_v[min_abundance_lim]) / ((double)hmap_min_contigs.size())) * 100) << " %)" << endl;
}
