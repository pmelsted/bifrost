#include "ContigMapper.hpp"
#include "CompressedSequence.hpp"
#include "KmerIterator.hpp"

#include <string>
#include <iterator>
#include <algorithm>
#include <fstream>


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

// user: i = cm.contigCount()
// pre:
// post: i is the number of contigs in the mapper
size_t ContigMapper::contigCount() const { return v_contigs.size() + hmap_kmer_contigs.size(); }

// use:  cm.mapBloomFilter(bf)
// pre:  bf != null
// post: uses the bloom filter bf to map reads
void ContigMapper::mapBloomFilter(const BlockedBloomFilter *bf) { this->bf = bf; }

// use:  cm.mapRead(km,pos,cc)
// pre:  cc is a reference to a current contig in cm, km maps to cc
// post: the coverage information in cc has been updated
void ContigMapper::mapRead(const ContigMap& cc) {

    if (cc.isEmpty) return; // nothing maps, move on

    CompressedCoverage& ccov = cc.isShort ? hmap_kmer_contigs.find(cc.pos_contig)->second : v_contigs[cc.pos_contig]->ccov;
    ccov.cover(cc.dist, cc.dist + cc.len - 1);
}

// use: b = cm.addContig(km,read)
// pre:
// post: either contig string containsin has been added and b == true
//       or it was present and the coverage information was updated, b == false
//       NOT Threadsafe!
bool ContigMapper::addContigSequence(Kmer km, const string& read, size_t pos, const string& seq) {
  // find the contig string to add

    string s;
    bool selfLoop = false;

    if (!seq.empty()) s = seq;
    else findContigSequence(km, s, selfLoop);

    size_t k = Kmer::k;

    if (selfLoop) {

        // ok, check if any other k-mer is mapped
        for (KmerIterator it(s.c_str()), it_end; it != it_end; ++it) {

            //ContigMap cm = find(it->first);
            mapRead(find(it->first));

            /*if (!cm.isEmpty) {

                CompressedCoverage& ccov = cm.isShort ? hmap_kmer_contigs.find(cm.pos_contig)->second : v_contigs[cm.pos_contig]->ccov;
                int loopSize = cm.isShort ? k : v_contigs[cm.pos_contig]->length(); // loop size in bps
                int fwMatch = stringMatch(s, read, pos) - k + 1; // length of match in k-mers
                int matchPos = (int) it->second; // position of matching k-mer within the string
                int readStart = cm.dist - matchPos;

                if (cm.strand) {

                    if (readStart < 0) {

                        readStart += loopSize;
                        int matchSize = std::min(loopSize - readStart, fwMatch);
                        ccov.cover(readStart, readStart+matchSize -1);
                        fwMatch -= matchSize;
                        readStart = 0;
                    }

                    if (fwMatch > 0) ccov.cover(readStart, readStart + fwMatch - 1);
                }
                else {

                    if (readStart < 0) {

                        readStart += loopSize;
                        int matchSize = std::min(readStart,fwMatch);
                        ccov.cover(readStart - matchSize, loopSize-1);
                        fwMatch -= matchSize;
                        readStart = loopSize -1;
                    }

                    if (fwMatch > 0) ccov.cover(readStart - fwMatch + 1, readStart);
                }

                return true;
            }*/
        }

        return true;
    }

    ContigMap cm = this->find(km);
    bool found = !cm.isEmpty;

    if (!found){

        if (s.length() == k) addShortContig(s);
        else addContig(s, v_contigs.size());
    }

    cm = findContig(km, read, pos);
    mapRead(cm);

    return found;
}

// use:  cm.findContigSequence(km, s, selfLoop)
// pre:  km is in the bloom filter
// post: s is the contig containing the kmer km
//       and the first k-mer in s is smaller (wrt. < operator)
//       than the last kmer
//       selfLoop is true of the contig is a loop or hairpin
size_t ContigMapper::findContigSequence(Kmer km, string& s, bool& selfLoop) {
    //cout << " s = " << s << endl;
    string fw_s;
    Kmer end = km;
    Kmer last = end;
    Kmer twin = km.twin();
    selfLoop = false;
    char c;
    size_t j = 0;
    size_t dummy;
    //cout << end.toString();
    while (fwBfStep(end,end,c,dummy)) {

        if (end == km) {
            selfLoop = true;
            break;
        }
        else if (end == twin) break;
        else if (end == last.twin()) break;

        j++;
        fw_s.push_back(c);
        last = end;
    }

    string bw_s;
    Kmer front = km;
    Kmer first = front;

    if (!selfLoop) {

        while (bwBfStep(front,front,c,dummy)) {

            if (front == km) {
                selfLoop = true;
                break;
            }
            else if (front == twin) break;
            else if (front == first.twin()) break;

            bw_s.push_back(c);
            first = front;
        }

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

    assert(bf != NULL);

    // need to check if we find it right away, need to treat this common case
    ContigMap cc = this->find(km);

    if (!cc.isEmpty && !cc.isShort){

        const CompressedSequence& seq = v_contigs[cc.pos_contig]->seq;
        size_t km_dist = cc.dist;
        size_t jlen = 0;
        size_t k = Kmer::k;

        if (cc.strand) jlen = seq.jump(s.c_str(), pos, cc.dist, false) - k + 1;
        else {

            jlen = seq.jump(s.c_str(), pos, cc.dist + k - 1, true) - k + 1; // match s_fw to comp(seq)_bw
            km_dist -= jlen - 1;
        }

        return ContigMap(cc.pos_contig, km_dist, jlen, cc.size, false, cc.strand);
    }

    return cc;
}

ContigMap ContigMapper::findContig(const Kmer& km, const string& s, size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h) const {

    assert(bf != NULL);

    // need to check if we find it right away, need to treat this common case
    ContigMap cc = this->find(km, pos, it_min_h);

    if (!cc.isEmpty && !cc.isShort){

        const CompressedSequence& seq = v_contigs[cc.pos_contig]->seq;
        size_t km_dist = cc.dist;
        size_t jlen = 0;
        size_t k = Kmer::k;

        if (cc.strand) jlen = seq.jump(s.c_str(), pos, cc.dist, false) - k + 1;
        else {

            jlen = seq.jump(s.c_str(), pos, cc.dist + k - 1, true) - k + 1; // match s_fw to comp(seq)_bw
            km_dist -= jlen - 1;
        }

        return ContigMap(cc.pos_contig, km_dist, jlen, cc.size, false, cc.strand);
    }

    return cc;
}

// use:  cc = cm.find(km)
// pre:
// post: cc is not empty if there is some info about km
//       in the contig map.
ContigMap ContigMapper::find(const Kmer& km) const {

    int k = Kmer::k;

    Kmer km_rep = km.rep();

    hmap_kmer_contigs_t::const_iterator it_km = hmap_kmer_contigs.find(km_rep);

    if (it_km == hmap_kmer_contigs.end()){

        size_t contig_id;
        int64_t pos_match;

        int g = Minimizer::g;

        char km_tmp[k + 1];
        km.toString(km_tmp); // Set k-mer to look-up in string version

        preAllocMinHashIterator<RepHash> it_min(km_tmp, k, k, g, RepHash(), /*false*/true);
        preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

        if (km_rep == km) km_rep = km.twin();

        while (it_it_min != it_it_min_end){

            const minHashResult& min_h_res = *it_it_min;
            Minimizer minz = Minimizer(&km_tmp[min_h_res.pos]).rep();

            it_it_min++;

            hmap_min_contigs_t::const_iterator it = hmap_min_contigs.find(minz); // Look for the minimizer in the hash table

            if (it != hmap_min_contigs.end()){ // If the minimizer is found

                for (auto contig_id_pos : it->second){ // For each contig associated with the minimizer

                    contig_id = contig_id_pos >> 32;
                    pos_match = (contig_id_pos & LOWER_32_MASK) - min_h_res.pos;

                    Contig* contig = v_contigs[contig_id];

                    if ((pos_match >= 0) && (pos_match <= contig->length() - k)){

                        Kmer km_contig = contig->seq.getKmer(pos_match);

                        if (km_contig == km) return ContigMap(contig_id, pos_match, 1, contig->length(), false, true);
                    }

                    pos_match = (contig_id_pos & LOWER_32_MASK) - k + g + min_h_res.pos;

                    if ((pos_match >= 0) && (pos_match <= contig->length() - k)){

                        Kmer km_contig = contig->seq.getKmer(pos_match);

                        if (km_contig == km_rep) return ContigMap(contig_id, pos_match, 1, contig->length(), false, false);
                    }
                }
            }
        }

        return ContigMap();
    }

    return ContigMap(it_km.getHash(), 0, 1, k, true, km == it_km->first);
}

ContigMap ContigMapper::find(const Kmer& km, bool extremities_only) const {

    int k = Kmer::k;

    Kmer km_rep = km.rep();

    hmap_kmer_contigs_t::const_iterator it_km = hmap_kmer_contigs.find(km_rep);

    if (it_km == hmap_kmer_contigs.end()){

        size_t contig_id;
        int64_t pos_match;

        int g = Minimizer::g;

        char km_tmp[k + 1];
        km.toString(km_tmp); // Set k-mer to look-up in string version

        preAllocMinHashIterator<RepHash> it_min(km_tmp, k, k, g, RepHash(), /*false*/true);
        preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

        if (km_rep == km) km_rep = km.twin();

        while (it_it_min != it_it_min_end){

            const minHashResult& min_h_res = *it_it_min;
            Minimizer minz = Minimizer(&km_tmp[min_h_res.pos]).rep();

            hmap_min_contigs_t::const_iterator it = hmap_min_contigs.find(minz); // Look for the minimizer in the hash table

            if (it != hmap_min_contigs.end()){ // If the minimizer is found

                for (auto contig_id_pos : it->second){ // For each contig associated with the minimizer

                    contig_id = contig_id_pos >> 32;
                    pos_match = (contig_id_pos & LOWER_32_MASK) - min_h_res.pos;

                    Contig* contig = v_contigs[contig_id];
                    size_t contig_len = contig->length();

                    if (!extremities_only || (pos_match == 0) || (pos_match == contig_len - k)){

                        if ((pos_match >= 0) && (pos_match <= contig_len - k)){

                            Kmer km_contig = contig->seq.getKmer(pos_match);

                            if (km_contig == km) return ContigMap(contig_id, pos_match, 1, contig_len, false, true);
                        }
                    }

                    pos_match = (contig_id_pos & LOWER_32_MASK) - k + g + min_h_res.pos;

                    if (!extremities_only || (pos_match == 0) || (pos_match == contig_len - k)){

                        if ((pos_match >= 0) && (pos_match <= contig_len - k)){

                            Kmer km_contig = contig->seq.getKmer(pos_match);

                            if (km_contig == km_rep) return ContigMap(contig_id, pos_match, 1, contig->length(), false, false);
                        }
                    }
                }
            }

            it_it_min++;
        }

        return ContigMap();
    }

    return ContigMap(it_km.getHash(), 0, 1, k, true, km == it_km->first);
}

ContigMap ContigMapper::find(const Kmer& km, const size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h) const {

    int k = Kmer::k;

    Kmer km_rep = km.rep();

    hmap_kmer_contigs_t::const_iterator it_km = hmap_kmer_contigs.find(km_rep);

    if (it_km == hmap_kmer_contigs.end()){

        size_t contig_id;
        int64_t pos_match;

        int g = Minimizer::g;

        preAllocMinHashResultIterator<RepHash> it_it_min = *it_min_h, it_it_min_end;

        if (km_rep == km) km_rep = km.twin();

        while (it_it_min != it_it_min_end){

            const minHashResult& min_h_res = *it_it_min;
            Minimizer minz = Minimizer(&it_min_h.s[min_h_res.pos]).rep();

            hmap_min_contigs_t::const_iterator it = hmap_min_contigs.find(minz); // Look for the minimizer in the hash table

            if (it != hmap_min_contigs.end()){ // If the minimizer is found

                for (auto contig_id_pos : it->second){ // For each contig associated with the minimizer

                    contig_id = contig_id_pos >> 32;
                    pos_match = (contig_id_pos & LOWER_32_MASK) - (min_h_res.pos - pos);

                    Contig* contig = v_contigs[contig_id];

                    if ((pos_match >= 0) && (pos_match <= contig->length() - k)){

                        Kmer km_contig = contig->seq.getKmer(pos_match);

                        if (km_contig == km) return ContigMap(contig_id, pos_match, 1, contig->length(), false, true);
                    }

                    pos_match = (contig_id_pos & LOWER_32_MASK) - k + g + (min_h_res.pos - pos);

                    if ((pos_match >= 0) && (pos_match <= contig->length() - k)){

                        Kmer km_contig = contig->seq.getKmer(pos_match);

                        if (km_contig == km_rep) return ContigMap(contig_id, pos_match, 1, contig->length(), false, false);
                    }
                }
            }

            it_it_min++;
        }

        return ContigMap();
    }

    return ContigMap(it_km.getHash(), 0, 1, k, true, km == it_km->first);
}

/*size_t ContigMapper::find(const size_t pos_kmer, const preAllocMinHashIterator<RepHash>& it_min_h) const {

    size_t last_pos = pos_kmer;

    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min_h, it_it_min_end;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;

        hmap_min_contigs_t::const_iterator it = hmap_min_contigs.find(Minimizer(&it_min_h.s[min_h_res.pos]).rep());

        if (it != hmap_min_contigs.end()) return last_pos - pos_kmer;// If the minimizer is found

        last_pos = min_h_res.pos;

        it_it_min++;
    }

    return last_pos - pos_kmer;
}*/

// use:  b = cm.bwBfStep(km,front,c,deg)
// pre:  km is in the bloom filter
// post: b is true if km is inside a contig, in that
//       case end is the bw link and c is the nucleotide used for the link.
//       if b is false, front and c are not updated
//       if km is an isolated self link (e.g. 'AAA') i.e. end == km then returns false
//       deg is the backwards degree of the front
bool ContigMapper::bwBfStep(Kmer km, Kmer& front, char& c, size_t& deg) const {

    size_t i,j = -1, k = Kmer::k;
    size_t nb_neigh = 0;

    const bool neighbor_hash = true;

    char km_tmp[k + 1];

    front.backwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, Minimizer::g, RepHash(), neighbor_hash).getHash();
    uint64_t* block = bf->getBlock(it_min_h);

    RepHash rep_h(k - 1), rep_h_cpy;
    rep_h.init(km_tmp + 1);

    for (i = 0; i < 4; ++i) {

        km_tmp[0] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendBW(alpha[i]);

        if (bf->contains(rep_h_cpy.hash(), block)){

            j = i;
            ++nb_neigh;

            if (nb_neigh > 1) break;
        }
    }

    if (nb_neigh != 1) {

        deg = nb_neigh;
        return false;
    }

    // only one k-mer in the bw link
    deg = 1;
    nb_neigh = 0;

    Kmer bw = front.backwardBase(alpha[j]);

    bw.forwardBase('A').toString(km_tmp);

    it_min_h = minHashKmer<RepHash>(km_tmp, k, Minimizer::g, RepHash(), neighbor_hash).getHash();
    block = bf->getBlock(it_min_h);

    rep_h.init(km_tmp);

    for (i = 0; i < 4; ++i) {

        km_tmp[k - 1] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendFW(alpha[i]);

        if (bf->contains(rep_h_cpy.hash(), block)){

            ++nb_neigh;
            if (nb_neigh > 1) break;
        }
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
}

// use:  b = cm.fwBfStep(km,end,c,deg)
// pre:  km is in the bloom filter
// post: b is true if km is inside a contig, in that
//       case end is the fw link and c is the nucleotide used for the link.
//       if b is false, end and c are not updated
//       if km is an isolated self link (e.g. 'AAA') i.e. end == km then returns false
//       deg is the degree of the end
bool ContigMapper::fwBfStep(Kmer km, Kmer& end, char& c, size_t& deg) const {

    size_t i,j = -1, k = Kmer::k, g = Minimizer::g;
    size_t nb_neigh = 0;

    const bool neighbor_hash = true;

    char km_tmp[k + 1];

    end.forwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), neighbor_hash).getHash();
    uint64_t* block = bf->getBlock(it_min_h);

    RepHash rep_h(k - 1), rep_h_cpy;
    rep_h.init(km_tmp);

    for (i = 0; i < 4; ++i) {

        km_tmp[k - 1] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendFW(alpha[i]);

        if (bf->contains(rep_h_cpy.hash(), block)){

            j = i;
            ++nb_neigh;

            if (nb_neigh > 1) break;
        }
    }

    if (nb_neigh != 1) {

        deg = nb_neigh;
        return false;
    }
    // only one k-mer in fw link
    deg = 1;
    nb_neigh = 0;

    Kmer fw = end.forwardBase(alpha[j]);

    // check bw from fw link
    fw.backwardBase('A').toString(km_tmp);

    it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), neighbor_hash).getHash();
    block = bf->getBlock(it_min_h);

    rep_h.init(km_tmp + 1);

    for (i = 0; i < 4; ++i) {

        km_tmp[0] = alpha[i];

        rep_h_cpy = rep_h;
        rep_h_cpy.extendBW(alpha[i]);

        if (bf->contains(rep_h_cpy.hash(), block)){

            ++nb_neigh;

            if (nb_neigh > 1) break;
        }
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
}

// NOT THREAD SAFE
void ContigMapper::addContig(const string& str_contig, const size_t id_contig){

    assert(id_contig <= LOWER_32_MASK);
    assert(id_contig <= v_contigs.size());

    unsigned int k = Kmer::k;

    const char* c = str_contig.c_str();
    size_t len = str_contig.size();
    size_t pos_id_contig = id_contig << 32;

    if (id_contig == v_contigs.size()) v_contigs.push_back(new Contig(c)); //Push contig to list of contigs
    else v_contigs[id_contig] = new Contig(c);

    minHashIterator<RepHash> it_min(c, len, k, Minimizer::g, RepHash(), /*false*/true), it_min_end;

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ //Link minimizers of contig to contig

        if (last_pos_min < it_min.getPosition()){

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

            while (it_it_min != it_it_min_end){

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&c[min_h_res.pos]).rep();
                pos_id_contig = (pos_id_contig & UPPER_32_MASK) | ((size_t) min_h_res.pos);

                it_it_min++;

                std::pair<hmap_min_contigs_t::iterator, bool> p_insert = hmap_min_contigs.insert(make_pair(minz_rep, tiny_vector<size_t,tiny_vector_sz>()));

                hmap_min_contigs_t::iterator it_h_contig_ids = p_insert.first;
                tiny_vector<size_t,tiny_vector_sz>& v = it_h_contig_ids->second;

                if (p_insert.second) v.push_back(pos_id_contig);
                else {

                    tiny_vector<size_t,tiny_vector_sz>::iterator it_v = std::find(v.begin(), v.end(), pos_id_contig);
                    if (it_v == v.end()) v.push_back(pos_id_contig);
                }

                last_pos_min = min_h_res.pos;
            }
        }
    }
}

KmerHashTable<CompressedCoverage>::iterator ContigMapper::addShortContig(const string& str_contig){

    assert(str_contig.size() == Kmer::k);

    return hmap_kmer_contigs.insert(make_pair(Kmer(str_contig.c_str()).rep(), CompressedCoverage(1))).first;
}

KmerHashTable<CompressedCoverage>::iterator ContigMapper::deleteShortContig(const size_t id_contig){
    return hmap_kmer_contigs.erase(id_contig);
}

void ContigMapper::deleteContig(const size_t id_contig){

    assert(id_contig < v_contigs.size());

    const string& str = v_contigs[id_contig]->seq.toString();
    const char* s = str.c_str();
    size_t len = str.size();
    size_t pos_id_contig = id_contig << 32;

    minHashIterator<RepHash> it_min(s, len, Kmer::k, Minimizer::g, RepHash(), /*false*/true), it_min_end;

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ // Iterate over minimizers of contig to delete

        if (last_pos_min < it_min.getPosition()){ // If a new minimizer hash is found in contig to delete

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in contig to delete

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(&s[min_h_res.pos]).rep(); // Get canonical minimizer

                it_it_min++;

                hmap_min_contigs_t::iterator it_h = hmap_min_contigs.find(minz_rep); // Look for the minimizer in the hash table

                if (it_h != hmap_min_contigs.end()){ // If the minimizer is found

                    tiny_vector<size_t, tiny_vector_sz> v_id_contigs_tmp;
                    tiny_vector<size_t, tiny_vector_sz>& v_id_contigs = it_h->second;
                    tiny_vector<size_t, tiny_vector_sz>::iterator it_v_c = v_id_contigs.begin();

                    for (; it_v_c != v_id_contigs.end(); it_v_c++){

                        if ((*it_v_c & UPPER_32_MASK) != pos_id_contig) v_id_contigs_tmp.push_back(*it_v_c);
                    }

                    if (v_id_contigs_tmp.size() == 0) hmap_min_contigs.erase(minz_rep);
                    else v_id_contigs = v_id_contigs_tmp;
                }

                last_pos_min = min_h_res.pos;
            }
        }
    }

    // The contig is deleted but its space in the contig vector is not because:
    // 1 - It would change indices in the minimizer hash table
    delete v_contigs[id_contig];
    v_contigs[id_contig] = NULL;
}

void ContigMapper::swapContigs(const size_t id_contig_a, const size_t id_contig_b){

    assert(id_contig_a < v_contigs.size());
    assert(id_contig_b < v_contigs.size());

    // Swap the contig pointers in v_contigs
    std::swap(v_contigs[id_contig_a], v_contigs[id_contig_b]);

    // Swap the contig IDs in the minimizer hash table
    string str = v_contigs[id_contig_a]->seq.toString();
    const char* s = str.c_str();
    size_t len = str.size();

    size_t shift_id_contig_a = id_contig_a << 32;
    size_t shift_id_contig_b = id_contig_b << 32;

    vector<Minimizer> v_min_a;

    minHashIterator<RepHash> it_min(s, len, Kmer::k, Minimizer::g, RepHash(), /*false*/true), it_min_end;

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ // Iterate over minimizers of contig

        if (last_pos_min < it_min.getPosition()){ // If a new minimizer is found in contig

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in contig

                const minHashResult& min_h_res = *it_it_min;

                v_min_a.push_back(Minimizer(&s[min_h_res.pos]).rep()); //Add minimizer to list of minimizers

                it_it_min++;
                last_pos_min = min_h_res.pos;
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
                    if ((*it_v_c & UPPER_32_MASK) == shift_id_contig_b) *it_v_c = shift_id_contig_a | (*it_v_c & LOWER_32_MASK);
                    else if ((*it_v_c & UPPER_32_MASK) == shift_id_contig_a) *it_v_c = shift_id_contig_b | (*it_v_c & LOWER_32_MASK);
                }
            }
        }
    }

    vector<Minimizer> v_min_b;

    str = v_contigs[id_contig_b]->seq.toString();
    s = str.c_str();
    len = str.size();

    it_min = minHashIterator<RepHash>(s, len, Kmer::k, Minimizer::g, RepHash(), /*false*/true);

    for (int64_t last_pos_min = -1; it_min != it_min_end; it_min++){ // Iterate over minimizers of contig

        if (last_pos_min < it_min.getPosition()){ // If a new minimizer is found in contig

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in contig

                const minHashResult& min_h_res = *it_it_min;

                v_min_b.push_back(Minimizer(&s[min_h_res.pos]).rep()); //Add minimizer to list of minimizers

                it_it_min++;
                last_pos_min = min_h_res.pos;
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
                    if ((*it_v_c & UPPER_32_MASK) == shift_id_contig_a) *it_v_c = shift_id_contig_b | (*it_v_c & LOWER_32_MASK);
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

    Kmer km;

    size_t k = Kmer::k;
    size_t split = 0, deleted = 0;

    for (hmap_kmer_contigs_t::iterator it = hmap_kmer_contigs.begin(); it != hmap_kmer_contigs.end();) {

        if (!it->second.isFull()){

            it = hmap_kmer_contigs.erase(it);
            deleted++;
        }
        else it++;
    }

    size_t v_contigs_size = v_contigs.size();
    size_t nxt_pos_insert = v_contigs.size();

    typedef vector<pair<int,int>> split_vector_t;

    for (size_t i = 0; i < v_contigs_size;) { // Iterate over contigs created so far

        if (!v_contigs[i]->ccov.isFull()) { //Coverage not full, contig must be splitted

            Contig* contig = v_contigs[i];

            pair<size_t, size_t> lowpair = contig->ccov.lowCoverageInfo();
            size_t lowcount = lowpair.first;
            size_t lowsum = lowpair.second;
            size_t totalcoverage = contig->coveragesum - lowsum;
            size_t ccov_size = contig->ccov.size();

            // remember pieces
            split_vector_t sp = contig->ccov.splittingVector();

            bool first_long_contig = true;

            if (!sp.empty()){

                const string& str = contig->seq.toString();

                for (split_vector_t::iterator sit = sp.begin(); sit != sp.end(); ++sit) { //Iterate over created split contigs

                    size_t pos = sit->first;
                    size_t len = sit->second - pos;

                    string split_str = str.substr(pos, len + k - 1); // Split contig sequence
                    uint64_t cov_tmp = (totalcoverage * len) / (ccov_size - lowcount); // Split contig coverage

                    if (split_str.length() == k){

                        KmerHashTable<CompressedCoverage>::iterator it = addShortContig(split_str);
                        it->second.setFull(); // We don't care about the coverage per k-mer anymore
                    }
                    else if (first_long_contig){

                        // The contig is deleted but its space in the contig vector is not because:
                        // 1 - It would change indices in the minimizer hash table
                        // 2 - It is going to be reused for one split contig (more efficient than deleting)
                        deleteContig(i);

                        addContig(split_str, i);

                        v_contigs[i]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                        v_contigs[i]->coveragesum = cov_tmp;

                        first_long_contig = false;
                    }
                    else {

                        addContig(split_str, nxt_pos_insert);

                        v_contigs[nxt_pos_insert]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                        v_contigs[nxt_pos_insert]->coveragesum = cov_tmp;

                        nxt_pos_insert++;
                    }
                }

                sp.clear();

                split++;
            }
            else deleted++;

            if (first_long_contig){

                nxt_pos_insert--; //Position of the last contig in the vector which is not NULL

                if (i != nxt_pos_insert){ // Do not proceed to swap if swap positions are the same

                    swapContigs(i, nxt_pos_insert); // Swap contigs

                    // If the swapped contig, previously in position nxt_pos_insert, was a split contig
                    // created in this method, do not try to split it again
                    if (nxt_pos_insert >= v_contigs_size) i++;

                }

                deleteContig(nxt_pos_insert); // Contig swapped is deleted

                v_contigs_size--;
            }
            else i++;
        }
        else i++;
    }

    if (nxt_pos_insert < v_contigs.size()) v_contigs.resize(nxt_pos_insert);

    return make_pair(split,deleted);
}

// use:  joined = mapper.joinAllContigs()
// pre:  no short contigs exist in sContigs.
// post: all contigs that could be connected have been connected
//       joined is the number of joined contigs
size_t ContigMapper::joinAllContigs() {

    size_t joined = 0;
    size_t k = Kmer::k;
    size_t v_contigs_size = v_contigs.size();

    // a and b are candidates for joining
    typedef pair<Kmer, Kmer> Join_t;
    vector<Join_t> joins;

    for (size_t i = 0; i != v_contigs.size(); i++) {

        const CompressedSequence& seq = v_contigs[i]->seq;

        Kmer head_twin = seq.getKmer(0).twin();
        Kmer tail = seq.getKmer(seq.size()-k);

        Kmer fw, bw;
        bool fw_dir = true, bw_dir = true;

        //if (checkJoin(tail, fw, fw_dir)) joins.push_back(make_pair(tail, fw));
        //if (checkJoin(head_twin, bw, bw_dir)) joins.push_back(make_pair(head_twin, bw));

        const ContigMap cm(i, 0, 1, seq.size(), false, true);

        if (checkJoin(tail, cm, fw, fw_dir)) joins.push_back(make_pair(tail, fw));
        if (checkJoin(head_twin, cm, bw, bw_dir)) joins.push_back(make_pair(head_twin, bw));
    }

    for (hmap_kmer_contigs_t::iterator it = hmap_kmer_contigs.begin(); it != hmap_kmer_contigs.end(); it++) {

        Kmer head_twin = it->first.twin();
        Kmer fw, bw;

        bool fw_dir = true, bw_dir = true;

        //if (checkJoin(it->first, fw, fw_dir)) joins.push_back(make_pair(it->first, fw));
        //if (checkJoin(head_twin, bw, bw_dir)) joins.push_back(make_pair(head_twin, bw));

        const ContigMap cm(it.getHash(), 0, 1, k, true, true);

        if (checkJoin(it->first, cm, fw, fw_dir)) joins.push_back(make_pair(it->first, fw));
        if (checkJoin(head_twin, cm, bw, bw_dir)) joins.push_back(make_pair(head_twin, bw));
    }

    for (vector<Join_t>::iterator it = joins.begin(); it != joins.end(); ++it) {

        const Kmer& head = it->first;
        const Kmer& tail = it->second;

        const ContigMap cmHead = find(head, true);
        ContigMap cmTail = find(tail, true);

        if (!cmHead.isEmpty && !cmTail.isEmpty) {

            const Kmer cmHead_head = cmHead.isShort ? hmap_kmer_contigs.find(cmHead.pos_contig)->first : Kmer(v_contigs[cmHead.pos_contig]->seq.getKmer(0));
            const Kmer cmTail_head = cmTail.isShort ? hmap_kmer_contigs.find(cmTail.pos_contig)->first : Kmer(v_contigs[cmTail.pos_contig]->seq.getKmer(0));

            if (cmHead_head != cmTail_head) { // can't join a sequence with itself, either hairPin, loop or mobius loop

                // both kmers are still end-kmers
                bool headDir = true;

                if (cmHead.isShort && head == cmHead_head) headDir = true;
                else if (!cmHead.isShort && (head == v_contigs[cmHead.pos_contig]->seq.getKmer(v_contigs[cmHead.pos_contig]->numKmers()-1))) headDir = true;
                else if (cmHead.isShort || (head.twin() == cmHead_head)) headDir = false;
                else continue; // can't join up

                bool tailDir = true;

                if (tail == cmTail_head) tailDir = true;
                else if (cmTail.isShort || (tail.twin() == v_contigs[cmTail.pos_contig]->seq.getKmer(v_contigs[cmTail.pos_contig]->numKmers()-1))) tailDir = false;
                else continue; // can't join up

                //Compute join sequence
                string joinSeq, tailSeq;

                if (headDir) joinSeq = cmHead.isShort ? cmHead_head.toString() : v_contigs[cmHead.pos_contig]->seq.toString();
                else joinSeq = cmHead.isShort ? cmHead_head.twin().toString() : v_contigs[cmHead.pos_contig]->seq.rev().toString();

                if (tailDir) tailSeq = cmTail.isShort ? cmTail_head.toString() : v_contigs[cmTail.pos_contig]->seq.toString();
                else tailSeq = cmTail.isShort ? cmTail_head.twin().toString() : v_contigs[cmTail.pos_contig]->seq.rev().toString();

                assert(joinSeq.substr(joinSeq.size()-k+1) == tailSeq.substr(0,k-1));

                joinSeq.append(tailSeq, k-1, string::npos);

                //Compute new coverage
                uint64_t covsum;

                if (cmHead.isShort){

                    CompressedCoverage& ccov = hmap_kmer_contigs.find(cmHead.pos_contig)->second;
                    covsum = (ccov.isFull() ? ccov.cov_full : ccov.covAt(0));
                }
                else covsum = v_contigs[cmHead.pos_contig]->coveragesum;

                if (cmTail.isShort){

                    CompressedCoverage& ccov = hmap_kmer_contigs.find(cmTail.pos_contig)->second;
                    covsum += (ccov.isFull()? ccov.cov_full : ccov.covAt(0));
                }
                else covsum += v_contigs[cmTail.pos_contig]->coveragesum;

                Contig* contig;

                if (cmHead.isShort && cmTail.isShort){

                    deleteShortContig(cmHead.pos_contig);
                    deleteShortContig(cmTail.pos_contig);

                    addContig(joinSeq, v_contigs_size);
                    contig = v_contigs[v_contigs_size];

                    v_contigs_size++;
                }
                else if (!cmHead.isShort && !cmTail.isShort){

                    v_contigs_size--;

                    if (cmHead.pos_contig != v_contigs_size){

                        swapContigs(cmHead.pos_contig, v_contigs_size);

                        if (v_contigs_size == cmTail.pos_contig) cmTail.pos_contig = cmHead.pos_contig;
                    }

                    deleteContig(v_contigs_size);
                    deleteContig(cmTail.pos_contig);

                    addContig(joinSeq, cmTail.pos_contig);
                    contig = v_contigs[cmTail.pos_contig];
                }
                else if (cmHead.isShort){

                    deleteShortContig(cmHead.pos_contig);
                    deleteContig(cmTail.pos_contig);

                    addContig(joinSeq, cmTail.pos_contig);
                    contig = v_contigs[cmTail.pos_contig];
                }
                else{

                    deleteShortContig(cmTail.pos_contig);
                    deleteContig(cmHead.pos_contig);

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

    return joined;
}

// use:  r = mapper.checkJoin(a,b,dir)
// pre:  a is and endpoint
// pos:  r is true iff a->b (dir is true) or a->~b (dir is false)
//       and this is the only such pair with a or b in it
/*bool ContigMapper::checkJoin(Kmer a, Kmer& b, bool& dir) {

    size_t k = Kmer::k;
    size_t fw_count = 0, bw_count = 0;

    bool fw_dir, bw_dir;

    Kmer fw_cand, bw_cand;

    for (size_t i = 0; i < 4; i++) {

        Kmer fw = a.forwardBase(alpha[i]);

        if (checkEndKmer(fw, fw_dir)) {

            fw_count++;
            fw_cand = fw;
        }
    }

    if (fw_count == 1) {

        ContigMap cand = find(fw_cand);
        ContigMap ac = find(a);

        Kmer cand_head = cand.isShort ? hmap_kmer_contigs.find(cand.pos_contig)->first : v_contigs[cand.pos_contig]->seq.getKmer(0);
        Kmer ac_head = ac.isShort ? hmap_kmer_contigs.find(ac.pos_contig)->first : v_contigs[ac.pos_contig]->seq.getKmer(0);

        if (cand_head != ac_head) { // not a self loop or hair-pin

            // no self-loop
            for (size_t j = 0; j < 4; j++) {

                Kmer bw = fw_cand.backwardBase(alpha[j]);

                if (checkEndKmer(bw.twin(), bw_dir)) {

                    bw_count++;
                    bw_cand = bw;
                }
            }

            if (bw_count == 1) {
                // join up
                Kmer candLast;

                if (cand.isShort) candLast = cand_head;
                else {
                    CompressedSequence& candSeq = v_contigs[cand.pos_contig]->seq;
                    candLast = candSeq.getKmer(candSeq.size()-k);
                }

                if (cand_head == fw_cand) {

                    b = fw_cand;
                    dir = true;

                    return true;
                }

                if (candLast.twin() == fw_cand) {

                    b = fw_cand;
                    dir = false;

                    return true;
                }

                return true;
            }
        }
        else return false;
    }

    return false;
}*/

bool ContigMapper::checkJoin(const Kmer& a, const ContigMap& cm, Kmer& b, bool& dir) {

    size_t k = Kmer::k, g = Minimizer::g;
    size_t fw_count = 0, bw_count = 0;

    Kmer fw_cand;

    ContigMap cm_cand, cm_cand_tmp;

    char km_tmp[k + 1];

    a.forwardBase('A').toString(km_tmp);

    uint64_t it_min_h = minHashKmer<RepHash>(km_tmp, k, g, RepHash(), true).getHash();
    uint64_t* block = bf->getBlock(it_min_h);

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

        Kmer cand_head = cm_cand.isShort ? hmap_kmer_contigs.find(cm_cand.pos_contig)->first : v_contigs[cm_cand.pos_contig]->seq.getKmer(0);
        Kmer ac_head = cm.isShort ? hmap_kmer_contigs.find(cm.pos_contig)->first : v_contigs[cm.pos_contig]->seq.getKmer(0);

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
                    dir = true;

                    return true;
                }

                Kmer candLast = cm_cand.isShort? cand_head : v_contigs[cm_cand.pos_contig]->seq.getKmer(v_contigs[cm_cand.pos_contig]->seq.size()-k);

                if (candLast.twin() == fw_cand) {

                    b = fw_cand;
                    dir = false;

                    return true;
                }

                return true;
            }
        }
    }

    return false;
}

// use: r = mapper.checkEndKmer(b,dir)
// pre:
// post: true iff b is an end contig in mapper and r is
//       set to true if beginning or false if b is the end
// TODO: test checkEndKmer, w.r.t. circular mapping
//       probably some problem there
/*bool ContigMapper::checkEndKmer(Kmer b, bool& dir) {

    ContigMap cand = find(b, true);

    if (cand.isEmpty) return false;

    size_t seqSize = cand.isShort ? 1 : v_contigs[cand.pos_contig]->numKmers();

    if (cand.dist == 0) {

        dir = true;
        return true;
    }
    else if (cand.dist == seqSize-1) {

        dir = false;
        return true;
    }

    return false;
}*/

void ContigMapper::writeGFA(int count1, string graphfilename) {

    size_t id = 1, id_limit_short_contig;
    size_t k = Kmer::k;
    size_t v_contigs_sz = v_contigs.size();

    Contig* contig = NULL;

    ofstream graphfile;
    ostream graph(0);

    graphfile.open(graphfilename.c_str());
    graph.rdbuf(graphfile.rdbuf());
    assert(!graphfile.fail());

    // gfa header
    graph << "H\tVN:Z:1.0\n";

    KmerHashTable<size_t> idmap(hmap_kmer_contigs.pop);

    for (; id <= v_contigs_sz; id++) {

        contig = v_contigs[id - 1];

        graph << "S\t" << id << "\t" << contig->seq.toString() << "\tLN:i:" <<
        contig->seq.size() << "\tXC:i:" << contig->coveragesum << "\n";
    }

    for (hmap_kmer_contigs_t::iterator it = hmap_kmer_contigs.begin(); it != hmap_kmer_contigs.end(); it++) {

        id++;
        idmap.insert({it->first, id});

        graph << "S\t" << id << "\t" << it->first.toString() << "\tLN:i:" << k << "\tXC:i:" <<
        (it->second.isFull() ? it->second.cov_full : it->second.covAt(0)) << "\n";
    }

    for (KmerHashTable<size_t>::iterator it = idmap.begin(); it != idmap.end(); it++) {

        size_t labelA = it->second;
        size_t labelB = 0;

        for (auto a : alpha) {

            Kmer b = it->first.backwardBase(a);
            ContigMap cand = find(b);

            if (!cand.isEmpty) {

                if (cand.isShort) labelB = idmap.find(b)->second;
                else labelB = cand.pos_contig + 1;

                if (!cand.strand && (labelA < labelB))
                    graph << "L\t" << labelA << "\t-\t" << labelB << "\t+\t" << (k-1) << "M\n";
            }
        }

        for (auto a : alpha) {

            Kmer b = it->first.forwardBase(a);
            ContigMap cand = find(b);

            if (!cand.isEmpty) {

                if (cand.isShort) labelB = idmap.find(b)->second;
                else labelB = cand.pos_contig + 1;

                if (cand.strand) graph << "L\t" << labelA << "\t+\t" << labelB << "\t+\t" << (k-1) << "M\n";
                else if (labelA <= labelB) graph << "L\t" << labelA << "\t+\t" << labelB << "\t-\t" << (k-1) << "M\n";
            }
        }
    }

    // We need to deal with the tail of long contigs
    for (size_t labelA = 1; labelA <= v_contigs_sz; labelA++) {

        Contig* contig = v_contigs[labelA - 1];

        Kmer head = contig->seq.getKmer(0);
        Kmer tail = contig->seq.getKmer(contig->seq.size() - k);

        size_t labelB = 0;

        for (auto a : alpha) {

            Kmer b = head.backwardBase(a);
            ContigMap cand = find(b);

            if (!cand.isEmpty) {

                if (cand.isShort) labelB = idmap.find(b)->second;
                else labelB = cand.pos_contig + 1;

                if (!cand.strand && (labelA < labelB))
                    graph << "L\t" << labelA << "\t-\t" << labelB << "\t+\t" << (k-1) << "M\n";
            }
        }

        for (auto a : alpha) {

            Kmer b = tail.forwardBase(a);
            ContigMap cand = find(b);

            if (!cand.isEmpty) {

                if (cand.isShort) labelB = idmap.find(b)->second;
                else labelB = cand.pos_contig + 1;

                if (cand.strand) graph << "L\t" << labelA << "\t+\t" << labelB << "\t+\t" << (k-1) << "M\n";
                else if (labelA <= labelB) graph << "L\t" << labelA << "\t+\t" << labelB << "\t-\t" << (k-1) << "M\n";
            }
        }
    }

    graphfile.close();
}

void ContigMapper::checkIntegrity(){

    for (hmap_min_contigs_t::iterator it = hmap_min_contigs.begin(); it != hmap_min_contigs.end(); it++) {

        tiny_vector<size_t,tiny_vector_sz>& v_id_contigs = it->second;

        sort(v_id_contigs.begin(), v_id_contigs.end()); // O(N log N)

        if ((adjacent_find(v_id_contigs.begin(), v_id_contigs.end()) == v_id_contigs.end()) == false){
            cerr << "Non unique list of contig IDs for a minimizer" << endl;
            exit(EXIT_FAILURE);
        }

        for (auto id_pos: v_id_contigs){

            if ((id_pos >> 32) >= v_contigs.size()){
                cerr << "(id = " << (id_pos >> 32) << ") >= (v_contigs.size = " << v_contigs.size() << ")" << endl;
                exit(EXIT_FAILURE);
            }

            if (v_contigs[id_pos >> 32] == NULL){
                cerr << "v_contigs[id] == NULL" << endl;
                exit(EXIT_FAILURE);
            }

            if ((id_pos & LOWER_32_MASK) >= v_contigs[id_pos >> 32]->length()){
                cerr << "(pos = " << (id_pos & LOWER_32_MASK) << ") >= (seq.length() = " << v_contigs[id_pos >> 32]->length() << ")" << endl;
                exit(EXIT_FAILURE);
            }

            Minimizer& minz = it->first;
            Minimizer minz_twin = minz.twin();

            string str = v_contigs[id_pos >> 32]->seq.toString().substr(id_pos & LOWER_32_MASK, Minimizer::g);

            if ((str != minz.toString()) && (str != minz_twin.toString())){
                cerr << "Couldn't find minimizer in contig at specified position" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    size_t nb_contig_length_k = 0;

    for (auto contig: v_contigs){

        if (contig == NULL){
            cerr << "contig == NULL" << endl;
            exit(EXIT_FAILURE);
        }

        nb_contig_length_k += (contig->length() == Kmer::k);
    }

    if (nb_contig_length_k != 0){
        cerr << nb_contig_length_k << " contigs of length k are in the vector of long contigs" << endl;
        exit(EXIT_FAILURE);
    }
}
