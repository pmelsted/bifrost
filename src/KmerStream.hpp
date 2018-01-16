#ifndef BFG_KMERSTREAM_HPP
#define BFG_KMERSTREAM_HPP

#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <sys/stat.h>
#include <fstream>
#include <stdint.h>
#include <sstream>
#include <bitset>

#include <thread>
#include <atomic>

#include "fastq.hpp"
#include "RepHash.hpp"
#include "StreamCounter.hpp"

using namespace std;

struct KmerStream_Build_opt {

    vector<int> klist;
    vector<size_t> q_cutoff;
    vector<string> files;
    bool verbose;
    double e;
    size_t q_base;
    size_t threads;
    size_t chunksize;

    KmerStream_Build_opt() :  verbose(false), e(0.01), threads(1), chunksize(100000), q_base(33) {}
};

class ReadQualityHasher;

class ReadHasher {

    friend class ReadQualityHasher;

    public:

        ReadHasher(const double e) : k(0), hf(), sc(e) {}

        ReadHasher(const ReadHasher& o) : k(o.k), hf(o.hf), sc(o.sc) {}

        void setK(const size_t _k) {

            k = _k;
            hf.setK(k);
        }

        // create hashes for all k-mers
        // operate on hashes
        void operator()(const char *s, const size_t l, const char *q, const size_t ql) {

            size_t i = 0, j = 0;
            bool last_valid = false;

            if (l < k) return;

            while (j < l) { // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed

                const char c = s[j];

                if (c != 'N' && c != 'n') {

                    if (last_valid) {
                        // s[i..j-1] was a valid k-mer k-mer, update
                        hf.update(s[i],s[j]);
                        ++i;
                    }
                    else if (i + k -1 == j) {

                        hf.init(s+i); // start the k-mer at position i
                        last_valid = true;
                    }

                    ++j;
                }
                else { // invalid character, restart

                    ++j;
                    i = j;
                    last_valid = false;
                }

                if (last_valid) handle(hf.hash());
            }
        }

        inline void handle(const uint64_t val) { sc(val); }

        inline void setQualityCutoff(const size_t q) {}

        inline string report() { return sc.report(); }

        inline string humanReport() { return sc.humanReport(); }

        inline string report(const bool tsv) { return sc.report(tsv); }

        inline size_t F0() const { return sc.F0(); }

        inline size_t F1() const { return sc.F1(); }

        inline size_t f1() const { return sc.f1(); }

        inline bool join(const ReadHasher& o) { return sc.join(o.sc); }

        inline bool join(const ReadQualityHasher& o);

        inline bool writeBinary(const std::string& fn) { return sc.writeBinary(fn); }

    private:

        size_t k;
        RepHash hf;
        StreamCounter sc;
};

class ReadQualityHasher {

    friend class ReadHasher;

    public:

        ReadQualityHasher(const double e_, const size_t q_base_) : k(0), hf(), sc(e_), q_cutoff(0), q_base(q_base_) {}

        ReadQualityHasher(const ReadQualityHasher& o) : k(o.k), q_cutoff(o.q_cutoff), q_base(o.q_base), hf(o.hf), sc(o.sc) {}

        void setK(size_t _k) {

            k = _k;
            hf.setK(k);
        }

        // create hashes for all k-mers
        // operate on hashes
        void operator()(const char* s, const size_t l, const char* q, const size_t ql) {

            size_t i = 0, j = 0;
            bool last_valid = false;

            if (l < k) return;

            while (j < l) {
                // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed
                const char c = s[j];

                if (c != 'N' && c != 'n' && (q[j] >= (char) (q_base+q_cutoff))) {

                    if (last_valid) {
                        // s[i..j-1] was a valid k-mer k-mer, update
                        hf.update(s[i],s[j]);
                        ++i;
                    }
                    else if (i + k -1 == j) {

                        hf.init(s+i); // start the k-mer at position i
                        last_valid = true;
                    }

                    ++j; // move along
                }
                else {
                    // invalid character, restart
                    ++j;
                    i = j;
                    last_valid = false;
                }

                if (last_valid) handle(hf.hash());
            }
        }

        inline void handle(const uint64_t val) { sc(val); }

        inline void setQualityCutoff(const size_t q) { q_cutoff = q; }

        inline string humanReport() { return sc.humanReport(); }

        inline string report(const bool tsv) { return sc.report(tsv); }

        inline size_t F0() const { return sc.F0(); }

        inline size_t F1() const { return sc.F1(); }

        inline size_t f1() const { return sc.f1(); }

        inline bool join(const ReadQualityHasher& o) { return sc.join(o.sc); }

        inline bool join(const ReadHasher& o);

        inline bool writeBinary(const std::string& fn) { return sc.writeBinary(fn); }

    private:

        size_t q_cutoff, q_base;
        size_t k;
        RepHash hf;
        StreamCounter sc;
};

inline bool ReadQualityHasher::join(const ReadHasher& o) { return sc.join(o.sc); }

inline bool ReadHasher::join(const ReadQualityHasher& o) { return sc.join(o.sc); }

class KmerStream {

    public:

        class KmerStream_const_iterator : public std::iterator<std::forward_iterator_tag, pair<size_t, size_t>>{

            private:

                size_t it;
                size_t sz;
                bool invalid;
                const KmerStream* kms;

            public:

                KmerStream_const_iterator() : kms(nullptr), it(0), sz(0), invalid(true) {}

                KmerStream_const_iterator(const KmerStream* kms_) : kms(kms_), it(0), sz(0), invalid(false) {

                    if ((kms == nullptr) || kms->invalid) invalid = true;
                    else sz = kms->klist.size() * kms->q_cutoff.size();
                }

                KmerStream_const_iterator(const KmerStream_const_iterator& o) : kms(o.kms), it(o.it), sz(o.sz), invalid(o.invalid) {}

                KmerStream_const_iterator& operator=(const KmerStream_const_iterator& o) {

                    kms = o.kms;
                    it = o.it;
                    sz = o.sz;
                    invalid = o.invalid;

                    return *this;
                }

                KmerStream_const_iterator operator++(int) {

                    KmerStream_const_iterator tmp(*this);
                    operator++();
                    return tmp;
                }

                KmerStream_const_iterator& operator++() {

                    if (!invalid){

                        if (it + 1 >= sz) invalid = true;
                        else ++it;
                    }

                    return *this;
                }

                bool operator==(const KmerStream_const_iterator& o) {

                    if (invalid || o.invalid) return invalid && o.invalid;
                    return  (it == o.it) && (kms == o.kms);
                }

                bool operator!=(const KmerStream_const_iterator& o) { return !operator==(o); }

                pair<size_t, size_t> operator*() const { return make_pair(F0(), f1()); }

                size_t F0() const {

                    if (invalid) return 0;
                    if (kms->v_fa_res.size() == 0) return kms->v_fq_res[it].F0();
                    return kms->v_fa_res[it].F0();
                }

                size_t F1() const {

                    if (invalid) return 0;
                    if (kms->v_fa_res.size() == 0) return kms->v_fq_res[it].F1();
                    return kms->v_fa_res[it].F1();
                }

                size_t f1() const {

                    if (invalid) return 0;
                    if (kms->v_fa_res.size() == 0) return kms->v_fq_res[it].f1();
                    return kms->v_fa_res[it].f1();
                }
        };

        typedef KmerStream_const_iterator const_iterator;

        KmerStream(const KmerStream_Build_opt& opt) :   klist(opt.klist), q_cutoff(opt.q_cutoff), verbose(opt.verbose), e(opt.e),
                                                        q_base(opt.q_base), nb_threads(opt.threads), chunksize(opt.chunksize), invalid(false) {

            size_t max_threads = std::thread::hardware_concurrency();

            if (nb_threads <= 0){

                cerr << "KmerStream::KmerStream(): Number of threads cannot be less than or equal to 0" << endl;
                invalid = true;
            }

            if (nb_threads > max_threads){

                cerr << "KmerStream::KmerStream(): Number of threads cannot be greater than or equal to " << max_threads << endl;
                invalid = true;
            }

            if (klist.size() == 0){

                cerr << "KmerStream::KmerStream(): No length k of k-mers given" << endl;
                invalid = true;
            }
            else {

                for (const auto& k : klist){

                    if (k <= 0){

                        cerr << "KmerStream::KmerStream(): Length k of k-mers cannot be less than or equal to 0" << endl;
                        invalid = true;
                    }
                }
            }

            if (q_cutoff.size() == 0){

                cerr << "KmerStream::KmerStream(): No quality cutoff given" << endl;
                invalid = true;
            }

            if (e <= 0){

                cerr << "KmerStream::KmerStream(): Guaranteed error rate cannot be less than or equal to 0" << endl;
                invalid = true;
            }

            if ((q_base != 33) && (q_base != 64)){

                cerr << "KmerStream::KmerStream(): Quality score can only be PHREAD+64 (q_base=64) or PHREAD+33 (q_base=33)" << endl;
                invalid = true;
            }

            if (opt.files.size() == 0) {

                cerr << "KmerStream::KmerStream(): Missing FASTA/FASTQ input files" << endl;
                invalid = true;
            }
            else {

                struct stat stFileInfo;
                int intStat;

                for (const auto& s : opt.files) {

                    intStat = stat(s.c_str(), &stFileInfo);

                    if (intStat != 0) {

                        cerr << "KmerStream::KmerStream(): File not found, " << s << endl;
                        invalid = true;
                    }
                    else {

                        const string s_ext = s.substr(s.find_last_of(".") + 1);

                        if ((s_ext == "fasta") || (s_ext == "fasta.gz") || (s_ext == "fa") || (s_ext == "fa.gz")) {

                            files_fasta.push_back(s);
                        }
                        else if ((s_ext == "fastq") || (s_ext == "fastq.gz") || (s_ext == "fq") || (s_ext == "fq.gz")) {

                            files_fastq.push_back(s);
                        }
                        else {

                            cerr << "KmerStream::KmerStream(): Input files must be in FASTA (*.fasta, *.fa, *.fasta.gz, *.fa.gz) or " <<
                            "FASTQ format (*.fastq, *.fq, *.fastq.gz, *.fq.gz)" << endl;

                            invalid = true;
                        }
                    }
                }
            }

            if (invalid) exit(1);

            if (files_fastq.size() != 0) v_fq_res = nb_threads > 1 ? RunThreadedFastqStream() : RunFastqStream();
            if (files_fasta.size() != 0) v_fa_res = nb_threads > 1 ? RunThreadedFastaStream() : RunFastaStream();

            size_t i = 0;

            if ((files_fastq.size() != 0) && (files_fasta.size() != 0)){

                for (auto& sp : v_fa_res){

                    sp.join(v_fq_res[i]);
                    i++;
                }

                v_fq_res.clear();
            }
        }

        const_iterator begin() { return const_iterator(this); }

        const_iterator end() { return const_iterator(); }

    private:

        vector<ReadQualityHasher> RunFastqStream() {

            std::ios_base::sync_with_stdio(false);

            FastqFile FQ(files_fastq);

            const kseq_t* kseq = FQ.get_kseq();

            size_t nreads = 0;

            const size_t qsize = q_cutoff.size();
            const size_t ksize = klist.size();

            vector<ReadQualityHasher> sps(qsize * ksize, ReadQualityHasher(e, q_base));

            for (size_t i = 0; i < qsize; ++i) {

                for (size_t j = 0; j < ksize; ++j) {

                    sps[i * ksize + j].setQualityCutoff(q_cutoff[i]);
                    sps[i * ksize + j].setK(klist[j]);
                }
            }

            while (FQ.read_next() >= 0){

                ++nreads;
                for (size_t i = 0; i < qsize * ksize; ++i) sps[i](kseq->seq.s, kseq->seq.l, kseq->qual.s, kseq->qual.l);
            }

            FQ.close();

            return sps;
        }

        vector<ReadQualityHasher> RunThreadedFastqStream() {

            std::ios_base::sync_with_stdio(false);

            size_t readCount;
            const size_t chunk = chunksize;
            const size_t qsize = q_cutoff.size();
            const size_t ksize = klist.size();
            const size_t nb_k_q = qsize * ksize;
            const size_t threads = nb_threads / nb_k_q;

            read_fastq_t reads;

            vector<ReadQualityHasher> sps(nb_k_q * threads, ReadQualityHasher(e, q_base));

            for (size_t i = 0; i < qsize; ++i) {

                for (size_t j = 0; j < ksize; ++j) {

                    for (size_t t = 0; t < threads; ++t){

                        sps[i * ksize + j * threads + t].setQualityCutoff(q_cutoff[i]);
                        sps[i * ksize + j * threads + t].setK(klist[j]);
                    }
                }
            }

            auto worker_function = [](read_fastq_t::const_iterator it_start, read_fastq_t::const_iterator it_end, ReadQualityHasher* sp) {

                for (auto it = it_start; it != it_end; ++it){

                    const pair<string, string>& p = *it;

                    (*sp)(p.first.c_str(), p.first.size(), p.second.c_str(), p.second.size());
                }
            };

            FastqFile FQ(files_fastq);

            const kseq_t* kseq = FQ.get_kseq();

            int round = 0;

            bool done = false;

            while (!done) {

                size_t readCount = 0;

                while (readCount < chunk) {

                    if (FQ.read_next() >= 0){

                        reads.push_back(make_pair(string(kseq->seq.s), string(kseq->qual.s)));
                        ++readCount;
                    }
                    else {

                        done = true;
                        break;
                    }
                }

                ++round;

                if (verbose) cout << "KmerStream::RunThreadedFastqStream(): Starting round " << round << endl;

                // run parallel code
                vector<thread> workers;

                const size_t batch_size = reads.size() / threads;
                const size_t leftover = reads.size() % threads;

                for (size_t i = 0; i < nb_k_q; ++i){

                    auto rit = reads.begin();

                    for (size_t j = 0; j < threads; ++j) {

                        const size_t jump = batch_size + (j < leftover ? 1 : 0);
                        auto rit_end(rit);

                        advance(rit_end, jump);
                        workers.push_back(thread(worker_function, rit, rit_end, &sps[i * threads + j]));

                        rit = rit_end;
                    }
                }

                for (auto& t : workers) t.join();

                reads.clear();
            }

            FQ.close();

            vector<ReadQualityHasher> res;

            for (size_t j = 0; j < nb_k_q; ++j){

                ReadQualityHasher& sp = sps[j * threads];

                for (size_t i = 1; i < threads; ++i) sp.join(sps[j * threads + i]);

                res.push_back(sp);
            }

            return res;
        };

        vector<ReadHasher> RunFastaStream() {

            std::ios_base::sync_with_stdio(false);

            FastqFile FQ(files_fasta);

            const kseq_t* kseq = FQ.get_kseq();

            size_t nreads = 0;

            const size_t qsize = q_cutoff.size();
            const size_t ksize = klist.size();

            vector<ReadHasher> sps(qsize * ksize, ReadHasher(e));

            for (size_t i = 0; i < qsize; ++i) {

                for (size_t j = 0; j < ksize; ++j) sps[i * ksize + j].setK(klist[j]);
            }

            while (FQ.read_next() >= 0){

                ++nreads;
                // seq->seq.s is of length seq->seq.l
                for (size_t i = 0; i < qsize * ksize; ++i) sps[i](kseq->seq.s, kseq->seq.l, nullptr, 0);
            }

            FQ.close();

            return sps;
        }

        vector<ReadHasher> RunThreadedFastaStream() {

            std::ios_base::sync_with_stdio(false);

            size_t readCount;
            const size_t chunk = chunksize;
            const size_t qsize = q_cutoff.size();
            const size_t ksize = klist.size();
            const size_t nb_k_q = qsize * ksize;
            const size_t threads = nb_threads / nb_k_q;

            read_fasta_t reads;

            vector<ReadHasher> sps(nb_k_q * threads, ReadHasher(e));

            for (size_t i = 0; i < qsize; ++i) {

                for (size_t j = 0; j < ksize; ++j) {

                    for (size_t t = 0; t < threads; ++t) sps[i * ksize + j * threads + t].setK(klist[j]);
                }
            }

            auto worker_function = [](read_fasta_t::const_iterator it_start, read_fasta_t::const_iterator it_end, ReadHasher* sp) {

                for (auto it = it_start; it != it_end; ++it) (*sp)(it->c_str(), it->size(), nullptr, 0);
            };

            FastqFile FQ(files_fasta);

            const kseq_t* kseq = FQ.get_kseq();

            int round = 0;

            bool done = false;

            while (!done) {

                size_t readCount = 0;

                while (readCount < chunk) {

                    if (FQ.read_next() >= 0){

                        reads.push_back(string(kseq->seq.s));
                        ++readCount;
                    }
                    else {

                        done = true;
                        break;
                    }
                }

                ++round;

                if (verbose) cout << "KmerStream::RunThreadedFastaStream(): Starting round " << round << endl;

                // run parallel code
                vector<thread> workers;

                const size_t batch_size = reads.size() / threads;
                const size_t leftover = reads.size() % threads;

                for (size_t i = 0; i < nb_k_q; ++i){

                    auto rit = reads.begin();

                    for (size_t j = 0; j < threads; ++j) {

                        const size_t jump = batch_size + (j < leftover ? 1 : 0);
                        auto rit_end(rit);

                        advance(rit_end, jump);
                        workers.push_back(thread(worker_function, rit, rit_end, &sps[i * threads + j]));

                        rit = rit_end;
                    }
                }

                for (auto& t : workers) t.join();

                reads.clear();
            }

            FQ.close();

            vector<ReadHasher> res;

            for (size_t j = 0; j < nb_k_q; ++j){

                ReadHasher& sp = sps[j * threads];

                for (size_t i = 1; i < threads; ++i) sp.join(sps[j * threads + i]);

                res.push_back(sp);
            }

            return res;
        }

        typedef vector<pair<string, string>> read_fastq_t;
        typedef vector<string> read_fasta_t;

        vector<int> klist;
        vector<size_t> q_cutoff;
        vector<string> files_fasta;
        vector<string> files_fastq;
        vector<ReadQualityHasher> v_fq_res;
        vector<ReadHasher> v_fa_res;
        bool verbose;
        bool invalid;
        double e;
        size_t q_base;
        size_t nb_threads;
        size_t chunksize;
};

#endif
