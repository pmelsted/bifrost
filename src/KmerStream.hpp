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

#include "File_Parser.hpp"
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

                        cerr << "KmerStream::KmerStream(): File not found: " << s << endl;
                        invalid = true;
                    }
                    else {

                        string s_ext = s.substr(s.find_last_of(".") + 1);

                        if ((s_ext == "gz")){

                            s_ext = s_ext.substr(s_ext.find_last_of(".") + 1);

                            if ((s_ext == "fasta") || (s_ext == "fa")) files_no_quality.push_back(s);
                            else if ((s_ext == "fastq") || (s_ext == "fq")) files_with_quality.push_back(s);
                            else {

                                cerr << "KmerStream::KmerStream(): Input files must be in FASTA (*.fasta, *.fa, *.fasta.gz, *.fa.gz) or " <<
                                "FASTQ (*.fastq, *.fq, *.fastq.gz, *.fq.gz) or GFA (*.gfa) format" << endl;

                                invalid = true;
                            }
                        }
                        else {

                            if ((s_ext == "fasta") || (s_ext == "fa") || (s_ext == "gfa")) files_no_quality.push_back(s);
                            else if ((s_ext == "fastq") || (s_ext == "fq")) files_with_quality.push_back(s);
                            else {

                                cerr << "KmerStream::KmerStream(): Input files must be in FASTA (*.fasta, *.fa, *.fasta.gz, *.fa.gz) or " <<
                                "FASTQ (*.fastq, *.fq, *.fastq.gz, *.fq.gz) or GFA (*.gfa) format" << endl;

                                invalid = true;
                            }
                        }
                    }
                }
            }

            if (invalid) exit(1);

            if (files_with_quality.size() != 0) v_fq_res = nb_threads > 1 ? RunThreadedFastqStream() : RunFastqStream();
            if (files_no_quality.size() != 0) v_fa_res = nb_threads > 1 ? RunThreadedFastaStream() : RunFastaStream();

            size_t i = 0;

            if ((files_with_quality.size() != 0) && (files_no_quality.size() != 0)){

                for (auto& sp : v_fa_res){

                    sp.join(v_fq_res[i]);
                    ++i;
                }

                v_fq_res.clear();
            }
        }

        const_iterator begin() { return const_iterator(this); }

        const_iterator end() { return const_iterator(); }

    private:

        vector<ReadQualityHasher> RunFastqStream() {

            //std::ios_base::sync_with_stdio(false);

            FileParser fp(files_with_quality);

            size_t nreads = 0;
            size_t file_id = 0;

            const size_t qsize = q_cutoff.size();
            const size_t ksize = klist.size();

            string seq;

            vector<ReadQualityHasher> sps(qsize * ksize, ReadQualityHasher(e, q_base));

            for (size_t i = 0; i < qsize; ++i) {

                for (size_t j = 0; j < ksize; ++j) {

                    sps[i * ksize + j].setQualityCutoff(q_cutoff[i]);
                    sps[i * ksize + j].setK(klist[j]);
                }
            }

            while (fp.read(seq, file_id)){

                ++nreads;

                const char* qss = fp.getQualityScoreString();

                for (size_t i = 0; i < qsize * ksize; ++i) sps[i](seq.c_str(), seq.length(), qss, strlen(qss));
            }

            fp.close();

            return sps;
        }

        vector<ReadQualityHasher> RunThreadedFastqStream() {

            //std::ios_base::sync_with_stdio(false);

            size_t file_id = 0;

            const size_t chunk_size = 1000;

            const size_t qsize = q_cutoff.size();
            const size_t ksize = klist.size();
            const size_t nb_k_q = qsize * ksize;

            string seq;

            vector<ReadQualityHasher> sps(nb_k_q * nb_threads, ReadQualityHasher(e, q_base));

            for (size_t i = 0; i < qsize; ++i) {

                for (size_t j = 0; j < ksize; ++j) {

                    for (size_t t = 0; t < nb_threads; ++t){

                        sps[(i * ksize + j) * nb_threads + t].setQualityCutoff(q_cutoff[i]);
                        sps[(i * ksize + j) * nb_threads + t].setK(klist[j]);
                    }
                }
            }

            FileParser fp(files_with_quality);

            auto reading_function = [&](vector<pair<string, string>>& readv) {

                size_t readCount = 0;

                const size_t chunk = chunk_size * 100;

                while (readCount < chunk) {

                    if (fp.read(seq, file_id)){

                        readv.push_back(make_pair(seq, string(fp.getQualityScoreString())));

                        readCount += seq.length();
                    }
                    else return true;
                }

                return false;
            };

            {
                vector<thread> workers; // need to keep track of threads so we can join them
                vector<vector<pair<string, string>>> readvs(nb_threads);

                bool stop = false;

                mutex mutex_file;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&, t]{

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    if (stop) return;

                                    stop = reading_function(readvs[t]);
                                }

                                for (size_t i = 0; i < qsize; ++i) {

                                    for (size_t j = 0; j < ksize; ++j){

                                        const size_t id = (i * ksize + j) * nb_threads + t;

                                        for (const auto& r : readvs[t]) sps[id](r.first.c_str(), r.first.size(), r.second.c_str(), r.second.size());
                                    }
                                }

                                readvs[t].clear();
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }

            fp.close();

            vector<ReadQualityHasher> res;

            for (size_t j = 0; j < nb_k_q; ++j){

                ReadQualityHasher& sp = sps[j * nb_threads];

                for (size_t i = 1; i < nb_threads; ++i) sp.join(sps[j * nb_threads + i]);

                res.push_back(sp);
            }

            return res;
        }

        vector<ReadHasher> RunFastaStream() {

            //std::ios_base::sync_with_stdio(false);

            size_t nreads = 0;
            size_t file_id = 0;

            const size_t qsize = q_cutoff.size();
            const size_t ksize = klist.size();

            string seq;

            vector<ReadHasher> sps(qsize * ksize, ReadHasher(e));

            for (size_t i = 0; i < qsize; ++i) {

                for (size_t j = 0; j < ksize; ++j) sps[i * ksize + j].setK(klist[j]);
            }

            FileParser fp(files_no_quality);

            while (fp.read(seq, file_id)){

                ++nreads;

                for (size_t i = 0; i < qsize * ksize; ++i) sps[i](seq.c_str(), seq.length(), nullptr, 0);
            }

            fp.close();

            return sps;
        }

        vector<ReadHasher> RunThreadedFastaStream() {

            //std::ios_base::sync_with_stdio(false);

            size_t file_id;

            const size_t chunk_size = 1000;

            const size_t qsize = q_cutoff.size();
            const size_t ksize = klist.size();
            const size_t nb_k_q = qsize * ksize;

            string seq;

            vector<ReadHasher> sps(nb_k_q * nb_threads, ReadHasher(e));

            for (size_t i = 0; i < qsize; ++i) {

                for (size_t j = 0; j < ksize; ++j) {

                    for (size_t t = 0; t < nb_threads; ++t) sps[(i * ksize + j) * nb_threads + t].setK(klist[j]);
                }
            }

            FileParser fp(files_no_quality);

            auto reading_function = [&](vector<string>& readv) {

                size_t readCount = 0;

                const size_t chunk = chunk_size * 100;

                while (readCount < chunk) {

                    if (fp.read(seq, file_id)){

                        readv.push_back(seq);

                        readCount += seq.length();
                    }
                    else {

                        return true;
                    }
                }

                return false;
            };

            {
                vector<thread> workers; // need to keep track of threads so we can join them
                vector<vector<string>> readvs(nb_threads);

                bool stop = false;

                mutex mutex_file;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&, t]{

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    if (stop) return;

                                    stop = reading_function(readvs[t]);
                                }

                                for (size_t i = 0; i < qsize; ++i) {

                                    for (size_t j = 0; j < ksize; ++j){

                                        const size_t id = (i * ksize + j) * nb_threads + t;

                                        for (const auto& r : readvs[t]) sps[id](r.c_str(), r.size(), nullptr, 0);
                                    }
                                }

                                readvs[t].clear();
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }

            fp.close();

            vector<ReadHasher> res;

            for (size_t j = 0; j < nb_k_q; ++j){

                ReadHasher& sp = sps[j * nb_threads];

                for (size_t i = 1; i < nb_threads; ++i) sp.join(sps[j * nb_threads + i]);

                res.push_back(sp);
            }

            return res;
        }

        vector<int> klist;
        vector<size_t> q_cutoff;
        vector<string> files_no_quality;
        vector<string> files_with_quality;
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
