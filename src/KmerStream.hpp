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

    vector<string> files;

    bool verbose;

    double e;

    size_t k;
    size_t q;
    size_t q_base;
    size_t threads;
    size_t chunksize;

    KmerStream_Build_opt() : q_base(33), q(0), k(31), verbose(false), e(0.01), threads(1), chunksize(64) {}
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

                const char c = s[j] & 0xDF;

                if ((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) {

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

        void operator()(const char* seq_buf, const size_t seq_buf_sz) {

            const char* str = seq_buf;
            const char* str_end = &seq_buf[seq_buf_sz];

            while (str < str_end) { // for each input

                const int sl = strlen(str);

                if (sl >= k){

                    size_t i = 0, j = 0;

                    bool last_valid = false;

                    while (j < sl) { // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed

                        const char c = str[j] & 0xDF;

                        if ((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) {

                            if (last_valid) {
                                // s[i..j-1] was a valid k-mer k-mer, update
                                hf.update(str[i], str[j]);
                                ++i;
                            }
                            else if (i + k -1 == j) {

                                hf.init(str + i); // start the k-mer at position i
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

                str += sl + 1;
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

            if (l < k) return;

            size_t i = 0, j = 0;

            bool last_valid = false;

            const char q_base_cut = (char) (q_base + q_cutoff);

            while (j < l) {
                // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed
                const char c = s[j] & 0xDF;

                if (((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) && (q[j] >= q_base_cut)) {

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

        void operator()(const char* seq_buf, const char* qual_buf, const size_t buf_sz) {

            const char q_base_cut = (char) (q_base + q_cutoff);

            const char* str = seq_buf;
            const char* q_str = qual_buf;
            const char* str_end = &seq_buf[buf_sz];

            while (str < str_end) { // for each input

                const int sl = strlen(str);

                if (sl >= k){

                    size_t i = 0, j = 0;

                    bool last_valid = false;

                    while (j < sl) {
                        // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed
                        const char c = str[j] & 0xDF;

                        if (((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) && (q_str[j] >= q_base_cut)) {

                            if (last_valid) {
                                // s[i..j-1] was a valid k-mer k-mer, update
                                hf.update(str[i], str[j]);
                                ++i;
                            }
                            else if (i + k - 1 == j) {

                                hf.init(str + i); // start the k-mer at position i
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

                str += sl + 1;
                q_str += sl + 1;
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

        KmerStream(const KmerStream_Build_opt& opt) :   k(opt.k), q(opt.q), q_base(opt.q_base), e(opt.e), rqh(e, q_base), rsh(e),
                                                        nb_threads(opt.threads), chunksize(opt.chunksize), invalid(false), verbose(opt.verbose) {

            const size_t max_threads = std::thread::hardware_concurrency();

            if (nb_threads <= 0){

                cerr << "KmerStream::KmerStream(): Number of threads cannot be less than or equal to 0" << endl;
                invalid = true;
            }

            if (nb_threads > max_threads){

                cerr << "KmerStream::KmerStream(): Number of threads cannot be greater than or equal to " << max_threads << endl;
                invalid = true;
            }

            if (k == 0){

                cerr << "KmerStream::KmerStream(): Length k of k-mers cannot be less than or equal to 0" << endl;
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

                        const size_t last_point = s.find_last_of(".");

                        string s_ext = s.substr(last_point + 1);

                        if ((s_ext == "gz")){

                            s_ext = s.substr(s.find_last_of(".", last_point - 1) + 1);

                            if ((s_ext == "fasta.gz") || (s_ext == "fa.gz") || (s_ext == "fna.gz")) files_no_quality.push_back(s);
                            else if ((s_ext == "fastq.gz") || (s_ext == "fq.gz")) files_with_quality.push_back(s);
                            else {

                                cerr << "KmerStream::KmerStream(): Input files must be in FASTA (*.fasta, *.fa, *.fna, *.fasta.gz, *.fa.gz, *.fna.gz)" <<
                                " or FASTQ (*.fastq, *.fq, *.fastq.gz, *.fq.gz) or GFA (*.gfa) format" << endl;
                                cerr << "KmerStream::KmerStream(): Erroneous file is " << s << endl;

                                invalid = true;
                            }
                        }
                        else {

                            if ((s_ext == "fasta") || (s_ext == "fa") || (s_ext == "fna") || (s_ext == "gfa")) files_no_quality.push_back(s);
                            else if ((s_ext == "fastq") || (s_ext == "fq")) files_with_quality.push_back(s);
                            else {

                                cerr << "KmerStream::KmerStream(): Input files must be in FASTA (*.fasta, *.fa, *.fna, *.fasta.gz, *.fa.gz, *.fna.gz)" <<
                                " or FASTQ (*.fastq, *.fq, *.fastq.gz, *.fq.gz) or GFA (*.gfa) format" << endl;
                                cerr << "KmerStream::KmerStream(): Erroneous file is " << s << endl;

                                invalid = true;
                            }
                        }
                    }
                }
            }

            if (invalid) exit(1);

            rqh.setQualityCutoff(q);

            rqh.setK(k);
            rsh.setK(k);

            if (verbose) cout << "KmerStream::KmerStream(): Start computing k-mer cardinality estimations" << endl;

            if (files_with_quality.size() != 0) nb_threads > 1 ? RunThreadedFastqStream() : RunFastqStream();
            if (files_no_quality.size() != 0) nb_threads > 1 ? RunThreadedFastaStream() : RunFastaStream();

            rsh.join(rqh);
        }

        inline size_t F0() const { return rsh.F0(); }

        inline size_t F1() const { return rsh.F1(); }

        inline size_t f1() const { return rsh.f1(); }

    private:

        void RunFastqStream() {

            FileParser fp(files_with_quality);

            size_t file_id = 0;

            string seq;

            while (fp.read(seq, file_id)){

                const char* qss = fp.getQualityScoreString();

                std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

                rqh(seq.c_str(), seq.length(), qss, strlen(qss));
            }

            fp.close();
        }

        void RunThreadedFastqStream() {

            size_t pos_read = 0;
            size_t len_read = 0;

            const size_t max_len_seq = 1024;
            const size_t thread_seq_buf_sz = chunksize * max_len_seq;

            string seq, qual;

            vector<ReadQualityHasher> sps(nb_threads, rqh);

            FileParser fp(files_with_quality);

            auto reading_function = [&](char* seq_buf, char* qual_buf, size_t& buf_sz) {

                size_t file_id = 0;

                const size_t seq_buf_sz = thread_seq_buf_sz - k;

                const char* s_str = seq.c_str();
                const char* q_str = qual.c_str();

                buf_sz = 0;

                while (buf_sz < seq_buf_sz) {

                    const bool new_reading = (pos_read >= len_read);

                    if (!new_reading || fp.read(seq, file_id)) {

                        if (new_reading) qual = fp.getQualityScoreString();

                        //pos_read = (new_reading ? 0 : pos_read);
                        pos_read &= static_cast<size_t>(new_reading) - 1;

                        len_read = seq.length();
                        s_str = seq.c_str();
                        q_str = qual.c_str();

                        if (len_read >= k){

                            if ((thread_seq_buf_sz - buf_sz - 1) < (len_read - pos_read)){

                                strncpy(&seq_buf[buf_sz], &s_str[pos_read], thread_seq_buf_sz - buf_sz - 1);
                                strncpy(&qual_buf[buf_sz], &q_str[pos_read], thread_seq_buf_sz - buf_sz - 1);

                                seq_buf[thread_seq_buf_sz - 1] = '\0';

                                pos_read += seq_buf_sz - buf_sz;
                                buf_sz = thread_seq_buf_sz;

                                break;
                            }
                            else {

                                strcpy(&seq_buf[buf_sz], &s_str[pos_read]);
                                strcpy(&qual_buf[buf_sz], &q_str[pos_read]);

                                buf_sz += (len_read - pos_read) + 1;
                                pos_read = len_read;
                            }
                        }
                        else pos_read = len_read;
                    }
                    else return true;
                }

                return false;
            };

            {
                vector<thread> workers; // need to keep track of threads so we can join them

                bool stop = false;

                mutex mutex_file;

                char* buffer_seq = new char[nb_threads * thread_seq_buf_sz];
                char* buffer_qual = new char[nb_threads * thread_seq_buf_sz];
                size_t* buffer_sz = new size_t[nb_threads];

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&, t]{

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    if (stop) return;

                                    stop = reading_function(&buffer_seq[t * thread_seq_buf_sz], &buffer_qual[t * thread_seq_buf_sz], buffer_sz[t]);
                                }

                                for (char* s = &buffer_seq[t * thread_seq_buf_sz]; s != &buffer_seq[(t + 1) * thread_seq_buf_sz]; ++s) *s &= 0xDF;

                                sps[t](&buffer_seq[t * thread_seq_buf_sz], &buffer_qual[t * thread_seq_buf_sz], buffer_sz[t]);
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();

                delete[] buffer_seq;
                delete[] buffer_qual;
                delete[] buffer_sz;
            }

            fp.close();

            for (size_t t = 0; t != nb_threads; ++t) rqh.join(sps[t]);
        }

        void RunFastaStream() {

            size_t file_id = 0;

            string seq;

            FileParser fp(files_no_quality);

            while (fp.read(seq, file_id)){

                std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

                rsh(seq.c_str(), seq.length(), nullptr, 0);
            }

            fp.close();
        }

        void RunThreadedFastaStream() {

            size_t pos_read = 0;
            size_t len_read = 0;

            const size_t max_len_seq = 1024;
            const size_t thread_seq_buf_sz = chunksize * max_len_seq;

            vector<ReadHasher> sps(nb_threads, rsh);

            FileParser fp(files_no_quality);

            string s;

            auto reading_function = [&](char* seq_buf, size_t& seq_buf_sz) {

                size_t file_id = 0;

                const size_t sz_buf = thread_seq_buf_sz - k;

                const char* s_str = s.c_str();

                seq_buf_sz = 0;

                while (seq_buf_sz < sz_buf) {

                    const bool new_reading = (pos_read >= len_read);

                    if (!new_reading || fp.read(s, file_id)) {

                        //pos_read = (new_reading ? 0 : pos_read);
                        pos_read &= static_cast<size_t>(new_reading) - 1;

                        len_read = s.length();
                        s_str = s.c_str();

                        if (len_read >= k){

                            if ((thread_seq_buf_sz - seq_buf_sz - 1) < (len_read - pos_read)){

                                strncpy(&seq_buf[seq_buf_sz], &s_str[pos_read], thread_seq_buf_sz - seq_buf_sz - 1);

                                seq_buf[thread_seq_buf_sz - 1] = '\0';

                                pos_read += sz_buf - seq_buf_sz;
                                seq_buf_sz = thread_seq_buf_sz;

                                break;
                            }
                            else {

                                strcpy(&seq_buf[seq_buf_sz], &s_str[pos_read]);

                                seq_buf_sz += (len_read - pos_read) + 1;
                                pos_read = len_read;
                            }
                        }
                        else pos_read = len_read;
                    }
                    else return true;
                }

                return false;
            };

            {
                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mutex_file;

                bool stop = false;

                char* buffer_seq = new char[nb_threads * thread_seq_buf_sz];
                size_t* buffer_seq_sz = new size_t[nb_threads];

                for (size_t t = 0; t != nb_threads; ++t){

                    workers.emplace_back(

                        [&, t]{

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    if (stop) return;

                                    stop = reading_function(&buffer_seq[t * thread_seq_buf_sz], buffer_seq_sz[t]);
                                }

                                for (char* s = &buffer_seq[t * thread_seq_buf_sz]; s != &buffer_seq[(t + 1) * thread_seq_buf_sz]; ++s) *s &= 0xDF;

                                sps[t](&buffer_seq[t * thread_seq_buf_sz], buffer_seq_sz[t]);
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();

                delete[] buffer_seq;
                delete[] buffer_seq_sz;
            }

            fp.close();

            for (size_t t = 0; t != nb_threads; ++t) rsh.join(sps[t]);
        }

        size_t k;
        size_t q;
        size_t q_base;

        double e;

        ReadQualityHasher rqh;
        ReadHasher rsh;

        vector<string> files_no_quality;
        vector<string> files_with_quality;

        bool verbose;
        bool invalid;

        size_t nb_threads;
        size_t chunksize;
};

#endif
