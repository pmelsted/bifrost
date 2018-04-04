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

    KmerStream_Build_opt() : q_base(33), q(0), k(31), verbose(false), e(0.01), threads(1), chunksize(10000) {}
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
                const char c = s[j] & 0xDF;

                if (((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) && (q[j] >= (char) (q_base+q_cutoff))) {

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

                        string s_ext = s.substr(s.find_last_of(".") + 1);

                        if ((s_ext == "gz")){

                            s_ext = s_ext.substr(s_ext.find_last_of(".") + 1);

                            if ((s_ext == "fasta") || (s_ext == "fa")) files_no_quality.push_back(s);
                            else if ((s_ext == "fastq") || (s_ext == "fq")) files_with_quality.push_back(s);
                            else {

                                cerr << "KmerStream::KmerStream(): Input files must be in FASTA (*.fasta, *.fa, *.fasta.gz, *.fa.gz) or " <<
                                "FASTQ (*.fastq, *.fq, *.fastq.gz, *.fq.gz) or GFA (*.gfa) format" << endl;
                                cerr << "KmerStream::KmerStream(): Erroneous file is " << s << endl;

                                invalid = true;
                            }
                        }
                        else {

                            if ((s_ext == "fasta") || (s_ext == "fa") || (s_ext == "gfa")) files_no_quality.push_back(s);
                            else if ((s_ext == "fastq") || (s_ext == "fq")) files_with_quality.push_back(s);
                            else {

                                cerr << "KmerStream::KmerStream(): Input files must be in FASTA (*.fasta, *.fa, *.fasta.gz, *.fa.gz) or " <<
                                "FASTQ (*.fastq, *.fq, *.fastq.gz, *.fq.gz) or GFA (*.gfa) format" << endl;
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

            size_t file_id = 0;

            size_t pos_read = k - 1;
            size_t len_read = 0;

            string seq, qual;

            vector<ReadQualityHasher> sps(nb_threads, rqh);

            FileParser fp(files_with_quality);

            auto reading_function = [&](vector<pair<string, string>>& readv) {

                size_t reads_now = 0;

                while ((pos_read < len_read) && (reads_now < chunksize)){

                    pos_read -= k - 1;

                    readv.push_back(make_pair(seq.substr(pos_read, 1000), qual.substr(pos_read, 1000)));

                    pos_read += 1000;

                    ++reads_now;
                }

                while (reads_now < chunksize) {

                    if (fp.read(seq, file_id)) {

                        qual = fp.getQualityScoreString();

                        len_read = seq.length();
                        pos_read = len_read;

                        if (len_read > 1000){

                            pos_read = k - 1;

                            while ((pos_read < len_read) && (reads_now < chunksize)){

                                pos_read -= k - 1;

                                readv.push_back(make_pair(seq.substr(pos_read, 1000), qual.substr(pos_read, 1000)));

                                pos_read += 1000;

                                ++reads_now;
                            }
                        }
                        else {

                            readv.push_back(make_pair(seq, qual));

                            ++reads_now;
                        }
                    }
                    else {

                        for (auto& p : readv) std::transform(p.first.begin(), p.first.end(), p.first.begin(), ::toupper);

                        return true;
                    }
                }

                for (auto& p : readv) std::transform(p.first.begin(), p.first.end(), p.first.begin(), ::toupper);

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

                                for (const auto& r : readvs[t]) sps[t](r.first.c_str(), r.first.size(), r.second.c_str(), r.second.size());

                                readvs[t].clear();
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
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

            size_t file_id = 0;

            size_t pos_read = k - 1;
            size_t len_read = 0;

            string seq;

            vector<ReadHasher> sps(nb_threads, rsh);

            FileParser fp(files_no_quality);

            auto reading_function = [&](vector<string>& readv) {

                size_t reads_now = 0;

                while ((pos_read < len_read) && (reads_now < chunksize)){

                    pos_read -= k - 1;

                    readv.emplace_back(seq.substr(pos_read, 1000));

                    pos_read += 1000;

                    ++reads_now;
                }

                while (reads_now < chunksize) {

                    if (fp.read(seq, file_id)) {

                        len_read = seq.length();
                        pos_read = len_read;

                        if (len_read > 1000){

                            pos_read = k - 1;

                            while ((pos_read < len_read) && (reads_now < chunksize)){

                                pos_read -= k - 1;

                                readv.emplace_back(seq.substr(pos_read, 1000));

                                pos_read += 1000;

                                ++reads_now;
                            }
                        }
                        else {

                            readv.emplace_back(seq);

                            ++reads_now;
                        }
                    }
                    else {

                        for (auto& s : readv) std::transform(s.begin(), s.end(), s.begin(), ::toupper);

                        return true;
                    }
                }

                for (auto& s : readv) std::transform(s.begin(), s.end(), s.begin(), ::toupper);

                return false;
            };

            {
                vector<thread> workers; // need to keep track of threads so we can join them
                vector<vector<string>> readvs(nb_threads);

                bool stop = false;

                mutex mutex_file;

                for (size_t t = 0; t != nb_threads; ++t){

                    workers.emplace_back(

                        [&, t]{

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    if (stop) return;

                                    stop = reading_function(readvs[t]);
                                }

                                for (const auto& r : readvs[t]) sps[t](r.c_str(), r.size(), nullptr, 0);

                                readvs[t].clear();
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
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
