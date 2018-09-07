#ifndef BFG_FILE_PARSER_HPP
#define BFG_FILE_PARSER_HPP

#include <sstream>

#include "FASTX_Parser.hpp"
#include "GFA_Parser.hpp"

class FileParser {

    public:

        FileParser(const vector<string>& filenames) :   files_it(0), files_fastx_it(0), files_gfa_it(0),
                                                        reading_fastx(false), invalid(false) {

            if (filenames.size() == 0) {

                cerr << "FileParser::FileParser(): Missing input files" << endl;
                invalid = true;
            }
            else {

                struct stat stFileInfo;

                files = filenames;

                for (const auto& s : files) {

                    const int intStat = stat(s.c_str(), &stFileInfo);

                    if (intStat != 0) {

                        cerr << "FileParser::FileParser(): File not found: " << s << endl;
                        invalid = true;
                    }
                    else {

                        const size_t last_point = s.find_last_of(".");

                        string s_ext = s.substr(last_point + 1);

                        if ((s_ext == "gz")){

                            s_ext = s.substr(s.find_last_of(".", last_point - 1) + 1);

                            if ((s_ext == "fasta.gz") || (s_ext == "fa.gz") || (s_ext == "fastq.gz") || (s_ext == "fq.gz")) files_fastx.push_back(s);
                            else {

                                cerr << "FileParser::FileParser(): Input files must be in FASTA (*.fasta, *.fa, *.fasta.gz, *.fa.gz) or " <<
                                "FASTQ (*.fastq, *.fq, *.fastq.gz, *.fq.gz) or GFA (*.gfa) format" << endl;

                                invalid = true;
                            }
                        }
                        else {

                            if ((s_ext == "fasta") || (s_ext == "fa") || (s_ext == "fastq") || (s_ext == "fq")) files_fastx.push_back(s);
                            else if (s_ext == "gfa") files_gfa.push_back(s);
                            else {

                                cerr << "FileParser::FileParser(): Input files must be in FASTA (*.fasta, *.fa, *.fasta.gz, *.fa.gz) or " <<
                                "FASTQ (*.fastq, *.fq, *.fastq.gz, *.fq.gz) or GFA (*.gfa) format" << endl;

                                invalid = true;
                            }
                        }
                    }
                }
            }

            if (!invalid){

                if (files_fastx.size() != 0){

                    ff = FastqFile(files_fastx);
                    reading_fastx = (files[0] == files_fastx[0]);
                }

                if (files_gfa.size() != 0){

                    gfap = GFA_Parser(files_gfa);
                    invalid = !gfap.open_read();
                    reading_fastx = (files[0] != files_gfa[0]);
                }
            }
        }

        bool read(string& seq, size_t& file_id){

            if (!invalid){

                bool new_file;

                if (reading_fastx){

                    const int ret = ff.read_next(seq, files_fastx_it, new_file);

                    if (new_file || (ret == -1)){ // Need to open next file

                        invalid = ((files_it + 1) >= files.size()); // Invalidate this object if no more file to read

                        if (!invalid){ // If still some files to read

                            ++files_it; //Increment iterator to next file to read

                            // Check if next file to read is FASTA/FASTQ format
                            reading_fastx = ((ret != -1) && (files[files_it] == files_fastx[files_fastx_it]));

                            return read(seq, file_id); // We read the next line of the file
                        }
                    }
                }
                else {
                    // Read first line of next GFA file, skip edge lines
                    GFA_Parser::GFA_line r = gfap.read(files_gfa_it, new_file, true);

                    if (new_file || ((r.first == nullptr) && (r.second == nullptr))){ // Need to open next file

                        invalid = ((files_it + 1) >= files.size()); // Invalidate this object if no more file to read

                        if (!invalid){ // If still some files to read

                            ++files_it; //Increment iterator to next file to read

                            // Check if next file to read is FASTA/FASTQ format
                            reading_fastx = (files[files_it] == files_fastx[files_fastx_it]);

                            return read(seq, file_id); // We read the next line of the file
                        }
                    }

                    if (r.first != nullptr) seq = r.first->seq;
                }

                file_id = files_it;
            }

            return !invalid;
        }

        bool read(stringstream& ss, size_t& file_id){

            if (!invalid){

                bool new_file;

                if (reading_fastx){

                    const int ret = ff.read_next(ss, files_fastx_it, new_file);

                    if (new_file || (ret == -1)){ // Need to open next file

                        invalid = ((files_it + 1) >= files.size()); // Invalidate this object if no more file to read

                        if (!invalid){ // If still some files to read

                            ++files_it; //Increment iterator to next file to read

                            // Check if next file to read is FASTA/FASTQ format
                            reading_fastx = ((ret != -1) && (files[files_it] == files_fastx[files_fastx_it]));

                            return read(ss, file_id); // We read the next line of the file
                        }
                    }
                }
                else {
                    // Read first line of next GFA file, skip edge lines
                    GFA_Parser::GFA_line r = gfap.read(files_gfa_it, new_file, true);

                    if (new_file || ((r.first == nullptr) && (r.second == nullptr))){ // Need to open next file

                        invalid = ((files_it + 1) >= files.size()); // Invalidate this object if no more file to read

                        if (!invalid){ // If still some files to read

                            ++files_it; //Increment iterator to next file to read

                            // Check if next file to read is FASTA/FASTQ format
                            reading_fastx = (files[files_it] == files_fastx[files_fastx_it]);

                            return read(ss, file_id); // We read the next line of the file
                        }
                    }

                    if (r.first != nullptr) ss << r.first->seq;
                }

                file_id = files_it;
            }

            return !invalid;
        }

        const char* getQualityScoreString() const {

            if (invalid || !reading_fastx) return nullptr;
            return ff.get_kseq()->qual.s;
        }

        void close(){

            ff.close();
            gfap.close();
        }

    private:

        bool invalid;
        bool reading_fastx;

        size_t files_it;
        size_t files_fastx_it;
        size_t files_gfa_it;

        vector<string> files;
        vector<string> files_fastx;
        vector<string> files_gfa;

        FastqFile ff;
        GFA_Parser gfap;
};

#endif
