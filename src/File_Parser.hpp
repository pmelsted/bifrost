#ifndef BFG_FILE_PARSER_HPP
#define BFG_FILE_PARSER_HPP

#include "FASTX_Parser.hpp"
#include "GFA_Parser.hpp"

/*class FileParser {

    public:

        FileParser(const vector<string>& filenames) : file_no(0), invalid(false) {

            if (filenames.size() == 0) {

                cerr << "FileParser::FileParser(): Missing input files" << endl;
                invalid = true;
            }
            else {

                struct stat stFileInfo;

                for (const auto& s : filenames) {

                    const int intStat = stat(s.c_str(), &stFileInfo);

                    if (intStat != 0) {

                        cerr << "FileParser::FileParser(): File not found: " << s << endl;
                        invalid = true;
                    }
                    else {

                        string s_ext = s.substr(s.find_last_of(".") + 1);

                        if ((s_ext == "gz")){

                            s_ext = s_ext.substr(s_ext.find_last_of(".") + 1);

                            if ((s_ext == "fasta") || (s_ext == "fa") || (s_ext == "fastq") || (s_ext == "fq")) files_fastx.push_back(s);
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

                if (files_fastx.size() != 0) ff = FastqFile(files_fastx);
                else if (files_gfa.size() != 0) gfap = GFA_Parser(files_gfa);
                else invalid = true;
            }
        }

    private:

        bool invalid;

        size_t file_no;

        vector<string> files_fastx;
        vector<string> files_gfa;

        FastqFile ff;
        GFA_Parser gfap;
};*/

#endif
