#ifndef BFG_GFA_PARSER_HPP
#define BFG_GFA_PARSER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <sys/stat.h>
#include <fstream>
#include <stdint.h>
#include <sstream>

class GFA_Parser {

    struct Sequence {

        string id;
        string seq;
        size_t len;
        int64_t coverage;

        Sequence() : id(), seq(), len(0), coverage(-1) {};
        Sequence(const string& id_, const string& seq_, const size_t len_, const int64_t coverage_) :   id(id_), seq(seq_), len(len_),
                                                                                                        coverage(coverage_) {};
    };

    struct Edge {

        int64_t edge_id;

        string vertexA_id;
        size_t pos_start_overlapA;
        size_t pos_end_overlapA;
        bool strand_overlapA;

        string vertexB_id;
        size_t pos_start_overlapB;
        size_t pos_end_overlapB;
        bool strand_overlapB;

        Edge() :    edge_id(-1), vertexA_id(), pos_start_overlapA(0), pos_end_overlapA(0), strand_overlapA(true),
                    vertexB_id(), pos_start_overlapB(0), pos_end_overlapB(0), strand_overlapB(true) {};

        Edge(const string vertexA_id_, const size_t pos_start_overlapA_, const size_t pos_end_overlapA_, const bool strand_overlapA_,
             const string vertexB_id_, const size_t pos_start_overlapB_, const size_t pos_end_overlapB_, const bool strand_overlapB_,
             const int64_t edge_id_ = -1) : edge_id(edge_id_), vertexA_id(vertexA_id_), pos_start_overlapA(pos_start_overlapA_),
             pos_end_overlapA(pos_end_overlapA_), strand_overlapA(strand_overlapA_), vertexB_id(vertexB_id_), pos_start_overlapB(pos_start_overlapB_),
             pos_end_overlapB(pos_end_overlapB_), strand_overlapB(strand_overlapB_) {};
    };

    struct GFA_line {

        const Sequence* seq;
        const Edge* edge;

        GFA_line() : seq(nullptr), edge(nullptr) {};
        GFA_line(const Sequence* seq_, const Edge* edge_) : seq(seq_), edge(edge_) {};
    };

    public:

        GFA_Parser(const string& filename, const size_t buffer_size = 100000) : file_open_write(false), file_open_read(false),
                                                                                graph_out(nullptr), graph_in(nullptr),
                                                                                buff_sz(buffer_size), buffer(nullptr) {

            graph_filename = filename;

            size_t pos_match_point = filename.find_last_of(".");

            if (pos_match_point == string::npos) graph_filename.append(".gfa");
            else {

                const string ext = filename.substr(pos_match_point + 1);

                if (ext != "gfa") graph_filename.append(".gfa");
            }

            buffer = new char[buff_sz]();
        }

        ~GFA_Parser() {

            if (buffer != nullptr) delete[] buffer;
        }

        bool open_write(const size_t version_GFA) {

            FILE* fp = fopen(graph_filename.c_str(), "w");

            if ((file_open_write = (fp != NULL)) == true) {

                fclose(fp);

                if (remove(graph_filename.c_str()) != 0) cerr << "GFA_Parser::open_write(): Could not remove temporary file " << graph_filename << endl;
            }
            else cerr << "GFA_Parser::open_write(): Could not open file " << graph_filename << " for writing graph" << endl;

            if (version_GFA > 2){

                cerr << "GFA_Parser::open_write(): Only supports GFA format version 1 and 2" << endl;
                file_open_write = false;
            }
            else v_gfa = version_GFA;

            if (file_open_write){

                graphfile.open(graph_filename.c_str());
                graph_out.rdbuf(graphfile.rdbuf());

                graph_out << "H\tVN:Z:" << (v_gfa == 1 ? "1" : "2") << ".0\n";
            }

            return file_open_write;
        }

        bool open_read() {

            FILE* fp = fopen(graph_filename.c_str(), "r");

            if ((file_open_read = (fp != NULL)) == true) fclose(fp);
            else cerr << "GFA_Parser::open_read(): Could not open file " << graph_filename << " for reading graph" << endl;

            if (file_open_read) {

                graphfile.open(graph_filename.c_str());
                graph_in.rdbuf(graphfile.rdbuf());

                graph_in >> buffer;

                const string header(buffer);

                if (header == "H\tVN:Z:1.0\n") v_gfa = 1;
                else if (header == "H\tVN:Z:2.0\n") v_gfa = 2;
                else {

                    cerr << "GFA_Parser::open_read(): Wrong GFA header or unsupported GFA format version (GFA_Parser only supports version 1 and 2)" << endl;
                    close();
                }
            }

            return file_open_read;
        }

        void close(){

            if (file_open_write || file_open_read){

                graphfile.close();

                file_open_write = false;
                file_open_read = false;
            }
        }

        bool write_sequence(const string& seq, const string& id, const size_t len, const int64_t coverage){

            if (file_open_write){

                if (v_gfa == 1) graph_out << "S" << "\t" << id << "\t" << seq << "\t" << "LN:i:" << len << "\t" << "XC:i:" << coverage << "\n";
                else graph_out << "S" << "\t" << id << "\t" << len << "\t" << seq << "\t" << "XC:i:" << coverage << "\n";
            }
            else {

                cerr << "GFA_Parser::write_sequence(): Input file is not open in writing mode" << endl;
                return false;
            }

            return true;
        }

        bool write_edge(const string vertexA_id, const size_t pos_start_overlapA, const size_t pos_end_overlapA, const bool strand_overlapA,
                        const string vertexB_id, const size_t pos_start_overlapB, const size_t pos_end_overlapB, const bool strand_overlapB,
                        const int64_t edge_id = -1) {

            if (file_open_write){

                if (pos_start_overlapA >= pos_end_overlapA){

                    cerr << "GFA_Parser::write_edge(): Vertex A overlap start position greater than or equal to vertex A overlap end position" << endl;
                    close();
                    return false;
                }

                if (pos_start_overlapB >= pos_end_overlapB){

                    cerr << "GFA_Parser::write_edge(): Vertex B overlap start position greater than or equal to vertex B overlap end position" << endl;
                    close();
                    return false;
                }

                if (v_gfa == 1){

                    if ((pos_end_overlapB - pos_start_overlapB) != (pos_end_overlapA - pos_start_overlapA)){

                        cerr << "GFA_Parser::write_edge(): Overlap lengths must be the same for vertex A and B in GFA format version 1" << endl;
                        close();
                        return false;
                    }

                    graph_out << "L" << "\t" <<
                    vertexA_id << "\t" << (strand_overlapA ? "+" : "-") << "\t" <<
                    vertexB_id << "\t" << (strand_overlapB ? "+" : "-") << "\t" <<
                    (pos_end_overlapA - pos_start_overlapA) << "M\n";
                }
                else {

                    graph_out << "E" << "\t";

                    if (edge_id == -1) graph_out << "*";
                    else graph_out << edge_id;

                    graph_out << "\t" <<
                    vertexA_id << (strand_overlapA ? "+" : "-") << "\t" <<
                    vertexB_id << (strand_overlapB ? "+" : "-") << "\t" <<
                    pos_start_overlapA << "\t" << pos_end_overlapA << "\t" <<
                    pos_start_overlapB << "\t" << pos_end_overlapB << "\t" <<
                    "*" << "\n";
                }
            }
            else {

                cerr << "GFA_Parser::write_edge(): Input file is not open in writing mode" << endl;
                return false;
            }

            return true;
        }

        const GFA_line read() {

            if (file_open_read){

                vector<string> line_fields;

                graph_in >> buffer;

                if (buffer[0] == 'S'){ // Segment line

                    stringstream ss(&buffer[2]); // Skip the first 2 char. of the line "S\t"
                    string substr;

                    while (ss.good()){ // Split line based on tabulation

                        getline(ss, substr, '\t');
                        line_fields.push_back(substr);
                    }

                    const size_t line_fields_sz = line_fields.size();
                    size_t i = 0;

                    s = Sequence();

                    if (v_gfa == 1){ // GFA format version 1

                        for (i = 0; i != line_fields_sz; ++i){

                            if (i == 0) s.id = line_fields[i];
                            else if (i == 1) s.seq = line_fields[i];
                            /*else {

                                const string sub = (*it).substr(0, 2);

                                if (sub == "LN") sscanf(&(sub.c_str()[4]), "%zu", &(s.len)); // Sequence length
                                else if ((sub != "RC") && (sub != "FC") && (sub != "KC") && (sub != "SH") && (sub != "UR") && (sub.substr(2, 2) == ":i")){
                                    // We assume any made-up tag that is an integer type is the coverage info
                                    s.coverage = atoi(&(sub.c_str()[4]));
                                }
                            }*/
                        }

                        if (i < 2){

                            cerr << "GFA_Parser::read(): Missing fields in Segment line" << endl;
                            close();
                        }
                    }
                    else {

                        for (i = 0; i != line_fields_sz; ++i){

                            if (i == 0) s.id = line_fields[i];
                            else if (i == 1) s.len = sscanf(line_fields[i].c_str(), "%zu", &(s.len));
                            else if (i == 2) s.seq = line_fields[i];
                            /*else if ((sub != "LN") && (sub != "RC") && (sub != "FC") && (sub != "KC") && (sub != "SH") && (sub != "UR") && (sub.substr(2, 2) == ":i")){
                                // We assume any made-up tag that is an integer type is the coverage info
                                s.coverage = atoi(&(sub.c_str()[4]));
                            }*/
                        }

                        if (i < 3){

                            cerr << "GFA_Parser::read(): Missing fields in Segment line" << endl;
                            close();
                        }
                    }

                    return GFA_line(&s, nullptr);
                }
                else if ((v_gfa == 1) && (buffer[0] == 'L')){ // Link line, only GFA v1

                    stringstream ss(&buffer[2]); // Skip the first 2 char. of the line "L\t"
                    string substr;

                    while (ss.good()){ // Split line based on tabulation

                        getline(ss, substr, '\t');
                        line_fields.push_back(substr);
                    }

                    const size_t line_fields_sz = line_fields.size();
                    size_t i = 0;

                    e = Edge();

                    for (i = 0; i != line_fields_sz; ++i){

                        if (i == 0) e.vertexA_id = line_fields[i];
                        else if (i == 1){

                            if (line_fields[i] == "+") e.strand_overlapA = true;
                            else if (line_fields[i] == "-") e.strand_overlapA = false;
                            else {

                                cerr << "GFA_Parser::read(): Orientation of Segment A on Link line is not + or -" << endl;
                                close();
                            }
                        }
                        else if (i == 2) e.vertexB_id = line_fields[i];
                        else if (i == 3){

                            if (line_fields[i] == "+") e.strand_overlapB = true;
                            else if (line_fields[i] == "-") e.strand_overlapB = false;
                            else {

                                cerr << "GFA_Parser::read(): Orientation of Segment B on Link line is not + or -" << endl;
                                close();
                            }
                        }
                        // We don't consider the rest (yet)
                    }

                    if (i < 4){

                        cerr << "GFA_Parser::read(): Missing fields in Link line" << endl;
                        close();
                    }

                    return GFA_line(nullptr, &e);
                }
                else if ((v_gfa == 2) && (buffer[0] == 'E')){ // Link line, only GFA v1
                }
            }

            return GFA_line();
        }

    private:

        string graph_filename;

        ofstream graphfile;

        ostream graph_out;
        istream graph_in;

        size_t v_gfa;
        size_t buff_sz;

        char* buffer;

        bool file_open_write;
        bool file_open_read;

        Sequence s;
        Edge e;
};

#endif
