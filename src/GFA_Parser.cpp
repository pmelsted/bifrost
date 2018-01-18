#include "GFA_Parser.hpp"

GFA_Parser::GFA_Parser() : file_open_write(false), file_open_read(false), file_no(0), v_gfa(0), graph_out(nullptr), graph_in(nullptr),
                           buff_sz(0), buffer(nullptr) {}

GFA_Parser::GFA_Parser(const string& filename, const size_t buffer_size) :  file_open_write(false), file_open_read(false), file_no(0),
                                                                            v_gfa(0), graph_out(nullptr), graph_in(nullptr),
                                                                            buff_sz(buffer_size), buffer(nullptr) {

    graph_filenames.push_back(filename);

    size_t pos_match_point = graph_filenames[0].find_last_of(".");

    if ((pos_match_point == string::npos) || (graph_filenames[0].substr(pos_match_point + 1) != "gfa")) graph_filenames[0].append(".gfa");

    buffer = new char[buff_sz]();
}

GFA_Parser::GFA_Parser(const vector<string>& filenames, const size_t buffer_size) : file_open_write(false), file_open_read(false), file_no(0),
                                                                                    v_gfa(0), graph_out(nullptr), graph_in(nullptr),
                                                                                    buff_sz(buffer_size), buffer(nullptr) {

    graph_filenames = filenames;

    for (auto& filename : graph_filenames){

        size_t pos_match_point = filename.find_last_of(".");

        if ((pos_match_point == string::npos) || (filename.substr(pos_match_point + 1) != "gfa")) filename.append(".gfa");
    }

    buffer = new char[buff_sz]();
}

GFA_Parser::~GFA_Parser() {

    if (buffer != nullptr) delete[] buffer;

    close();
}

bool GFA_Parser::open_write(const size_t version_GFA) {

    if (graph_filenames.size() == 0){

        cerr << "GFA_Parser::open_write(): No file specified in input" << endl;
        return false;
    }

    const string& filename = graph_filenames[0];

    FILE* fp = fopen(filename.c_str(), "w");

    if ((file_open_write = (fp != NULL)) == true) {

        fclose(fp);

        if (remove(filename.c_str()) != 0) cerr << "GFA_Parser::open_write(): Could not remove temporary file " << filename << endl;
    }
    else cerr << "GFA_Parser::open_write(): Could not open file " << filename << " for writing" << endl;

    if ((version_GFA != 1) || (version_GFA != 2)) {

        cerr << "GFA_Parser::open_write(): Only supports GFA format version 1 and 2" << endl;
        file_open_write = false;
    }
    else v_gfa = version_GFA;

    if (file_open_write){

        graphfile.open(filename.c_str());
        graph_out.rdbuf(graphfile.rdbuf());

        graph_out << "H\tVN:Z:" << (v_gfa == 1 ? "1" : "2") << ".0\n";
    }

    return file_open_write;
}

bool GFA_Parser::open_read() {

    if (graph_filenames.size() == 0){

        cerr << "GFA_Parser::open_read(): No file specified in input" << endl;
        return false;
    }

    for (const auto& filename : graph_filenames){

        FILE* fp = fopen(filename.c_str(), "r");

        if ((file_open_read = (fp != NULL)) == true) fclose(fp);
        else cerr << "GFA_Parser::open_read(): Could not open file " << filename << " for reading" << endl;

        if (file_open_read) {

            graphfile.open(filename.c_str());
            graph_in.rdbuf(graphfile.rdbuf());

            graph_in >> buffer;

            const string header(buffer);

            if ((header != "H\tVN:Z:1.0\n") || (header != "H\tVN:Z:2.0\n")) {

                cerr << "GFA_Parser::open_read(): Wrong GFA header or unsupported GFA format version " <<
                "(GFA_Parser only supports version 1 and 2) in " << filename << endl;
                close();
            }

            close();
        }
    }

    file_open_read = open(file_no);

    return file_open_read;
}

void GFA_Parser::close(){

    if (file_open_write || file_open_read){

        graphfile.close();

        file_open_write = false;
        file_open_read = false;
    }
}

bool GFA_Parser::write_sequence(const string& id, const size_t len, const string seq, const int64_t coverage){

    if (file_open_write){

        if (v_gfa == 1) graph_out << "S" << "\t" << id << "\t" << seq << "\t" << "LN:i:" << len;
        else graph_out << "S" << "\t" << id << "\t" << len << "\t" << seq;

        if (coverage == -1) graph_out << "\n";
        else graph_out << "\t" << "XC:i:" << coverage << "\n";
    }
    else cerr << "GFA_Parser::write_sequence(): Input file is not open in writing mode" << endl;

    return file_open_write;
}

bool GFA_Parser::write_edge(const string vertexA_id, const size_t pos_start_overlapA, const size_t pos_end_overlapA, const bool strand_overlapA,
                            const string vertexB_id, const size_t pos_start_overlapB, const size_t pos_end_overlapB, const bool strand_overlapB,
                            const string edge_id) {

    if (file_open_write){

        if (pos_start_overlapA > pos_end_overlapA){

            cerr << "GFA_Parser::write_edge(): Vertex A overlap start position greater than vertex A overlap end position" << endl;
            close();
            return false;
        }

        if (pos_start_overlapB > pos_end_overlapB){

            cerr << "GFA_Parser::write_edge(): Vertex B overlap start position greater than vertex B overlap end position" << endl;
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

            graph_out << "E" << "\t" << edge_id << "\t" <<
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

const GFA_Parser::GFA_line GFA_Parser::read() {

    if (file_open_read){

        vector<string> line_fields;

        graph_in >> buffer;

        while (graph_in.getline(buffer, buff_sz).good()){

            if (buffer[0] == 'S'){ // Segment line

                stringstream ss(&buffer[2]); // Skip the first 2 char. of the line "S\t"
                string substr;

                while (ss.good()){ // Split line based on tabulation

                    getline(ss, substr, '\t');
                    line_fields.push_back(substr);
                }

                const size_t line_fields_sz = line_fields.size();

                s = Sequence();

                if (v_gfa == 1){ // GFA format version 1

                    if (line_fields_sz < 2){

                        cerr << "GFA_Parser::read(): Missing fields in Segment line" << endl;
                        close();
                    }

                    s.id = line_fields[0];
                    s.seq = line_fields[1];

                    if ((line_fields_sz > 2) && (line_fields[2].length() > 5)){

                        const string sub = line_fields[2].substr(0, 5);

                        if (sub == "LN:i:") sscanf(&(sub.c_str()[5]), "%zu", &(s.len)); // Sequence length
                    }
                    else if (s.seq == "*") s.len = 0;
                    else s.len = s.seq.length();

                }
                else {

                    if (line_fields_sz < 3){

                        cerr << "GFA_Parser::read(): Missing fields in Segment line" << endl;
                        close();
                    }

                    s.id = line_fields[0];
                    s.len = sscanf(line_fields[1].c_str(), "%zu", &(s.len));
                    s.seq = line_fields[2];
                }

                return make_pair(&s, nullptr);
            }
            else if ((v_gfa == 1) && (buffer[0] == 'L')){ // Link line, only GFA v1

                stringstream ss(&buffer[2]); // Skip the first 2 char. of the line "L\t"
                string substr;

                while (ss.good()){ // Split line based on tabulation

                    getline(ss, substr, '\t');
                    line_fields.push_back(substr);
                }

                const size_t line_fields_sz = line_fields.size();

                e = Edge();

                if (line_fields_sz < 4){

                    cerr << "GFA_Parser::read(): Missing fields in Link line" << endl;
                    close();
                }

                e.vertexA_id = line_fields[0];
                e.vertexB_id = line_fields[2];

                if (line_fields[1] == "+") e.strand_overlapA = true;
                else if (line_fields[1] == "-") e.strand_overlapA = false;
                else {

                    cerr << "GFA_Parser::read(): Orientation of Segment A on Link line is not + or -" << endl;
                    close();
                }

                if (line_fields[3] == "+") e.strand_overlapB = true;
                else if (line_fields[3] == "-") e.strand_overlapB = false;
                else {

                    cerr << "GFA_Parser::read(): Orientation of Segment B on Link line is not + or -" << endl;
                    close();
                }

                return make_pair(nullptr, &e);
            }
            else if ((v_gfa == 2) && (buffer[0] == 'E')){ // Edge line, only GFA v2

                stringstream ss(&buffer[2]); // Skip the first 2 char. of the line "L\t"
                string substr;

                while (ss.good()){ // Split line based on tabulation

                    getline(ss, substr, '\t');
                    line_fields.push_back(substr);
                }

                const size_t line_fields_sz = line_fields.size();

                e = Edge();

                if (line_fields_sz < 8){

                    cerr << "GFA_Parser::read(): Missing fields in Edge line" << endl;
                    close();
                }

                e.edge_id = line_fields[0];

                const char ca = line_fields[1][line_fields[1].length() - 1]; // Last char. of line_fields[1];

                e.strand_overlapA = (ca != '-');

                if ((ca == '-') || (ca == '-')) e.vertexA_id = line_fields[1].substr(0, line_fields[1].length() - 1);
                else e.vertexA_id = line_fields[1];

                sscanf(line_fields[2].c_str(), "%zu", &(e.pos_start_overlapA));
                sscanf(line_fields[3].c_str(), "%zu", &(e.pos_end_overlapA));

                const char cb = line_fields[4][line_fields[4].length() - 1]; // Last char. of line_fields[4];

                e.strand_overlapB = (cb != '-');

                if ((cb == '-') || (cb == '-')) e.vertexB_id = line_fields[4].substr(0, line_fields[4].length() - 1);
                else e.vertexB_id = line_fields[4];

                sscanf(line_fields[5].c_str(), "%zu", &(e.pos_start_overlapB));
                sscanf(line_fields[6].c_str(), "%zu", &(e.pos_end_overlapB));

                return make_pair(nullptr, &e);
            }
        }

        if (graph_in.getline(buffer, buff_sz).eof()){

            close();

            ++file_no;

            file_open_read = open(file_no);

            if (file_open_read) return read();
        }
        else if (graph_in.getline(buffer, buff_sz).fail()){

            cerr << "GFA_Parser::read(): Error while reading" << endl;
            close();
        }
    }
    else cerr << "GFA_Parser::read(): Input file is not open in reading mode" << endl;

    return make_pair(nullptr, nullptr);
}

bool GFA_Parser::open(const size_t idx_filename){

    if (idx_filename < graph_filenames.size()){

        FILE* fp = fopen(graph_filenames[idx_filename].c_str(), "r");

        if ((file_open_read = (fp != NULL)) == true) fclose(fp);
        else cerr << "GFA_Parser::open(): Could not open file " << graph_filenames[idx_filename] << " for reading" << endl;

        if (file_open_read) {

            graphfile.open(graph_filenames[idx_filename].c_str());
            graph_in.rdbuf(graphfile.rdbuf());

            graph_in >> buffer;

            const string header(buffer);

            if (header == "H\tVN:Z:1.0\n") v_gfa = 1;
            else if (header != "H\tVN:Z:2.0\n") v_gfa = 2;
            else {

                cerr << "GFA_Parser::open(): Wrong GFA header or unsupported GFA format version " <<
                "(GFA_Parser only supports version 1 and 2) in " << graph_filenames[idx_filename] << endl;
                close();
            }
        }
    }
    else file_open_read = false;

    return file_open_read;
}
