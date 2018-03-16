#include <iostream>
#include <stack>
#include <queue>
#include <unordered_set>

#include <bifrost/ColoredCDBG.hpp>

using namespace std;

// Boolean class indicating if its associated unitig was traversed (b=true) or not (b=false)
// The class inherits from CCDBG_Data_t<myBool> to be used with ColoredCDBG and from CDBG_Data_t<myBool>
// to be used with CompactedDBG

class myBool : public CCDBG_Data_t<myBool>, CDBG_Data_t<myBool> {

    public:

        myBool() : b(NOT_VISITED_SEEN) {} // Initiate the boolean to "not visited"

        // Join method for ColoredCDBG
        static void join(const UnitigColorMap<myBool>& um_dest, const UnitigColorMap<myBool>& um_src){

            // When joining the unitig matching um_src to the unitig matching um_dest,
            // we set um_dest to "not visited" because it will be a new unitig in the graph.

            DataAccessor<myBool>* da = um_dest.getData(); // Get DataAccessor from unitig matching um_dest
            myBool* data = da->getData(um_dest); // Get boolean from DataAccessor

            data->set_not_seen_visited(); // Set the unitig to "not visited"
        }

        // Sub method for ColoredCDBG
        static void sub(const UnitigColorMap<myBool>& um_src, myBool* new_data, bool last_extraction){

            // This function creates a new unitig which is a sub-unitig from um_src
            // The new unitig created is set to "not visited" as a measure of precaution
            // (it is already initialed by default to "not visited" in the constructor)

            new_data->set_not_seen_visited();
        }

        // Join method for CompactedDBG
        static void join(const UnitigMap<myBool>& um_dest, const UnitigMap<myBool>& um_src){

            // When joining the unitig matching um_src to the unitig matching um_dest,
            // we set um_dest to "not visited" because it will be a new unitig in the graph.

            myBool* data = um_dest.getData(); // Get boolean directly from unitig matching um_dest

            data->set_not_seen_visited(); // Set the unitig to "not visited"
        }

        // Sub method for CompactedDBG
        static void sub(const UnitigMap<myBool>& um_src, myBool* new_data, bool last_extraction){

            // This function creates a new unitig which is a sub-unitig from um_src
            // The new unitig created is set to "not visited" as a measure of precaution
            // (it is already initialed by default to "not visited" in the constructor)

            new_data->set_not_seen_visited();
        }

        void toString() const {

            cout << "Unitig visited = " << (is_visited() ? "true" : "false") << endl;
            cout << "Unitig seen = " << (is_seen() ? "true" : "false") << endl;
        }

        inline void set_visited() { b = VISITED; } // Set the boolean to "visited"
        inline void set_seen() { b = SEEN; } // Set the boolean to "seen"
        inline void set_not_seen_visited() { b = NOT_VISITED_SEEN; } // Set the boolean to "not seen and not visited"

        inline bool is_visited() const { return (b == VISITED); } // return if the boolean is "visited"
        inline bool is_not_visited() const { return (b != VISITED); } // return if the boolean is "not visited"

        inline bool is_seen() const { return (b == SEEN); } // return if the boolean is "seen"
        inline bool is_not_seen() const { return (b != SEEN); } // return if the boolean is "not seen"

    private:

        const static uint8_t NOT_VISITED_SEEN = 0x0;
        const static uint8_t VISITED = 0x1;
        const static uint8_t SEEN = 0x2;

        uint8_t b;
};

void cleanMarking(const unordered_set<UnitigColorMap<myBool>, UnitigMapHash<DataAccessor<myBool>, DataStorage<myBool>, false>>& set_km_seen){

    unordered_set<UnitigColorMap<myBool>, UnitigMapHash<DataAccessor<myBool>, DataStorage<myBool>, false>>::const_iterator it, it_end;

    for (it = set_km_seen.begin(), it_end = set_km_seen.end(); it != it_end; ++it){

        const UnitigColorMap<myBool> ucm = *it;

        DataAccessor<myBool>* da_ucm = ucm.getData();
        myBool* data_ucm = da_ucm->getData(ucm);

        data_ucm->set_not_seen_visited();
    }
}

void cleanMarking(ColoredCDBG<myBool>& ccdbg){

    for (auto& unitig : ccdbg){

        DataAccessor<myBool>* da_ucm = unitig.getData();
        myBool* data_ucm = da_ucm->getData(unitig);

        data_ucm->set_not_seen_visited();
    }
}

void BFS_Recursive(const UnitigColorMap<myBool>& ucm){

    for (auto& successor : ucm.getSuccessors()){ // Iterate over successors of a unitig

        DataAccessor<myBool>* da = successor.getData(); // Get DataAccessor from unitig successor
        myBool* data = da->getData(successor); // Get boolean from DataAccessor

        if (data->is_not_visited()){ // If boolean indicates the unitig was not visited

            data->set_visited(); // Set boolean to indicate unitig was visited
        }
    }

    for (auto& predecessor : ucm.getPredecessors()){ // Iterate over predecessors of a unitig

        DataAccessor<myBool>* da = predecessor.getData(); // Get DataAccessor from unitig predecessor
        myBool* data = da->getData(predecessor); // Get boolean from DataAccessor

        if (data->is_not_visited()){ // If boolean indicates the unitig was not visited

            data->set_visited(); // Set boolean to indicate unitig was visited
        }
    }

    // Traverse successors
    for (auto& successor : ucm.getSuccessors()) BFS_Recursive(successor);
    // Traverse predecessors
    for (auto& predecessor : ucm.getPredecessors()) BFS_Recursive(predecessor);
}

void BFS_Iterative(const UnitigColorMap<myBool>& ucm){

    queue<UnitigColorMap<myBool>> q; // Create queue of unitig to traverse
    UnitigColorMap<myBool> ucm_tmp(ucm); // Create a non-const local copy of unitig given in parameter

    DataAccessor<myBool>* da = ucm_tmp.getData(); // Get DataAccessor from unitig
    myBool* data = da->getData(ucm_tmp); // Get boolean from DataAccessor

    data->set_visited(); // Set boolean to indicate unitig was visited

    q.push(ucm_tmp); // Push unitig to traverse on the stack

    while (!q.empty()){ // While they are unitigs to traverse in the stack

        ucm_tmp = q.front(); // Get unitig at the front of the queue

        q.pop(); // Delete unitig on the top of the stock from the stack

        for (auto& successor : ucm_tmp.getSuccessors()){ // Traverse successors

            DataAccessor<myBool>* da_succ = successor.getData(); // Get DataAccessor from successor
            myBool* data_succ = da_succ->getData(successor); // Get boolean from DataAccessor

            if (data_succ->is_not_visited()){ // If boolean indicates the successor was not visited

                data_succ->set_visited(); // Set boolean to indicate successor was visited

                q.push(successor); // Traverse neighbors of successor
            }
        }

        // Traverse predecessors
        for (auto& predecessor : ucm_tmp.getPredecessors()){

            DataAccessor<myBool>* da_pred = predecessor.getData(); // Get DataAccessor from predecessor
            myBool* data_pred = da_pred->getData(predecessor); // Get boolean from DataAccessor

            if (data_pred->is_not_visited()){ // If boolean indicates the predecessor was not visited

                data_pred->set_visited(); // Set boolean to indicate predecessor was visited

                q.push(predecessor); // Traverse neighbors of predecessor
            }
        }
    }
}

void DFS_Recursive(const UnitigColorMap<myBool>& ucm){

    for (auto& successor : ucm.getSuccessors()){ // Iterate over successors of a unitig

        DataAccessor<myBool>* da = successor.getData(); // Get DataAccessor from unitig successor
        myBool* data = da->getData(successor); // Get boolean from DataAccessor

        if (data->is_not_visited()){ // If boolean indicates the unitig was not visited

            data->set_visited(); // Set boolean to indicate unitig was visited

            DFS_Recursive(successor); // Traverse neighbors of successor
        }
    }

    for (auto& predecessor : ucm.getPredecessors()){ // Iterate over predecessors of a unitig

        DataAccessor<myBool>* da = predecessor.getData(); // Get DataAccessor from unitig predecessor
        myBool* data = da->getData(predecessor); // Get boolean from DataAccessor

        if (data->is_not_visited()){ // If boolean indicates the unitig was not visited

            data->set_visited(); // Set boolean to indicate unitig was visited

            DFS_Recursive(predecessor); // Traverse neighbors of predecessor
        }
    }
}

void DFS_Iterative(const UnitigColorMap<myBool>& ucm){

    stack<UnitigColorMap<myBool>> stck; // Create stack of unitig to traverse
    UnitigColorMap<myBool> ucm_tmp(ucm); // Create a non-const local copy of unitig given in parameter

    stck.push(ucm_tmp); // Push first unitig to traverse on the stack

    while (!stck.empty()){ // While they are unitigs to traverse in the stack

        ucm_tmp = stck.top(); // Get the unitig on top of the stack

        stck.pop(); // Delete unitig on the top of the stock from the stack

        DataAccessor<myBool>* da = ucm_tmp.getData(); // Get DataAccessor from unitig
        myBool* data = da->getData(ucm_tmp); // Get boolean from DataAccessor

        if (data->is_not_visited()){ // If boolean indicates the unitig was not visited

            data->set_visited(); // Set boolean to indicate unitig was visited

            // Add successors to stack of unitigs to traverse
            for (auto& successor : ucm_tmp.getSuccessors()) stck.push(successor);
            // Add predecessors to stack of unitigs to traverse
            for (auto& predecessor : ucm_tmp.getPredecessors()) stck.push(predecessor);
        }
    }
}

// - Parameter "iterative" indicates if you want to traverse the graph in a iterative (true)
// or recursive (false) way. Recursive version is limited to small components because of
// the stack size (recursive calls are piling up variables on your stack) which is not
// the case of the iterative version. Default is iterative.
void Traverse(ColoredCDBG<myBool>& ccdbg, const bool DFS = true, const bool iterative = true){

    for (auto& unitig : ccdbg){ // Iterate over unitigs of a colored de Bruijn graph

        DataAccessor<myBool>* da = unitig.getData(); // Get DataAccessor from unitig
        myBool* data = da->getData(unitig); // Get boolean from DataAccessor

        if (data->is_not_visited()){ // If boolean indicates the unitig was not visited

            if (iterative) DFS ? DFS_Iterative(unitig) : BFS_Iterative(unitig); // Traverse neighbors of unitig in an iterative manner
            else {

                data->set_visited(); // Set boolean to indicate unitig was visited

                DFS ? DFS_Recursive(unitig) : BFS_Recursive(unitig); // Traverse neighbors of unitig in a recursive manner
            }
        }
    }

    cleanMarking(ccdbg);
}

size_t getNbConnectedComponent(ColoredCDBG<myBool>& ccdbg){

    size_t nb_cc = 0; // Number of connected components

    for (auto& unitig : ccdbg){ // Iterate over unitigs of a colored de Bruijn graph

        DataAccessor<myBool>* da = unitig.getData(); // Get DataAccessor from unitig
        myBool* data = da->getData(unitig); // Get boolean from DataAccessor

        if (data->is_not_visited()){ // If boolean indicates the unitig was not visited

            ++nb_cc; // It's a new connected components

            DFS_Iterative(unitig); // Traverse neighbors of unitig in DFS
        }
    }

    cleanMarking(ccdbg);

    return nb_cc;
}

void getCoreGraph(const ColoredCDBG<myBool>& cdbg_in, CompactedDBG<myBool>& cdbg_out){

    size_t nb_colors = cdbg_in.getNbColors(); // Number of colors used in the graph

    for (const auto& unitig : cdbg_in){ // Iterate over unitigs of the colored de Bruijn graph

        // For each k-mer in the unitig
        for (size_t i = 0; i <= unitig.size - static_cast<size_t>(cdbg_in.getK()); ++i){

            const const_UnitigColorMap<myBool> ucm = unitig.getKmerMapping(i); // Get the mapping for this k-mer
            const UnitigColors<myBool>* uc = ucm.getData()->getUnitigColors(ucm); // Get the color set associated with this unitig

            bool isCore = true; // Is the k-mer a core k-mer (has all colors, default is yes)?

            for (size_t color_id = 0; (color_id != nb_colors) && isCore; ++color_id) isCore = uc->contains(ucm, color_id);

            // If k-mer is core, add it to a new compacted de Bruijn graph
            if (isCore) cdbg_out.add(unitig.getUnitigKmer(i).toString());
        }
    }
}

pair<UnitigColorMap<myBool>, UnitigColorMap<myBool>> extractSuperBubble(ColoredCDBG<myBool>& ccdbg, const UnitigColorMap<myBool>& s){

    vector<UnitigColorMap<myBool>> vertices_visit;

    pair<UnitigColorMap<myBool>, UnitigColorMap<myBool>> p;

    unordered_set<UnitigColorMap<myBool>, UnitigMapHash<DataAccessor<myBool>, DataStorage<myBool>, false>> set_km_seen;

    UnitigColorMap<myBool> v(s);

    vertices_visit.push_back(v);

    while (!vertices_visit.empty()){

        v = vertices_visit.back();

        vertices_visit.pop_back();

        DataAccessor<myBool>* da = v.getData(); // Get DataAccessor from unitig
        myBool* data = da->getData(v); // Get boolean from DataAccessor

        data->set_visited(); // Set boolean to indicate unitig was visited
        set_km_seen.insert(v);

        size_t nb_succ = 0;

        for (const auto& u : v.getSuccessors()) ++nb_succ;

        if (nb_succ == 0){ // Unitig is a tip

            cleanMarking(set_km_seen);
            return p;
        }

        for (auto& u : v.getSuccessors()){

            if (u == s){ // Cycle

                cleanMarking(set_km_seen);
                return p;
            }

            DataAccessor<myBool>* da_u = u.getData(); // Get DataAccessor from successor
            myBool* bool_u = da_u->getData(u); // Get boolean from DataAccessor

            bool_u->set_seen(); // Set boolean to indicate unitig was seen
            set_km_seen.insert(u);

            bool all_predecessor_visited = true;

            for (const auto& predecessor : u.getPredecessors()){

                const DataAccessor<myBool>* da_pred = predecessor.getData();
                const myBool* data_pred = da_pred->getData(predecessor);

                if (data_pred->is_not_visited()){

                    all_predecessor_visited = false;
                    break;
                }
            }

            if (all_predecessor_visited) vertices_visit.push_back(u);
        }

        if (vertices_visit.size() == 1){

            bool not_seen = true;

            unordered_set<UnitigColorMap<myBool>, UnitigMapHash<DataAccessor<myBool>, DataStorage<myBool>, false>>::const_iterator it, it_end;

            for (it = set_km_seen.begin(), it_end = set_km_seen.end(); it != it_end; ++it){

                const const_UnitigColorMap<myBool> cucm = *it;

                if (cucm != vertices_visit[0]){

                    const DataAccessor<myBool>* da_cucm = cucm.getData();
                    const myBool* data_cucm = da_cucm->getData(cucm);

                    if (data_cucm->is_seen()){

                        not_seen = false;
                        break;
                    }
                }
            }

            if (not_seen){

                for (const auto& successor : vertices_visit[0].getSuccessors()){

                    if (successor == s){ // cycle

                        cleanMarking(set_km_seen);
                        return p;
                    }
                }

                p.first = s;
                p.second = vertices_visit[0];

                cleanMarking(set_km_seen);
                return p;
            }
        }
    }

    cleanMarking(set_km_seen);
    return p;
}

int main(int argc, char *argv[])
{

    if (argc != 1){

        ColoredCDBG<myBool> ccdbg(31); // Create a new (empty) colored dBG of 31-mers
        CCDBG_Build_opt opt;

        opt.nb_threads = 4; // Use 4 threads when possible
        opt.verbose = true; // Print messages during execution

        // Read input filenames
        for (int i = 1; i != argc; ++i) opt.filename_seq_in.push_back(string(argv[i]));

        // Check if input files can be opened
        for (vector<string>::const_iterator it = opt.filename_seq_in.begin(); it != opt.filename_seq_in.end(); ++it){

            FILE* fp = fopen(it->c_str(), "r");

            if (fp == NULL) {

                cerr << "Could not open input FASTA/FASTQ file " << *it << endl;
                exit(1);
            }
            else fclose(fp);
        }

        cout << "=== Building graph ===" << endl;

        ccdbg.build(opt);

        cout << "=== Mapping colors to unitigs ===" << endl;

        ccdbg.mapColors(opt);

        cout << "=== Traversing graph: Depth First Search ===" << endl;

        Traverse(ccdbg, true);

        cout << "=== Traversing graph: Breadth First Search ===" << endl;

        Traverse(ccdbg, false);

        cout << "=== Computing number of connected components ===" << endl;

        const size_t nb_cc = getNbConnectedComponent(ccdbg);

        cout << nb_cc << " connected components found" << endl;

        cout << "=== Computing super bubbles ===" << endl;

        size_t nb_super_bubble = 0;
        size_t nb_unitig_processed = 0;

        for (const auto& unitig : ccdbg){

            const pair<UnitigColorMap<myBool>, UnitigColorMap<myBool>> p = extractSuperBubble(ccdbg, unitig);

            if (!p.first.isEmpty && !p.second.isEmpty) ++nb_super_bubble;

            ++nb_unitig_processed;

            if (nb_unitig_processed % 100000 == 0){

                cout << "Processed " << nb_unitig_processed << " unitigs: " << nb_super_bubble << " super bubbles so far" << endl;
            }
        }

        cout << nb_super_bubble << " super bubbles found ===" << endl;

        cout << "=== Extracting core graph ===" << endl;

        CompactedDBG<myBool> cdbg(31);

        getCoreGraph(ccdbg, cdbg);
    }
    else cerr << "No input FASTA/FASTQ file(s) provided" << endl;
}
