#ifndef COMPACTED_DBG_HPP
#define COMPACTED_DBG_HPP

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cstdio>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <vector>

#include <thread>
#include <atomic>
#include <mutex>

//#include <jemalloc/jemalloc.h>

#include "BlockedBloomFilter.hpp"
#include "Common.hpp"
#include "File_Parser.hpp"
#include "FASTX_Parser.hpp"
#include "GFA_Parser.hpp"
#include "Kmer.hpp"
#include "KmerHashTable.hpp"
#include "KmerIterator.hpp"
#include "KmerStream.hpp"
#include "minHashIterator.hpp"
#include "RepHash.hpp"
#include "TinyVector.hpp"
#include "Unitig.hpp"
#include "UnitigIterator.hpp"
#include "UnitigMap.hpp"

#define MASK_CONTIG_ID (0xffffffff00000000)
#define MASK_CONTIG_TYPE (0x80000000)
#define MASK_CONTIG_POS (0x7fffffff)
#define RESERVED_ID (0xffffffff)

#define DEFAULT_K 31
#define DEFAULT_G 23

/** @file src/CompactedDBG.hpp
* Interface for the Compacted de Bruijn graph API.
* Code snippets using this interface are provided in snippets/test.cpp.
*/

using namespace std;

/** @struct CDBG_Build_opt
* @brief Members of this structure are parameters for CompactedDBG<T>::build, except for:
* - CDBG_Build_opt::k and CDBG_Build_opt::g as they are parameters of the graph constructor.
* - CDBG_Build_opt::clipTips, CDBG_Build_opt::deleteIsolated and CDBG_Build_opt::useMercyKmers are used
* by CompactedDBG<T>::simplify
* - CDBG_Build_opt::prefixFilenameOut and CDBG_Build_opt::outputGFA are used by CompactedDBG<T>::write
* Most parameters have default values.
* An example of using such a structure is shown in src/Bifrost.cpp.
* @var CDBG_Build_opt::reference_mode
* Input are assembled genomes or graph files (true), no filtering step performed, all k-mers are used.
* Otherwise, input are reads (false), filtering step is performed, k-mers with coverage < 2 are removed.
* Default is false.
* @var CDBG_Build_opt::verbose
* Print information messages during execution if true. Default is false.
* @var CDBG_Build_opt::clipTips
* Clip short tips (length < 2k) of the graph (not used by CompactedDBG<T>::build). Default is false.
* @var CDBG_Build_opt::deleteIsolated
* Remove short isolated unitigs (length < 2k) of the graph (not used by CompactedDBG<T>::build).
* Default is false.
* @var CDBG_Build_opt::useMercyKmers
* Keep in the graph low coverage k-mers (cov=1) connecting tips of the graph. Default is false.
* @var CDBG_Build_opt::k
* Length of k-mers (not used by CompactedDBG<T>::build). Default is 31.
* @var CDBG_Build_opt::g
* Length of g-mers, the minimizers, such that g < k (not used by CompactedDBG<T>::build).
* Default is 23.
* @var CDBG_Build_opt::nb_threads
* Number of threads to use for building the graph. Default is 1.
* @var CDBG_Build_opt::read_chunksize
* Number of reads shared and processed by CDBG_Build_opt::nb_threads threads at the same time.
* Default is 10000.
* @var CDBG_Build_opt::unitig_size
* Maximum length of a unitig. Default is 100000.
* @var CDBG_Build_opt::nb_unique_kmers
* Estimated number (upper bound) of different k-mers in the input FASTA/FASTQ/GFA files of
* CDBG_Build_opt::filename_in. Default is KmerStream estimation.
* @var CDBG_Build_opt::nb_non_unique_kmers
* Estimated number (upper bound) of different k-mers occurring twice or more in the FASTA/FASTQ/GFA files
* of CDBG_Build_opt::filename_in. Default is a KmerStream estimation.
* @var CDBG_Build_opt::nb_bits_unique_kmers_bf
* Number of Bloom filter bits per k-mer occurring at least once in the FASTA/FASTQ/GFA files of
* CDBG_Build_opt::filename_in. Default is 14.
* @var CDBG_Build_opt::nb_bits_non_unique_kmers_bf
* Number of Bloom filter bits per k-mer occurring at least twice in the FASTA/FASTQ/GFA files of
* CDBG_Build_opt::filename_in. Default is 14.
* @var CDBG_Build_opt::prefixFilenameOut
* Prefix for the name of the file to which the graph must be written. Mandatory parameter.
* @var CDBG_Build_opt::inFilenameBBF
* String containing the name of a Bloom filter file that is generated by CompactedDBG<T>::filter.
* If empty, CompactedDBG<T>::filter is called. Otherwise, the Bloom filter is loaded from this file
* and CompactedDBG<T>::filter is not called. Note that you need such a Bloom filter even in reference
* mode. Default is no input file.
* @var CDBG_Build_opt::outFilenameBBF
* String containing the name of a Bloom filter file that will be generated by CompactedDBG<T>::filter.
* If empty, the file is not created. Otherwise, the Bloom filter is written to this file. Default is
* no output file.
* @var CDBG_Build_opt::filename_in
* Vector of strings, each string is the name of a FASTA/FASTQ/GFA file to use for the graph construction.
* Mandatory parameter.
* @var CDBG_Build_opt::outputGFA
* Boolean indicating if the graph is written to a GFA file (true) or if the unitigs are written to a
* FASTA file (false). Default is true.
*/
struct CDBG_Build_opt {

    bool reference_mode;
    bool verbose;

    size_t nb_threads;
    size_t read_chunksize;
    size_t unitig_size;
    size_t nb_unique_kmers;
    size_t nb_non_unique_kmers;
    size_t nb_bits_unique_kmers_bf;
    size_t nb_bits_non_unique_kmers_bf;

    string inFilenameBBF;
    string outFilenameBBF;

    vector<string> filename_in;

    // The following members are not used by CompactedDBG<T>::build
    // but you can set them to use them as parameters for other functions
    // such as CompactedDBG<T>::simplify or CompactedDBG<T>::write.

    size_t k, g;

    bool clipTips;
    bool deleteIsolated;
    bool useMercyKmers;

    bool outputGFA;

    string prefixFilenameOut;

    CDBG_Build_opt() :  nb_threads(1), k(DEFAULT_K), g(DEFAULT_G), nb_unique_kmers(0), nb_non_unique_kmers(0), nb_bits_unique_kmers_bf(14),
                        nb_bits_non_unique_kmers_bf(14), read_chunksize(10000), unitig_size(1000000), reference_mode(false),
                        verbose(false), clipTips(false), deleteIsolated(false), useMercyKmers(false), outputGFA(true), inFilenameBBF(""),
                        outFilenameBBF("") {}
};

/** @typedef const_UnitigMap
* @brief const_UnitigMap is a constant UnitigMap. The main difference in its usage with a UnitigMap object
* is when you call the method UnitigMap::getCompactedDBG(): with a const_UnitigMap, this method returns
* a pointer to a constant CompactedDBG (you can't modify it).
*/
template<typename U, typename G> using const_UnitigMap = UnitigMap<U, G, true>;

/** @class CDBG_Data_t
* @brief If data are to be associated with the unitigs of the compacted de Bruijn graph, those data
* must be wrapped into a class that inherits from the abstract class CDBG_Data_t. Otherwise it will
* not compile. To associate data of type myUnitigData to unitigs, class myUnitigData must be declared as follows:
* \code{.cpp}
* class myUnitigData : public CDBG_Data_t<myUnitigData, myGraphData> { ... };
* ...
* CompactedCDBG<myUnitigData, myGraphData> cdbg;
* \endcode
* CDBG_Data_t has two template parameters: the type of unitig data and the type of graph data. Indeed, if class myUnitigData
* is going to be used in combination with class myGraphData for a CompactedDBG, myUnitigData must be "aware" of myGraphData
* for the parameters of its mandatory functions (see below). If no graph data is used, you don't have to specify it:
* \code{.cpp}
* class myUnitigData : public CDBG_Data_t<myUnitigData> { ... };
* class myUnitigData : public CDBG_Data_t<myUnitigData, void> { ... }; // Equivalent to previous notation
* ...
* CompactedCDBG<myUnitigData> cdbg;
* \endcode
* Because CDBG_Data_t is an abstract class, methods CDBG_Data_t::join, CDBG_Data_t::split and CDBG_Data_t::serialize
* must be implemented in your wrapper. IMPORTANT: If you don't overload those methods, default ones that have no effects
* will be applied!
* An example of using such a structure is shown in snippets/test.cpp.
*/
template<typename Unitig_data_t, typename Graph_data_t = void> //Curiously Recurring Template Pattern (CRTP)
class CDBG_Data_t {

    public:

        /** Join data of two unitigs (each represented with a UnitigMap given as parameter) which are going to be concatenated.
        * Specifically, if A is the unitig represented by parameter um_dest and B is the unitig represented by parameter um_src
        * then, after the call to this function, A will become the concatenation of itself with B (A = AB) and B will be removed.
        * Be careful that if um_dest.strand = false, then the reverse-complement of A is going to be used in the concatenation.
        * Reciprocally, if um_src.strand = false, then the reverse-complement of B is going to be used in the concatenation.
        * The data of each unitig can be accessed through the method UnitigMap::getData(). Note that this method is static.
        * @param um_dest is a UnitigMap object representing a unitig (and its data) to which another unitig is going to be appended.
        * @param um_src is a UnitigMap object representing a unitig (and its data) that will be appended at the end of the unitig
        * represented by parameter um_dest.
        */
        static void join(const UnitigMap<Unitig_data_t, Graph_data_t>& um_dest, const UnitigMap<Unitig_data_t, Graph_data_t>& um_src){}

        /** Extract data from a unitig A to be associated with a unitig B which is the unitig mapping given by the UnitigMap object
        * um_src. Hence, B = A[um_src.dist, um_src.dist + um_src.len + k - 1]. Be careful that if um_src.strand = false, then B will
        * be extracted from the reverse-complement of A, i.e, B = rev(A[um_src.dist, um_src.dist + um_src.len + k - 1]). Note that
        * this method is static.
        * @param um_src is a UnitigMap object representing the mapping to a unitig A from which a new unitig B will be created, i.e,
        * B = A[um_src.dist, um_src.dist + um_src.len + k - 1] or B = rev(A[um_src.dist, um_src.dist + um_src.len + k - 1]) if
        * um_src.strand == false.
        * @param new_data is a pointer to a newly constructed object that you can fill in with new data to associate with unitig B.
        * @param last_extraction is a boolean indicating if this is the last call to this function on the reference unitig used for the
        * mapping given by um_src. If last_extraction is true, the reference unitig of um_src will be removed from the graph right after
        * this function returns.
        */
        static void sub(const UnitigMap<Unitig_data_t, Graph_data_t>& um_src, Unitig_data_t* new_data, bool last_extraction){}

        /** Serialize the data to a string. This function is used when the graph is written to disk in GFA format.
        * If the returned string is not empty, the string is appended to an optional field of the Segment line matching the unitig
        * of this data. If the returned string is empty, no optional field or string are appended to the Segment line matching the
        * unitig of this data.
        */
        string serialize() const { return string(); }
};

/** @class CompactedDBG
* @brief Represent a Compacted de Bruijn graph. The two template parameters of this class corresponds to the type of data
* to associate with the unitigs of the graph (unitig data) and the type of data to associate with the graph (graph data).
* If no template parameters are specified or if the types are void, no data are associated with the unitigs nor the graph and
* no memory will be allocated for such data.
* \code{.cpp}
* CompactedDBG<> cdbg_1; // No unitig data, no graph data
* CompactedDBG<void> cdbg_2; // Equivalent to previous notation
* CompactedDBG<void, void> cdbg_3; // Equivalent to previous notation
* CompactedDBG<myUnitigData> cdbg_4; // An object of type myUnitigData will be associated with each unitig, no graph data
* CompactedDBG<myUnitigData, void> cdbg_5; // Equivalent to previous notation
* CompactedDBG<void, myGraphData> cdbg_6; // No unitig data, an object of type myGraphData will be associated with the graph
* CompactedDBG<myUnitigData, myGraphData> cdbg_7; // Unitig data of type myUnitigData for each unitig, graph data of type myGraphData
* \endcode
* If data are to be associated with the unitigs, these data must be wrapped into a class that inherits from the abstract class
* CDBG_Data_t, such as in:
* \code{.cpp}
* class myUnitigData : public CDBG_Data_t<myUnitigData> { ... };
* CompactedDBG<myUnitigData> cdbg;
* \endcode
*/
template<typename Unitig_data_t = void, typename Graph_data_t = void>
class CompactedDBG {

    static_assert(is_void<Unitig_data_t>::value || is_base_of<CDBG_Data_t<Unitig_data_t, Graph_data_t>, Unitig_data_t>::value,
                  "Type of data associated with vertices of class CompactedDBG must be void (no data) or a class extending class CDBG_Data_t");

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    public:

        template<typename U, typename G, bool is_const> friend class UnitigMap;
        template<typename U, typename G, bool is_const> friend class unitigIterator;
        template<typename U, typename G, bool is_const> friend class neighborIterator;

        typedef unitigIterator<U, G, false> iterator; /**< An iterator for the unitigs of the graph. No specific order is assumed. */
        typedef unitigIterator<U, G, true> const_iterator; /**< A constant iterator for the unitigs of the graph. No specific order is assumed. */

        /** Constructor (set up an empty compacted dBG).
        * @param kmer_length is the length k of k-mers used in the graph (each unitig is of length at least k).
        * @param minimizer_length is the length g of minimizers (g < k) used in the graph.
        */
        CompactedDBG(int kmer_length = DEFAULT_K, int minimizer_length = DEFAULT_G);

        /** Copy constructor (copy a compacted de Bruijn graph).
        * This function is expensive in terms of time and memory as the content of a compacted
        * de Bruijn graph is copied.  After the call to this function, the same graph exists twice in memory.
        * @param o is a constant reference to the compacted de Bruijn graph to copy.
        */
        CompactedDBG(const CompactedDBG& o); // Copy constructor

        /** Move constructor (move a compacted de Bruijn graph).
        * The content of o is moved ("transfered") to a new compacted de Bruijn graph.
        * The compacted de Bruijn graph referenced by o will be empty after the call to this constructor.
        * @param o is a reference on a reference to the compacted de Bruijn graph to move.
        */
        CompactedDBG(CompactedDBG&& o); // Move constructor

        /** Destructor.
        */
        virtual ~CompactedDBG();

        /** Copy assignment operator (copy a compacted de Bruijn graph).
        * This function is expensive in terms of time and memory as the content of a compacted
        * de Bruijn graph is copied.  After the call to this function, the same graph exists twice in memory.
        * @param o is a constant reference to the compacted de Bruijn graph to copy.
        * @return a reference to the compacted de Bruijn which is the copy.
        */
        CompactedDBG<U, G>& operator=(const CompactedDBG& o);

        /** Move assignment operator (move a compacted de Bruijn graph).
        * The content of o is moved ("transfered") to a new compacted de Bruijn graph.
        * The compacted de Bruijn graph referenced by o will be empty after the call to this operator.
        * @param o is a reference on a reference to the compacted de Bruijn graph to move.
        * @return a reference to the compacted de Bruijn which has (and owns) the content of o.
        */
        CompactedDBG<U, G>& operator=(CompactedDBG&& o);

        /** Clear the graph: empty the graph and reset its parameters.
        */
        void clear();

        /** Empty the graph (does not reset its parameters).
        */
        void empty();

        /** Build the Compacted de Bruijn graph.
        * @param opt is a structure from which the members are parameters of this function. See CDBG_Build_opt.
        * @return boolean indicating if the graph has been built successfully.
        */
        bool build(CDBG_Build_opt& opt);

        /** Simplify the Compacted de Bruijn graph: clip short (< 2k length) tips and/or delete short (< 2k length) isolated unitigs.
        * @param delete_short_isolated_unitigs is a boolean indicating short isolated unitigs must be removed.
        * @param clip_short_tips is a boolean indicating short tips must be clipped.
        * @param verbose is a boolean indicating if information messages must be printed during the function execution.
        * @return boolean indicating if the graph has been simplified successfully.
        */
        bool simplify(const bool delete_short_isolated_unitigs = true, const bool clip_short_tips = true, const bool verbose = false);

        /** Write the Compacted de Bruijn graph to disk (GFA1 format).
        * @param output_filename is a string containing the name of the file in which the graph will be written.
        * @param nb_threads is a number indicating how many threads can be used to write the graph to disk.
        * @param GFA_output indicates if the graph will be output in GFA format (true) or FASTA format (false).
        * @param verbose is a boolean indicating if information messages must be printed during the function execution.
        * @return boolean indicating if the graph has been written successfully.
        */
        bool write(const string output_filename, const size_t nb_threads = 1, const bool GFA_output = true, const bool verbose = false);

        /** Find the unitig containing the queried k-mer in the Compacted de Bruijn graph.
        * @param km is the queried k-mer (see Kmer class). It does not need to be a canonical k-mer.
        * @param extremities_only is a boolean indicating if the k-mer must be searched only in the unitig heads and tails (extremities_only = true).
        * By default, the k-mer is searched everywhere (extremities_only = false) but is is slightly slower than looking only in the unitig heads and tails.
        * @return UnitigMap<U, G> object containing the k-mer mapping information to the unitig containing the queried k-mer (if present).
        * If the queried k-mer is not found, UnitigMap::isEmpty = true (see UnitigMap class).
        */
        UnitigMap<U, G> find(const Kmer& km, const bool extremities_only = false);

        /** Find the unitig containing the queried k-mer in the Compacted de Bruijn graph.
        * @param km is the queried k-mer (see Kmer class). It does not need to be a canonical k-mer.
        * @param extremities_only is a boolean indicating if the k-mer must be searched only in the unitig heads and tails (extremities_only = true).
        * By default, the k-mer is searched everywhere (extremities_only = false) but is is slightly slower than looking only in the unitig heads and tails.
        * @return const_UnitigMap<U, G> object containing the k-mer mapping information to the unitig having the queried k-mer (if present).
        * If the k-mer is not found, const_UnitigMap::isEmpty = true (see UnitigMap class).
        */
        const_UnitigMap<U, G> find(const Kmer& km, const bool extremities_only = false) const;

        /** Add a sequence to the Compacted de Bruijn graph. Non-{A,C,G,T} characters such as Ns are discarded.
        * The function automatically breaks the sequence into unitig(s). Those unitigs can be stored as the reverse-complement
        * of the input sequence.
        * @param seq is a string containing the sequence to insert.
        * @param verbose is a boolean indicating if information messages must be printed during the function execution.
        * @return a boolean indicating if the sequence was successfully inserted in the graph.
        */
        bool add(const string& seq, const bool verbose = false);

        /** Remove a unitig from the Compacted de Bruijn graph.
        * @param um is a UnitigMap object containing the information of the unitig to remove from the graph.
        * @param verbose is a boolean indicating if information messages must be printed during the execution of the function.
        * @return a boolean indicating if the unitig was successfully removed from the graph.
        */
        bool remove(const const_UnitigMap<U, G>& um, const bool verbose = false);

        /** Create an iterator to the first unitig of the Compacted de Bruijn graph (unitigs are NOT sorted lexicographically).
        * @return an iterator to the first unitig of the graph.
        */
        iterator begin();

        /** Create an constant iterator to the first unitig of the Compacted de Bruijn graph (unitigs are NOT sorted lexicographically).
        * @return a constant iterator to the first unitig of the graph.
        */
        const_iterator begin() const;

        /** Create an iterator to the "past-the-last" unitig of the Compacted de Bruijn graph (unitigs are NOT sorted lexicographically).
        * @return an iterator to the "past-the-last" unitig of the graph.
        */
        iterator end();

        /** Create a constant iterator to the "past-the-last" unitig of the Compacted de Bruijn graph (unitigs are NOT sorted lexicographically).
        * @return a constant iterator to the "past-the-last" unitig of the graph.
        */
        const_iterator end() const;

        /** Return a boolean indicating if the graph is invalid (wrong input parameters/files, error occurring during a method, etc.).
        * @return A boolean indicating if the graph is invalid.
        */
        inline bool isInvalid() const { return invalid; }

        /** Return the length of k-mers of the graph.
        * @return Length of k-mers of the graph.
        */
        inline int getK() const { return k_; }

        /** Return the number of unitigs in the graph.
        * @return Number of unitigs in the graph.
        */
        inline size_t size() const { return v_unitigs.size() + v_kmers.size() + h_kmers_ccov.size(); }

        /** Return a pointer to the graph data. Pointer is nullptr if type of graph data is void.
        * @return A pointer to the graph data. Pointer is nullptr if type of graph data is void.
        */
        inline G* getData() { return data.getData(); }

        /** Return a constant pointer to the graph data. Pointer is nullptr if type of graph data is void.
        * @return A constant pointer to the graph data. Pointer is nullptr if type of graph data is void.
        */
        inline const G* getData() const { return data.getData(); }

    private:

        bool filter(const CDBG_Build_opt& opt);
        bool construct(const CDBG_Build_opt& opt);

        bool addUnitigSequenceBBF(Kmer km, const string& seq, const size_t pos_match_km, const size_t len_match_km,
                                  vector<std::atomic_flag>& locks_mapping, vector<std::atomic_flag>& locks_unitig,
                                  const size_t thread_id);
        size_t findUnitigSequenceBBF(Kmer km, string& s, bool& isIsolated, vector<Kmer>& l_ignored_km_tip);
        bool bwStepBBF(Kmer km, Kmer& front, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;
        bool fwStepBBF(Kmer km, Kmer& end, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, bool check_fp_cand = true) const;

        UnitigMap<U, G> findUnitig(const Kmer& km, const string& s, size_t pos);
        UnitigMap<U, G> findUnitig(const Kmer& km, const string& s, size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h);

        bool addUnitig(const string& str_unitig, const size_t id_unitig);
        void deleteUnitig(const bool isShort, const bool isAbundant, const size_t id_unitig);
        void swapUnitigs(const bool isShort, const size_t id_a, const size_t id_b);
        template<bool is_void> typename std::enable_if<!is_void, bool>::type splitUnitig_(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs,
                                                                                          size_t& v_unitigs_sz, size_t& v_kmers_sz, const vector<pair<int,int>>& sp);
        template<bool is_void> typename std::enable_if<is_void, bool>::type splitUnitig_(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs,
                                                                                         size_t& v_unitigs_sz, size_t& v_kmers_sz, const vector<pair<int,int>>& sp);

        UnitigMap<U, G> find(const Kmer& km, const preAllocMinHashIterator<RepHash>& it_min_h);
        vector<const_UnitigMap<U, G>> findPredecessors(const Kmer& km, const bool extremities_only = false) const;
        vector<const_UnitigMap<U, G>> findSuccessors(const Kmer& km, const size_t limit = 4, const bool extremities_only = false) const;

        inline size_t find(const preAllocMinHashIterator<RepHash>& it_min_h) const {

            const int pos = it_min_h.getPosition();
            return (hmap_min_unitigs.find(Minimizer(&it_min_h.s[pos]).rep()) != hmap_min_unitigs.end() ? 0 : pos - it_min_h.p);
        }

        pair<size_t, size_t> splitAllUnitigs();
        template<bool is_void>
        typename std::enable_if<!is_void, size_t>::type joinUnitigs_(vector<Kmer>* v_joins = nullptr, const size_t nb_threads = 1);
        template<bool is_void>
        typename std::enable_if<is_void, size_t>::type joinUnitigs_(vector<Kmer>* v_joins = nullptr, const size_t nb_threads = 1);

        void createJoinHT(vector<Kmer>* v_joins, KmerHashTable<Kmer>& joins, const size_t nb_threads) const;
        bool checkJoin(const Kmer& a, const const_UnitigMap<U, G>& cm_a, Kmer& b) const;
        void check_fp_tips(KmerHashTable<bool>& ignored_km_tips);
        size_t removeUnitigs(bool rmIsolated, bool clipTips, vector<Kmer>& v);

        size_t joinTips(string filename_MBBF_uniq_kmers, const size_t nb_threads = 1, const bool verbose = false);
        vector<Kmer> extractMercyKmers(BlockedBloomFilter& bf_uniq_km, const size_t nb_threads = 1, const bool verbose = false);

        void writeGFA(string graphfilename, const size_t nb_threads = 1) const;
        void writeFASTA(string graphfilename) const;

        template<bool is_void>
        typename std::enable_if<!is_void, void>::type writeGFA_sequence_(GFA_Parser& graph, KmerHashTable<size_t>& idmap) const;
        template<bool is_void>
        typename std::enable_if<is_void, void>::type writeGFA_sequence_(GFA_Parser& graph, KmerHashTable<size_t>& idmap) const;

        void mapRead(const UnitigMap<U, G>& cc);

        void print() const;


        int k_;
        int g_;

        bool invalid;

        static const int tiny_vector_sz = 2;
        static const int min_abundance_lim = 15;
        static const int max_abundance_lim = 15;

        typedef KmerHashTable<CompressedCoverage_t<U>> h_kmers_ccov_t;
        typedef MinimizerHashTable_2Val hmap_min_unitigs_t;

        typedef typename hmap_min_unitigs_t::iterator hmap_min_unitigs_iterator;
        typedef typename hmap_min_unitigs_t::const_iterator hmap_min_unitigs_const_iterator;

        vector<Unitig<U>*> v_unitigs;
        vector<pair<Kmer, CompressedCoverage_t<U>>> v_kmers;

        hmap_min_unitigs_t hmap_min_unitigs;

        h_kmers_ccov_t h_kmers_ccov;

        BlockedBloomFilter bf;

        wrapperData<G> data;
};

#include "CompactedDBG.tcc"

#endif
