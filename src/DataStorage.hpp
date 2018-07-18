#ifndef BFG_DATASTORAGE_HPP
#define BFG_DATASTORAGE_HPP

#include <atomic>

#include "ColorSet.hpp"
#include "CompactedDBG.hpp"

#define BFG_COLOREDCDBG_FORMAT_VERSION 1

template<typename Unitig_data_t> class ColoredCDBG;
template<typename Unitig_data_t> class DataAccessor;
template<typename Unitig_data_t> class DataStorage;

template<typename U> using UnitigColorMap = UnitigMap<DataAccessor<U>, DataStorage<U>>;
template<typename U> using const_UnitigColorMap = const_UnitigMap<DataAccessor<U>, DataStorage<U>>;

template<typename Unitig_data_t = void>
class DataStorage {

        template<typename U> friend class ColoredCDBG;
        template<typename U> friend class DataAccessor;

    public:

        typedef Unitig_data_t U;

        DataStorage();
        DataStorage(const size_t nb_seeds_, const size_t nb_elem_, const vector<string>& color_names_);
        DataStorage(const DataStorage& o);
        DataStorage(DataStorage&& o);

        ~DataStorage();

        void clear();
        void empty();

        DataStorage<U>& operator=(const DataStorage& o);
        DataStorage<U>& operator=(DataStorage&& o);

        inline size_t getNbColors() const { return color_names.size(); }

        const U* getData(const const_UnitigColorMap<U>& um) const;
        U* getData(const UnitigColorMap<U>& um);

        const UnitigColors* getUnitigColors(const const_UnitigColorMap<U>& um) const;
        UnitigColors* getUnitigColors(const UnitigColorMap<U>& um);

        UnitigColors getSubUnitigColors(const const_UnitigColorMap<U>& um) const;
        vector<string> getSubUnitigColorNames(const const_UnitigColorMap<U>& um) const;

        bool write(const string prefix_output_filename, const size_t nb_threads, const bool verbose = false) const;
        bool read(const string& filename_colors, bool verbose = false);

        bool joinUnitigColors(const UnitigColorMap<U>& um_dest, const UnitigColorMap<U>& um_src);

        uint64_t getHash(const UnitigColorMap<U>& um) const;

        UnitigColors* insert();

        size_t getUnitigColorsSize(const size_t nb_threads = 1) const;

    private:

        size_t nb_seeds;

        size_t nb_cs;
        size_t sz_cs;
        size_t sz_free_cs;

        size_t sz_shared_cs;

        uint64_t seeds[256];

        UnitigColors* color_sets;
        UnitigColors::SharedUnitigColors* shared_color_sets;
        atomic<uint64_t>* unitig_cs_link;

        U* data;

        KmerHashTable<size_t> overflow;

        vector<string> color_names;
};

#endif
