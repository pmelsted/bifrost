#ifndef BFG_DATAACCESSOR_HPP
#define BFG_DATAACCESSOR_HPP

#include "CompactedDBG.hpp"
#include "ColorSet.hpp"

/** @file src/DataAccessor.hpp
* Interface for the class DataAccessor. The purpose of a DataAccessor object is to provide access
* to the colors and the data associated with a unitig of a ColoredCDBG. Code snippets using this
* interface are provided in snippets/test.cpp.
*/

template<typename U> using UnitigColorMap = UnitigMap<DataAccessor<U>, DataStorage<U>>;
template<typename U> using const_UnitigColorMap = const_UnitigMap<DataAccessor<U>, DataStorage<U>>;

/** @class DataAccessor
* @brief Interface to access the colors and the data associated with a unitig of a ColoredCDBG.
* The class as one template parameter: the type of data associated with the unitigs of the graph.
*/
template<typename Unitig_data_t>
class DataAccessor : public CDBG_Data_t<DataAccessor<Unitig_data_t>, DataStorage<Unitig_data_t>> {

    typedef Unitig_data_t U;

    public:

        /** Constructor (set up an empty DataAccessor).
        */
        DataAccessor(const uint8_t id = 0);

        /** Get the value of a DataAccessor.
        *
        */
        inline uint8_t get() const { return da_id; }
        inline void set(const uint8_t id) { da_id = id; }

        const U* getData(const const_UnitigColorMap<U>& um) const;
        U* getData(const UnitigColorMap<U>& um) const;

        const UnitigColors<U>* getUnitigColors(const const_UnitigColorMap<U>& um) const;
        UnitigColors<U>* getUnitigColors(const UnitigColorMap<U>& um) const;

        UnitigColors<U> getSubUnitigColors(const const_UnitigColorMap<U>& um) const;
        vector<string> getSubUnitigColorNames(const const_UnitigColorMap<U>& um) const;

        static void join(const UnitigColorMap<U>& um_dest, const UnitigColorMap<U>& um_src);
        static void sub(const UnitigColorMap<U>& um_src, DataAccessor<Unitig_data_t>* new_data, const bool last_extraction);

        string serialize() const;

    private:

        uint8_t da_id;
};

#endif
