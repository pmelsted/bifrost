#ifndef BFG_HASHID_HPP
#define BFG_HASHID_HPP

#include "CompactedDBG.hpp"
#include "ColorSet.hpp"

template<typename U> using UnitigColorMap = UnitigMap<DataAccessor<U>, DataStorage<U>>;
template<typename U> using const_UnitigColorMap = const_UnitigMap<DataAccessor<U>, DataStorage<U>>;

template<typename Unitig_data_t>
class DataAccessor : public CDBG_Data_t<DataAccessor<Unitig_data_t>, DataStorage<Unitig_data_t>> {

    typedef Unitig_data_t U;

    public:

        DataAccessor(const uint8_t id = 0);

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
