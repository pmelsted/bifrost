#ifndef BFG_HASHID_HPP
#define BFG_HASHID_HPP

#include "CompactedDBG.hpp"

class HashID : public CDBG_Data_t<HashID> {

    public:

        HashID(const uint8_t hid = 0);

        void join(const UnitigMap<HashID>& um_dest, const UnitigMap<HashID>& um_src);
        void sub(const UnitigMap<HashID>& um_src, HashID& new_data, const bool last_extraction) const;

        string serialize() const;

        inline uint8_t get() const { return hash_id; }
        inline void set(const uint8_t hid) { hash_id = hid; }

    private:

        uint8_t hash_id;
};

#endif
