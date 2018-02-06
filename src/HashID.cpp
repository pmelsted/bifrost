#include "HashID.hpp"
#include "ColorSet.hpp"
#include "ColoredCDBG.hpp"

HashID::HashID(const uint8_t hid) : hash_id(hid) {}

void HashID::join(const UnitigMap<HashID>& um_dest, const UnitigMap<HashID>& um_src){

    ColoredCDBG* colored_cdbg = static_cast<ColoredCDBG*>(um_dest.cdbg);

    ColorSet* cs_dest = colored_cdbg->getColorSet(um_dest);
    ColorSet* cs_src = colored_cdbg->getColorSet(um_src);

    const HashID hid(0);

    if (cs_dest != nullptr){ // If a colorset exists for um_dest

        const Kmer new_head = um_dest.strand ? um_dest.getHead() : um_dest.getTail().rep();

        colored_cdbg->joinColors(um_dest, um_src); // Join the color sets

        // TODO: Insert in tombstone if available
        // TODO: if new_head = head, do not insert in overflow

        // Insert new colorset with corresponding head into overflow of k-mers
        colored_cdbg->km_overflow.insert(new_head, ColorSet(*cs_dest));

        cs_dest->setUnoccupied(); // Set um_dest colorset as a tombstone

        um_dest.setData(&hid);
    }

    if (cs_src != nullptr){

        cs_src->setUnoccupied();
        // TODO: if new_head = head, do not insert in overflow
        um_src.setData(&hid);
    }
}

void HashID::sub(const UnitigMap<HashID>& um, HashID& data_dest, const bool last_extraction) const {

    ColoredCDBG* colored_cdbg = static_cast<ColoredCDBG*>(um.cdbg);

    ColorSet cs = colored_cdbg->extractColors(um);

    if (cs.size() != 0){

        const Kmer km = um.getKmer(um.dist);

        // TODO: Insert in tombstone if available
        colored_cdbg->km_overflow.insert(km, cs);

        if (last_extraction){

            const HashID hid(0);

            colored_cdbg->getColorSet(um)->setUnoccupied();

            um.setData(&hid);
        }
    }
}

string HashID::serialize() const { return std::to_string(hash_id); }
