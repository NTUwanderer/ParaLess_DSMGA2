#ifndef _MYHASH_
#define _HYHASH_

#include <list>
#include "global.h"
#include "dsmga2.h"

class MyHash {
public:
    MyHash (list<int> _mask, list<bool> _bits) {
        mask = _mask;
        bits = _bits;

        setKey();
    }

    void setKey() {
        key = 0;

        list<int>::const_iterator it1 = mask.begin();
        list<bool>::const_iterator it2 = bits.begin();

        while (it1 != mask.end() && it2 != bits.end()) {
            int index = *it1;
            bool bit = *it2;

            key ^= zKey[index];
            if (bit)
                key ^= zKey[index + 2 * Chromosome::length];
            else
                key ^= zKey[index + Chromosome::length];

            ++it1;
            ++it2;
        }
    }

    unsigned long getKey() const {
        return key;
    }


    list<int> mask;
    list<bool> bits;
    unsigned long key;
};

#endif

