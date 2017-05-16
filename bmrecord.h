
#ifndef _BMRECORD_H_
#define _BMRECORD_H_

#include <list>
#include "chromosome.h"


class BMRecord {

public:
    BMRecord(const Chromosome& _pattern, const list<int>& _mask, bool _eq, double _rate) {
        pattern = _pattern;
        mask = _mask;
        eq = _eq;
        rate = _rate;

        mask.sort();
    }

    bool operator< (const BMRecord& rec) const {
        //return mask.size() < rec.mask.size();
        return rate > rec.rate;
    } 

    bool operator== (const BMRecord& rec) const {
        if (mask.size() != rec.mask.size())
            return false;
        list<int>::const_iterator it1 = mask.begin();
        list<int>::const_iterator it2 = rec.mask.begin();

        bool result = true;
        while (it1 != mask.end() && it2 != rec.mask.end()) {
            if (*it1 != *it2) {
                result = false;
                break;
            }

            if (pattern.getVal(*it1) != rec.pattern.getVal(*it2)) {
                result = false;
                break;
            }
            ++it1;
            ++it2;
        }

        return result;
    }

public:
    Chromosome pattern;
    list<int> mask;
    bool eq;
    double rate;
};

#endif
