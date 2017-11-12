/***************************************************************************
 *   Copyright (C) 2011 by TEIL                                        *
 *                                                                         *
 ***************************************************************************/
#ifndef _CHROMOSOME_H
#define _CHROMOSOME_H

#include <unordered_map>
#include "global.h"
#include "nk-wa.h"
#include <vector>

using namespace std;

class Chromosome {

public:

    static enum Function {
        ONEMAX=0,
        MKTRAP=1,
        FTRAP=2,
        CYCTRAP=3,
        NK=4,
        SPINGLASS=5,
        SAT=6,
        FIX_SPIN=7
    } function;


    Chromosome ();
    Chromosome (const Chromosome& c) {
        count = 0;
        layer = 0;
		gene = NULL;
        init();
		bm_history = c.bm_history;
        *this = c;
    }

    ~Chromosome ();

    int layer;
    int count;
    bool justGHC;

    bool hasSeen() const;

    bool GHC();
    void steepestDescent();

    void init ();
    void init0 ();
    void initR ();

    bool tryFlipping (int index);

    int getHammingDistance (const Chromosome& c) const {
        int result = 0;
        for (int i=0; i<lengthLong; ++i)
            result += myBD.getHammingDistance(gene[i], c.gene[i]);
        return result;
    }

    int getVal (int index) const {
        assert (index >= 0 && index < length);

        int q = quotientLong(index);
        int r = remainderLong(index);

        if ( (gene[q] & (1lu << r)) == 0 )
            return 0;
        else
            return 1;
    }

    void setVal (int index, int val) {

        assert (index >= 0 && index < length);

        if (getVal(index) == val) return;

		++(countFlipped[index]);
        count++;
        justGHC = false;
        setValF(index, val);
        key ^= zKey[index];
    }

    unsigned long getKey () const {
        return key;
    }


    void setValF (int index, int val) {

        assert (index >= 0 && index < length);
        //outputErrMsg ("Index overrange in Chromosome::operator[]");

        int q = quotientLong(index);
        int r = remainderLong(index);

        if (val == 1)
            gene[q] |= (1lu<<r);
        else
            gene[q] &= ~(1lu<<r);

        evaluated = false;
    }

    void flip (int index) {
        assert (index >= 0 && index < length);

		++(countFlipped[index]);
        int q = quotientLong(index);
        int r = remainderLong(index);

        gene[q] ^= (1lu<<r);
        key ^= zKey[index];

        evaluated = false;
    }

	void resetCountFlipped () {
		for (int i = 0; i < length; ++i)
			countFlipped[i] = 0;
	}

    /** real evaluator */
    double evaluate ();

    bool isEvaluated () const;

    bool operator== (const Chromosome & c) const;
    Chromosome & operator= (const Chromosome & c);

    double getFitness ();
    double trap (int u, double high, double low, int trapK) const;
    double oneMax () const;
    double mkTrap (double high, double low) const;
    double cycTrap(double fHigh, double fLow) const;
    double fTrap () const;
    double spinGlass () const;
    double fixSpinGlass () const;
    double nkFitness() const;
    double satFitness() const;


    void printOut () const;
    void shortPrintOut () const;

    int getLength () const;
    void setLength ();

    double getMaxFitness () const;


public:
    static int length;
    static int lengthLong;
    static int nfe;
    static int lsnfe;
    static int hitnfe;
    static bool hit;
    static unordered_map<unsigned long, double> cache;

	vector<char> bm_history;
	int* countFlipped;
	double improveValue;
protected:

    unsigned long *gene;
    double fitness;
    bool evaluated;
    unsigned long long key;

};


#endif
