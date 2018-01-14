/*
 * dsmga2.h
 *
 *  Created on: May 2, 2011
 *      Author: tianliyu
 */

#ifndef _DSMGA2_H_
#define _DSMGA2_H_

#include <list>
#include "chromosome.h"
#include "statistics.h"
#include "trimatrix.h"
#include "doublelinkedlistarray.h"
#include "fastcounting.h"
#include "bmrecord.h"

class DSMGA2 {
public:
    DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff);

    ~DSMGA2 ();

    void selection ();
    /** tournament selection without replacement*/
    void tournamentSelection();
    void OrigSelection();

    void oneRun (bool output = true);
    int doIt (bool output = true);

    void buildGraph ();
    void buildOrigGraph ();
    void mixing ();
    bool OrigRM(Chromosome&);
    bool restrictedMixing(Chromosome&);
    bool restrictedMixing(Chromosome& ch, list<int>& mask);
    bool restrictedMixing(Chromosome& ch, list<int>& mask, int size);
    bool restrictedMixing(Chromosome& ch, int pos);
    int backMixing(Chromosome& source, list<int>& mask, Chromosome& des);
    int backMixingE(Chromosome& source, list<int>& mask, Chromosome& des);

    bool shouldTerminate ();

    bool foundOptima ();

    int getGeneration () const {
        return generation;
    }

    bool isInP(const Chromosome& ) const;
    bool isInOrigP(const Chromosome& ) const;
    void genOrderN();
    void genOrderELL();

    void showStatistics ();

    bool isSteadyState ();

    void printDist(string fileName) const;
    void checkNfe () const {
        if (initNfe + initBMNfe + rmNfe + bmNfe != Chromosome::nfe) {
            printf ("Not match nfe\n");
        }
    }

//protected:
public:

    bool ERASE;
    bool ADD;
    bool NEW;

    bool NOHIS;
    bool NOBM;

    int ell;                                  // chromosome length
    int nCurrent;                             // population size
    int nOrig;                             // population size
    int nPrev;                             // population size
    bool EQ;
    unordered_map<unsigned long, double> pHash; // to check if a chromosome is in the population
    unordered_map<unsigned long, double> pHashOrig; // to check if a chromosome is in the population

    vector<BMRecord> BMhistory;

    list<int> *masks;
    vector<double> *linkValues;
    list<int> *orig_masks;
    vector<int> selectionIndex;
    vector<int> orig_selectionIndex;
    vector<int> orderN;                             // for random order
    int *orderELL;                             // for random order
    int selectionPressure;
    int maxGen;
    int maxFe;
    int repeat;
    int generation;
    int bestIndex;

    int lastNfe;
    int initNfe;
    int initBMNfe;
    int rmNfe;
    int bmNfe;

    vector<Chromosome> population;
    vector<Chromosome> orig_popu;
    FastCounting* fastCounting;
    FastCounting* orig_fc;

    TriMatrix<double> graph;
    TriMatrix<double> orig_graph;


    double previousFitnessMean;
    Statistics stFitness;

    // methods
    double computeMI(double, double, double, double) const;


    void findClique(int startNode, list<int>& result);
    void findClique(int startNode, list<int>& result, vector<double>& linkValue);
    void buildFastCounting();
    int countXOR(int, int) const;
    int countOne(int) const;

    void findOrigClique(int startNode, list<int>& result);
    void buildOrigFastCounting();
    int countOrigXOR(int, int) const;
    int countOrigOne(int) const;

    size_t findSize(Chromosome&, list<int>&) const;
    size_t findOrigSize(Chromosome&, list<int>&) const;
    size_t findSize(Chromosome&, list<int>&, Chromosome&) const;

    void increaseOne();
    void decreaseOne(int index);


};


#endif /* _DSMGA2_H_ */
