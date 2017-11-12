/*
 * dsmga2.cpp
 *
 *  Created on: May 2, 2011
 *      Author: tianliyu
 */


#include <list>
#include <vector>
#include <algorithm>
#include <iterator>

#include <iostream>
#include "chromosome.h"
#include "dsmga2.h"
#include "fastcounting.h"
#include "statistics.h"


using namespace std;


DSMGA2::DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff) {

    NEW = false;
    NOHIS = false;
    NOBM = false;

    previousFitnessMean = -INF;
    ell = fffff == 7 ? n_ell - 1 : n_ell;
    nPrev = 0;
    nCurrent = n_nInitial;
    nOrig = nCurrent;

    Chromosome::length = ell;
    Chromosome::lengthLong = quotientLong(ell)+1;
    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0;
    Chromosome::lsnfe = 0;
    Chromosome::hitnfe = 0;
    Chromosome::hit = false;

    selectionPressure = 2;
    maxGen = n_maxGen;
    maxFe = n_maxFe;

    graph.init(ell);
    orig_graph.init(ell);

    bestIndex = 0;
    masks = new list<int>[ell];
    orig_masks = new list<int>[ell];
    orderELL = new int[ell];

    selectionIndex.resize(nCurrent);
    orig_selectionIndex.resize(nCurrent);

    orderN.resize(nCurrent);

    orig_fc = new FastCounting[ell];
    for (int i = 0; i < ell; i++)
        orig_fc[i].init(nCurrent);

    fastCounting = new FastCounting[ell];
    for (int i = 0; i < ell; i++)
        fastCounting[i].init(nCurrent);

    pHash.clear();
    pHashOrig.clear();

    for (int i=0; i<nCurrent; ++i) {
        Chromosome ch;
        population.push_back(ch);
        population[i].initR();
    }

    if (GHC) {
        for (int i=0; i < nCurrent; i++)
            population[i].GHC();
    }

    for (int i=0; i < nCurrent; i++) {
        orig_popu.push_back(population[i]);
        double f = population[i].getFitness();
        pHash[population[i].getKey()] = f;
        pHashOrig[population[i].getKey()] = f;
    }
}


DSMGA2::~DSMGA2 () {
    delete []masks;
    delete []orig_masks;
    delete []orderELL;
    delete []fastCounting;
    delete []orig_fc;
}



bool DSMGA2::isSteadyState () {

    if (stFitness.getNumber () <= 0)
        return false;

    if (previousFitnessMean < stFitness.getMean ()) {
        previousFitnessMean = stFitness.getMean () + 1e-6;
        return false;
    }

    return true;
}



int DSMGA2::doIt (bool output) {
    generation = 0;

    while (!shouldTerminate ()) {
        oneRun (output);
    }
    return generation;
}


void DSMGA2::oneRun (bool output) {

    /*
    if (CACHE)
        Chromosome::cache.clear();
        */

    mixing();


    double max = -INF;
    stFitness.reset ();

    for (int i = 0; i < nCurrent; ++i) {
        double fitness = population[i].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        stFitness.record (fitness);

        if (SHOW_POPULATION) {
            if (SHORT_HAND)
                population[i].shortPrintOut();
            else
                population[i].printOut();
            cout << endl;
        }
    }

    if (output)
        showStatistics ();

    ++generation;


}


bool DSMGA2::shouldTerminate () {
    bool
    termination = false;

    if (maxFe != -1) {
        if (Chromosome::nfe > maxFe)
            termination = true;
    }

    if (maxGen != -1) {
        if (generation > maxGen)
            termination = true;
    }


    if (population[0].getMaxFitness() <= stFitness.getMax() )
        termination = true;


    /*
    if (stFitness.getMax() - 1e-10 <= stFitness.getMean() )
        termination = true;
        */

    return termination;

}


bool DSMGA2::foundOptima () {
    return (stFitness.getMax() > population[0].getMaxFitness());
}


void DSMGA2::showStatistics () {

    printf ("Gen:%d  N:%d  Fitness:(Max/Mean/Min):%f/%f/%f ",
            generation, nCurrent, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin ());
    printf ("best chromosome:\n");
    population[bestIndex].printOut();
    printf ("\n");


    fflush(NULL);
}


void DSMGA2::buildOrigFastCounting() {

    if (SELECTION) {
        for (int i = 0; i < nOrig; i++)
            for (int j = 0; j < ell; j++) {
                orig_fc[j].setVal(i, orig_popu[orig_selectionIndex[i]].getVal(j));
            }

    } else {
        for (int i = 0; i < nOrig; i++) {
            for (int j = 0; j < ell; j++) {
                orig_fc[j].setVal(i, orig_popu[i].getVal(j));
            }
        }
    }

}

void DSMGA2::buildFastCounting() {

    if (SELECTION) {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[selectionIndex[i]].getVal(j));
            }

    } else {
        for (int i = 0; i < nCurrent; i++) {
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[i].getVal(j));
            }
        }
    }

}

int DSMGA2::countOrigOne(int x) const {

    int n = 0;

    for (int i=0; i<orig_fc[0].lengthLong; ++i) {
        unsigned long val = 0;

        val = orig_fc[x].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countOrigXOR(int x, int y) const {

    int n = 0;

    for (int i=0; i<orig_fc[0].lengthLong; ++i) {

        unsigned long val = 0;


        val = orig_fc[x].gene[i];

        val ^= orig_fc[y].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countOne(int x) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;

        val = fastCounting[x].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countXOR(int x, int y) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;


        val = fastCounting[x].gene[i];

        val ^= fastCounting[y].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


bool DSMGA2::restrictedMixing(Chromosome& ch, int pos) {

    list<int> mask = masks[pos];

    size_t size;
    if (NEW)
        size = findSize(ch, mask, population[nCurrent-1]);
    else
        size = findSize(ch, mask);

    if (size > (size_t)ell/2)
        size = ell/2;

    // prune mask to exactly size
    while (mask.size() > size)
        mask.pop_back();


    bool taken = restrictedMixing(ch, mask);

    EQ = true;
    //if (taken && mask.size()<log(nCurrent)/log(2)) {
    if (taken) {

        genOrderN();

        double sum = 0;
        if (!NOBM) {
            for (int i=0; i<nCurrent-1; ++i) {
                //if (2*population[orderN[i]].count > ell) continue;

                bool bm = false;
                if (EQ)
                    bm = backMixingE(ch, mask, population[orderN[i]]);
                else
                    bm = backMixing(ch, mask, population[orderN[i]]);

                if (bm) ++sum;
            }

            if (!NOHIS) BMhistory.push_back(BMRecord(ch, mask, EQ, (double)sum/(double)nCurrent));
        }
    }


    return taken;
}

bool DSMGA2::restrictedMixing(Chromosome& ch) {

    int r = myRand.uniformInt(0, ell-1);
    return restrictedMixing(ch, r);

}

bool DSMGA2::backMixing(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial = des;

    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    Chromosome& real = des;

    if (RTRBM) {
        int minDist = trial.getHammingDistance(des);
        for (int i=0; i<RTRW; ++i) {
            int r = myRand.uniformInt(0, nCurrent-1);
            int dist = trial.getHammingDistance(population[r]);
            if (minDist > dist) {
                minDist = dist;
                real = population[r];
            }
        }
    }


    if (trial.getFitness() > real.getFitness()) {
        pHash.erase(real.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        real = trial;
        real.layer++;
        real.bm_history.push_back('1');
        return true;
    }

    real.bm_history.push_back('0');
    return false;

}

bool DSMGA2::backMixingE(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial = des;

    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    Chromosome& real = des;

    if (RTRBM) {
        int minDist = trial.getHammingDistance(des);
        for (int i=0; i<RTRW; ++i) {
            int r = myRand.uniformInt(0, nCurrent-1);
            int dist = trial.getHammingDistance(population[r]);
            if (minDist > dist) {
                minDist = dist;
                real = population[r];
            }
        }
    }

    if (trial.getFitness() > real.getFitness()) {
        pHash.erase(real.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        EQ = false;
        real = trial;
        real.layer++;
        real.bm_history.push_back('1');
        return true;
    }

    if (trial.getFitness() >= real.getFitness()) {
        pHash.erase(real.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        real = trial;
        //real.layer++;
        real.bm_history.push_back('2');
        return true;
    }

    real.bm_history.push_back('0');
    return false;

}


bool DSMGA2::restrictedMixing(Chromosome& ch, list<int>& mask) {

    bool taken = false;
    size_t lastUB = 0;

    for (size_t ub = 1; ub <= mask.size(); ++ub) {

        size_t size = 1;
        Chromosome trial = ch;

        for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

            trial.flip(*it);

            ++size;
            if (size > ub) break;
        }


        if (isInP(trial))
            break;

        if (trial.getFitness() > ch.getFitness()) {
            pHash.erase(ch.getKey());
            pHash[trial.getKey()] = trial.getFitness();

            taken = true;
            ch = trial;
            ch.layer++;

            ch.bm_history.push_back('A');
        }
        else {
            if (trial.getFitness() >= ch.getFitness()) {
                pHash.erase(ch.getKey());
                pHash[trial.getKey()] = trial.getFitness();

                taken = true;
                ch = trial;
                //ch.layer++;
                ch.bm_history.push_back('B');
            }
        }

        if (taken) {
            lastUB = ub;
            break;
        }
    }

    // prune mask for backmixing
    if (lastUB != 0) {
        while (mask.size() > lastUB)
            mask.pop_back();
    }

    return taken;

}

size_t DSMGA2::findOrigSize(Chromosome& ch, list<int>& mask) const {

    DLLA candidate(nOrig);
    for (int i=0; i<nOrig; ++i)
        candidate.insert(i);

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

        int allele = ch.getVal(*it);

        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            if (orig_popu[*it2].getVal(*it) == allele)
                candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;

        ++size;
    }

    return size;


}


size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask) const {

    DLLA candidate(nCurrent);
    for (int i=0; i<nCurrent; ++i)
        candidate.insert(i);

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

        int allele = ch.getVal(*it);

        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            if (population[*it2].getVal(*it) == allele)
                candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;

        ++size;
    }

    return size;


}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask, Chromosome& ch2) const {

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
        if (ch.getVal(*it) == ch2.getVal(*it)) break;
        ++size;
    }
    return size;
}

void DSMGA2::mixing() {

    //bool ADD;
    //int timer = 0;
    //bool second = false;

    increaseOne();
    if (SELECTION)
        selection();
    // really learn model
    buildFastCounting();
    buildGraph();
    for (int i=0; i<ell; ++i)
        findClique(i, masks[i]);

    genOrderELL();

    for (int i=0; i<ell; ++i) {
        int pos = orderELL[i];
        double prevFitness = population[nCurrent-1].getFitness();
        bool taken = restrictedMixing(population[nCurrent-1], pos);
        if (taken && population[nCurrent-1].getFitness() > prevFitness) break;
    }

}

inline bool DSMGA2::isInOrigP(const Chromosome& ch) const {

    unordered_map<unsigned long, double>::const_iterator it = pHashOrig.find(ch.getKey());
    return (it != pHashOrig.end());
}

inline bool DSMGA2::isInP(const Chromosome& ch) const {

    unordered_map<unsigned long, double>::const_iterator it = pHash.find(ch.getKey());
    return (it != pHash.end());
}

inline void DSMGA2::genOrderN() {
    myRand.uniformArray(orderN, nCurrent, 0, nCurrent-1);
}

inline void DSMGA2::genOrderELL() {
    myRand.uniformArray(orderELL, ell, 0, ell-1);
}

void DSMGA2::buildOrigGraph() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOrigOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countOrigXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nOrig - n01 - n10 - n11;

            double p00 = (double)n00/(double)nOrig;
            double p01 = (double)n01/(double)nOrig;
            double p10 = (double)n10/(double)nOrig;
            double p11 = (double)n11/(double)nOrig;

            double linkage;
            linkage = computeMI(p00,p01,p10,p11);
            orig_graph.write(i,j,linkage);

        }
    }

    if (SHOW_GRAPH) {
        for (int i=0; i<ell; ++i) {
            for (int j=0; j<ell; ++j)
                printf("%3f ", graph(i,j));
            printf("\n");
        }
    }

    delete []one;

}
void DSMGA2::buildGraph() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;

            double linkage;
            linkage = computeMI(p00,p01,p10,p11);
            graph.write(i,j,linkage);
        }
    }

    if (SHOW_GRAPH) {
        for (int i=0; i<ell; ++i) {
            for (int j=0; j<ell; ++j)
                printf("%3f ", graph(i,j));
            printf("\n");
        }
    }

    delete []one;

}

// from 1 to ell, pick by max edge
void DSMGA2::findOrigClique(int startNode, list<int>& result) {


    result.clear();

    DLLA rest(ell);
    genOrderELL();
    for (int i=0; i<ell; ++i) {
        if (orderELL[i]==startNode)
            result.push_back(orderELL[i]);
        else
            rest.insert(orderELL[i]);
    }

    double *connection = new double[ell];

    for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
        connection[*iter] = orig_graph(startNode, *iter);

    while (!rest.isEmpty()) {

        double max = -INF;
        int index = -1;
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter) {
            if (max < connection[*iter]) {
                max = connection[*iter];
                index = *iter;
            }
        }

        rest.erase(index);
        result.push_back(index);

        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
            connection[*iter] += orig_graph(index, *iter);
    }


    delete []connection;

}

// from 1 to ell, pick by max edge
void DSMGA2::findClique(int startNode, list<int>& result) {


    result.clear();

    DLLA rest(ell);
    genOrderELL();
    for (int i=0; i<ell; ++i) {
        if (orderELL[i]==startNode)
            result.push_back(orderELL[i]);
        else
            rest.insert(orderELL[i]);
    }

    double *connection = new double[ell];

    for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
        connection[*iter] = graph(startNode, *iter);

    while (!rest.isEmpty()) {

        double max = -INF;
        int index = -1;
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter) {
            if (max < connection[*iter]) {
                max = connection[*iter];
                index = *iter;
            }
        }

        rest.erase(index);
        result.push_back(index);

        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
            connection[*iter] += graph(index, *iter);
    }


    delete []connection;

}

double DSMGA2::computeMI(double a00, double a01, double a10, double a11) const {

    double p0 = a00 + a01;
    double q0 = a00 + a10;
    double p1 = 1-p0;
    double q1 = 1-q0;

    double join = 0.0;
    if (a00 > 1e-10)
        join += a00*log(a00);
    if (a01 > 1e-10)
        join += a01*log(a01);
    if (a10 > 1e-10)
        join += a10*log(a10);
    if (a11 > 1e-10)
        join += a11*log(a11);

    double p = 0.0;
    if (p0 > 1e-10)
        p += p0*log(p0);
    if (p1 > 1e-10)
        p += p1*log(p1);


    double q = 0.0;
    if (q0 > 1e-10)
        q += q0*log(q0);
    if (q1 > 1e-10)
        q += q1*log(q1);

    return -p-q+join;

}


void DSMGA2::selection () {
    tournamentSelection ();
}


// tournamentSelection without replacement
void DSMGA2::OrigSelection () {
    int i, j;

    int randArray[selectionPressure * nOrig];

    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray (randArray + (i * nOrig), nOrig, 0, nOrig - 1);

    for (i = 0; i < nOrig; i++) {

        int winner = 0;
        double winnerFitness = -INF;

        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = orig_popu[challenger].getFitness();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }

        }
        orig_selectionIndex[i] = winner;
    }
}

// tournamentSelection without replacement
void DSMGA2::tournamentSelection () {
    int i, j;

    int randArray[selectionPressure * nCurrent];

    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray (randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);

    for (i = 0; i < nCurrent; i++) {

        int winner = 0;
        double winnerFitness = -INF;

        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = population[challenger].getFitness();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }

        }
        selectionIndex[i] = winner;
    }
}


void DSMGA2::increaseOne () {


    /*
    double max = -INF;
    int bestIndex = 0;

    genOrderN();
    for (int i = 0; i < nCurrent; ++i) {
        double fitness = population[orderN[i]].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = orderN[i];
        }
    }

    double *one = new double[ell];
    for (int i=0; i<ell; ++i) {
        one[i] = 0.0;
        for (int j=0; j<nCurrent; ++j)
            //if (orig_popu[j].getVal(i)==1)
            if (population[bestIndex].getVal(i)==1)
                one[i] += 1.0;
        //one[i] = nCurrent*0.5;
    }

    Chromosome ch;

    do {
        //ch.initR();
        do {
            for (int i=0; i<ell; ++i) 
                ch.setVal(i, myRand.flip(one[i]/(double)nCurrent *0.33333 + 0.33333)? 0:1);
        } while (isInOrigP(ch) || isInP(ch));
        //ch.setVal(i, myRand.flip((one[i]+1.0)/(double)(nCurrent+2.0))? 0:1);
        ch.GHC();
    } while (isInOrigP(ch) || isInP(ch));
    ch.resetCountFlipped();

    delete []one;
    */

    //} while (isInP(ch) || ch.hasSeen());

    /*
    do {
        for (int i=BMhistory.size()-1; i>BMhistory.size()/2; --i) {
            int r = i;
        if (BMhistory[r].eq)
            backMixingE(BMhistory[r].pattern, BMhistory[r].mask, ch);
        else
            backMixing(BMhistory[r].pattern, BMhistory[r].mask, ch);
        }
    } while (!isInP(ch) && !ch.hasSeen());
    */
    Chromosome ch;
    do {
        ch.initR();
        ch.GHC();
    } while (isInP(ch));

    ++nCurrent;
    ++nOrig;

    population.push_back(ch);

    orig_popu.push_back(population[nCurrent-1]);

    pHash[population[nCurrent-1].getKey()] = population[nCurrent-1].getFitness();
    pHashOrig[population[nCurrent-1].getKey()] = population[nCurrent-1].getFitness();

    selectionIndex.emplace_back();
    orig_selectionIndex.emplace_back();
    orderN.push_back(nCurrent-1);

    for (int i = 0; i < ell; i++) {
        fastCounting[i].init(nCurrent);
        orig_fc[i].init(nOrig);
    }
    
/*
            bool temp = SELECTION;
            SELECTION = false;
            OrigSelection();
            buildOrigFastCounting();
            buildOrigGraph();
            SELECTION = temp;


            for (int i=0; i<ell; ++i)
                findOrigClique(i, orig_masks[i]);


            genOrderELL();
            int ssize = nCurrent;
            if (ssize > ell/2) ssize = ell/2;
            for (int i=0; i<ssize; ++i) {

                int r = orderELL[i];

                list<int> mask = orig_masks[r];
                int size = findOrigSize(population[nCurrent-1], mask);
                if (size > (size_t)ell/2)
                    size = ell/2;

                // prune mask to exactly size
                while (mask.size() > size)
                    mask.pop_back();

                bool taken = restrictedMixing(population[nCurrent-1], mask);
            }
  */     

    int size = BMhistory.size();
    int *rrr = new int[size];
    myRand.uniformArray(rrr, size, 0, size-1);
    for (int i=0; i<size; ++i) {
        //int r = i;
        int r = rrr[i];
        if (BMhistory[r].eq)
            backMixingE(BMhistory[r].pattern, BMhistory[r].mask, population[nCurrent-1]);
        else
            backMixing(BMhistory[r].pattern, BMhistory[r].mask, population[nCurrent-1]);
    }

    //population[nCurrent-1].GHC();

    delete []rrr;

    /*
    //cout << "(" << BMhistory.size() << "," << nCurrent << ")" << endl;
    int size = BMhistory.size();
    //int size = (BMhistory.size()+1)/2;
    //int size = ell/2;
    if (size>BMhistory.size()) size = BMhistory.size();
    //int size = nCurrent/2;
    int *rrr = new int[size];
    myRand.uniformArray(rrr, size, 0, size-1);

    //std::sort(BMhistory.begin(), BMhistory.end());


    //cout << "=========================" << endl;
    //population[nCurrent-1].printOut();
    //cout << "-------------------------" << endl;


    for (int i = 0; i<size; ++i) {

        //int r = rrr[i];
        //int r = rrr[i] + size - 1;
        //if (r<0) r=0;
        //if (r>BMhistory.size()-1) r=BMhistory.size()-1;
        //cout << r << endl;
        //cout << BMhistory.size() << endl;
        int r = size - i - 1;
        //int r = i;
        //cout << BMhistory[r].mask.size() << " ";


        EQ = true;
        bool taken;
        if (BMhistory[r].eq)
            taken = backMixingE(BMhistory[r].pattern, BMhistory[r].mask, population[nCurrent-1]);
        else
            taken = backMixing(BMhistory[r].pattern, BMhistory[r].mask, population[nCurrent-1]);


        if (!EQ)
            BMhistory[r].eq = false;

        //double alpha = 1.0/(double)nCurrent;
        double alpha = 0.1;
        if (taken)
            BMhistory[r].rate = (1-alpha)*BMhistory[r].rate + alpha;
        else
            BMhistory[r].rate = (1-alpha)*BMhistory[r].rate;
    }

    delete []rrr;
    */


    /*
            for (int i=0; i<nCurrent; ++i) {
                for (int j=i+1; j<nCurrent; ++j) {
                    int dist = population[i].getHammingDistance(population[j]);
                    if (dist == 0) {
                        decreaseOne(j);
                        --j;
                    }
                }
            }
     */
}


void DSMGA2::decreaseOne(int index) {

    nCurrent--;
    pHash.erase(population[index].getKey());
    population.erase(population.begin()+index);

}
