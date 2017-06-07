#include "pyramid.h"

/*
    Initialization of the pyramid structure

    Build one choromsome in the population first.
*/
Pyramid :: Pyramid (int n_ell, int n_maxFe, int fffff)
{   ell = n_ell;
    maxFe = n_maxFe;
    
    Chromosome::length = n_ell;
    Chromosome::lengthLong = quotientLong(ell)+1;
    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0;
    Chromosome::lsnfe = 0;
    Chromosome::hitnfe = 0;
    Chromosome::hit = false;

    graph.init(ell);

    bestIndex = -1; // Not found
    masks = new list<int>[ell];
    orderELL = new int[ell];

    Chromosome ch;
    ch.initR();     // Random initialization
    population.push_back(ch);
    nCurrent = 1;
    
    double f = ch.getFitness();
    pHash[population[0].getKey()] = f;


    fastCounting = new FastCounting[ell];
    for (int i = 0; i < ell; ++i)
        fastCounting[i].init(nCurrent);

    pHash.clear();
}

Pyramid::~Pyramid()
{   delete[] masks;
    delete[] orderELL;
    delete[] fastCounting;
}


