#ifndef PYRA
#define PYRA

#include <list>

#include "chromosome.h"
#include "statistics.h"
#include "trimatrix.h"
#include "doublelinkedlistarray.h"
#include "fastcounting.h"

class Pyramid {
    public:
        Pyramid (int n_ell, int n_maxFe, int fffff);
        ~Pyramid ();

    private:
        int ell;                                    // Problem length
        int nCurrent;                               // Population size
        bool EQ;                                    // Equal covering
        int maxFe;                                  // maximum toleration of NFE
        int bestIndex;                              // The best performing chromosome population[bestIndex]

        unordered_map<unsigned long, double> pHash; // Population Hash
        FastCounting* fastCounting;     
        
        list<int> *masks;                           // Restricted Mixing masks
        vector<int> selectionIndex;                 // Record of tournament selection
        
        // Random order purpose
        vector<int> orderN;
        int *orderELL;

        vector<Chromosome>      population;         // The set of chromosomes.
        vecotr<int>             layers;             // The layer index of the correspoding chromosome in the population


        TriMatrix<double> graph;                    // TODO: Make an ILS graph for each layer.

        Statistics stFitness;                       // Statistics to record the fitness of the population
}

#endif
