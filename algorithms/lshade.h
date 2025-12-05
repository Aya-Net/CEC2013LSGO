//
// Created by GUO XIUYI on 2025/11/5.
//

#ifndef CEC2014_LSHADE_H
#define CEC2014_LSHADE_H

#include "../de.h"
#include <random>

class LSHADE: public searchAlgorithm {
    public:
        virtual Fitness run();
        void setSHADEParameters();
        void reducePopulationWithSort(vector<Individual> &pop, vector<Fitness> &fitness);
        void operateCurrentToPBest1BinWithArchive(const vector<Individual> &pop, Individual child, int &target, int &p_best_individual, variable &scaling_factor, variable &cross_rate, const vector<Individual> &archive, int &arc_ind_count, std::mt19937 &rng);

        int arc_size;
        double arc_rate;
        variable p_best_rate;
        int memory_size;
        int reduction_ind_num;
};

#endif //CEC2014_LSHADE_H