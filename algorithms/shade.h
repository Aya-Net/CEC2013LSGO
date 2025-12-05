//
// Created by GUO XIUYI on 2025/11/5.
//

#ifndef CEC2014_SHADE_H
#define CEC2014_SHADE_H

#include "../de.h"
#include <random>

class SHADE: public searchAlgorithm {
    public:
        virtual Fitness run();
        void setSHADEParameters();
        void operateCurrentToPBest1BinWithArchive(const vector<Individual> &pop, Individual child, int &target, int &p_best_individual, variable &scaling_factor, variable &cross_rate, const vector<Individual> &archive, int &arc_ind_count, mt19937 &rng);

        int arc_size;
        double arc_rate;
        variable p_best_rate;
        int memory_size;
};

#endif //CEC2014_SHADE_H