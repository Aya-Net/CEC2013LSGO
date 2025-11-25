//
// Created by GUO XIUYI on 2025/11/5.
//

#ifndef CEC2014_UADE_H
#define CEC2014_UADE_H

#include "../de.h"

#include <fstream>

class UADE: public searchAlgorithm {
    public:
        virtual Fitness run();
        void setSHADEParameters();
        void reducePopulationWithSort(vector<Individual> &pop, vector<Fitness> &fitness);
        //void reducePopulationWithSort(vector<Individual> &pop, vector<Fitness> &fitness, deque<Individual> &archive);
        void operateCurrentToPBest1BinWithArchive(const vector<pair<Fitness, Individual>> &pop, const vector<tuple<Fitness, Individual, int>> &pbests, Individual child, int &target, int &p_best_individual, int &r1, int &r2, variable &scaling_factor, variable &cross_rate);
        //void operateCurrentToPBest1BinWithArchive(const vector<Individual> &pop, Individual child, const deque<Individual> &archive, int &target, int &p_best_individual, int &r1, int &r2, variable &scaling_factor, variable &cross_rate);
        int arc_size;
        double arc_rate;
        variable p_best_rate;
        int memory_size;
        int reduction_ind_num;
        ofstream fout,pout;
};

#endif //CEC2014_UADE_H