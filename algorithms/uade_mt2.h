//
// Created by qq328 on 2025/11/17.
//

#ifndef CEC2013LSGO_UADE_MT_H
#define CEC2013LSGO_UADE_MT_H


#include "../de.h"

#include <fstream>
#include <random>

class UADE_MT: public searchAlgorithm {
public:
  virtual Fitness run();
  void setSHADEParameters();
  void reducePopulationWithSort(vector<Individual> &pop, vector<Fitness> &fitness);
  //void reducePopulationWithSort(vector<Individual> &pop, vector<Fitness> &fitness, deque<Individual> &archive);
  void operateCurrentToPBest1BinWithArchive(const vector<pair<Fitness, Individual>> &pop, const vector<tuple<Fitness, Individual, int>> &pbests, Individual child, int &target, int &p_best_individual, int &r1, int &r2, variable &scaling_factor, variable &cross_rate, std::mt19937 &generator);
  //void operateCurrentToPBest1BinWithArchive(const vector<Individual> &pop, Individual child, const deque<Individual> &archive, int &target, int &p_best_individual, int &r1, int &r2, variable &scaling_factor, variable &cross_rate);
  int arc_size;
  double arc_rate;
  variable p_best_rate;
  int memory_size;
  int reduction_ind_num;
  ofstream fout,pout;
};


#endif //CEC2013LSGO_UADE_MT_H