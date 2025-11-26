//
// Created by qq328 on 2025/11/17.
//

#include "Header.h"
#include "de.h"
#include "algorithms/uade_mt.h"
#include "eval_func.h"
#include <omp.h>

int g_function_number;
int g_problem_size = 1000;
unsigned int g_max_num_evaluations;
int g_pop_size = static_cast<int>(std::round(g_problem_size * 18));
int g_memory_size = 6;
double g_arc_rate = 0.0;
double g_p_best_rate = 0.11;


int main() {
  omp_set_num_threads(16);
  srand((unsigned)time(NULL));
  g_max_num_evaluations = 300000000;

  int num_runs = 1;

  for (int i = 0; i < 3; i++) {
    g_function_number = i + 1;
    set_func(g_function_number);
    cout << "\n-------------------------------------------------------" << endl;
    cout << "Function = " << g_function_number << ", Dimension size = " << g_problem_size << "\n" << endl;
    Fitness *bsf_fitness_array = (Fitness*)malloc(sizeof(Fitness) * num_runs);
    Fitness mean_bsf_fitness = 0;
    Fitness std_bsf_fitness = 0;

    for (int j = 0; j < num_runs; j++) {
      searchAlgorithm *alg = new UADE_MT();
      alg->setMinRegion(get_min_region());
      alg->setMaxRegion(get_max_region());
      // cout << get_min_region() << " " << get_max_region() << endl;
      bsf_fitness_array[j] = alg->run();
      cout << j + 1 << "th run, " << "error value = " << bsf_fitness_array[j] << endl;
    }

    for (int j = 0; j < num_runs; j++) mean_bsf_fitness += bsf_fitness_array[j];
    mean_bsf_fitness /= num_runs;

    for (int j = 0; j < num_runs; j++) std_bsf_fitness += pow((mean_bsf_fitness - bsf_fitness_array[j]), 2.0);
    std_bsf_fitness /= num_runs;
    std_bsf_fitness = sqrt(std_bsf_fitness);

    cout  << "\nmean = " << mean_bsf_fitness << ", std = " << std_bsf_fitness << endl;
    free(bsf_fitness_array);
  }
  return 0;
}