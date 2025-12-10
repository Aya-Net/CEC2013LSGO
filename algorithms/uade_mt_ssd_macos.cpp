//
// Created by qq328 on 2025/11/17.
//

#include "uade_mt_ssd_macos.h"
#include "../eval_func.h"
#include "../Benchmarks.h"
#include "../utils/multi_thread_race.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <atomic>
#include <iostream>
#include <chrono>
#include <algorithm>
/*
  SHADE 1.1 implemented by C++ for Special Session & Competition on Real-Parameter Single Objective Optimization at CEC-2014

  Version: 1.1   Date: 9/Jun/2014
  Written by Ryoji Tanabe (rt.ryoji.tanabe [at] gmail.com)
*/


#include <algorithm>
#include <functional>
#include <unordered_set>

int macos_fallocate(int fd, off_t offset, off_t len) {
  fstore_t store;
  store.fst_flags = F_ALLOCATEALL;
  store.fst_posmode = F_PEOFPOSMODE;
  store.fst_offset = offset;
  store.fst_length = len;
  store.fst_bytesalloc = 0;

  if (fcntl(fd, F_PREALLOCATE, &store) == -1) {
    return -1;
  }


  off_t new_size = offset + len;
  if (ftruncate(fd, new_size) == -1) {
    return -1;
  }

  return 0;
}

static bool preallocateFilePosix(const std::string &path, uint64_t size_bytes, int &outFd) {
  int fd = open(path.c_str(), O_CREAT | O_RDWR, 0644);
  if (fd < 0) return false;
  int rc = macos_fallocate(fd, 0, (off_t)size_bytes);
  if (rc != 0) {
    if (ftruncate(fd, (off_t)size_bytes) != 0) {
      close(fd);
      return false;
    }
  }
  outFd = fd;
  return true;
}


static ssize_t writeAtOffset(int fd, const void *buf, size_t bytesToWrite, uint64_t offset) {
  ssize_t w = pwrite(fd, buf, bytesToWrite, (off_t)offset);
  return w;
}

static ssize_t readAtOffset(int fd, void *buf, size_t bytesToRead, uint64_t offset) {
  ssize_t r = pread(fd, buf, bytesToRead, (off_t)offset);
  return r;
}


int addIndividualToDisk(const variable* data, int pop_fd,
                        std::atomic<std::size_t> &next_slot,
                        int problem_size) {
  std::size_t idx = next_slot.fetch_add(1, std::memory_order_relaxed);
  uint64_t pos = (uint64_t) idx * sizeof(variable) * problem_size;
  size_t bytes = sizeof(variable) * (size_t)problem_size;
  ssize_t w = writeAtOffset(pop_fd, data, bytes, pos);
  if (w != (ssize_t)bytes) return -1;

  return static_cast<int>(idx);
}

void loadIndividual(int idx, variable* buffer, int pop_fd, int problem_size) {
  uint64_t pos = (uint64_t) idx * sizeof(variable) * problem_size;
  size_t bytes = sizeof(variable) * (size_t)problem_size;
  ssize_t r = readAtOffset(pop_fd, buffer, bytes, pos);
  if (r != (ssize_t)bytes) {
    // optional: 处理读取错误
    memset(buffer, 0, bytes);
  }
}

void saveIndividual(int idx, const variable* buffer, int pop_fd, int problem_size) {
  uint64_t pos = (uint64_t) idx * sizeof(variable) * problem_size;
  size_t bytes = sizeof(variable) * (size_t)problem_size;
  ssize_t w = writeAtOffset(pop_fd, buffer, bytes, pos);
  (void)w;
}

Fitness UADE_MT::run()
{
  initializeParameters();
  setSHADEParameters();
  cout << scientific << setprecision(8);
  function_number -= 1;
  string fname = string("mt4_") + to_string(function_number / 5) + "_" + to_string((function_number % 5));
  fout.open(fname + "_f.txt");
  fout << scientific << setprecision(8);
  //cout << scientific << setprecision(8);
  function_number /= 5;
  string pop_file_path = "population.bin";
  if (access(pop_file_path.c_str(), F_OK) == 0) {
    if (unlink(pop_file_path.c_str()) != 0) {
      std::cerr << "Warning: failed to remove existing `population.bin`: " << strerror(errno) << "\n";
    }
  }

  uint64_t slot_bytes = (uint64_t)sizeof(variable) * problem_size;
  uint64_t PREALLOC_BYTES = 770ULL * 1024ULL * 1024ULL * 1024ULL; // 600GB
  atomic<std::size_t> next_slot(0);
  uint64_t capacity = PREALLOC_BYTES / slot_bytes;
  if (capacity == 0) capacity = 1;

  int pop_fd = -1;
  if (!preallocateFilePosix(pop_file_path, capacity * slot_bytes, pop_fd)) {
    std::cerr << "Failed to create/preallocate `population.bin`\n";
    return 0;
  }
  constexpr int batch_size = 100;
  vector<pair<Fitness, int>> pop;
  vector<Individual> children(batch_size, nullptr);
  // variable child[problem_size];
  vector<Fitness> children_fitness(pop_size, 0);
  // Fitness child_fitness;
  variable bsf_solution[problem_size];
  vector<mt19937> generators;
  vector<std::uniform_real_distribution<double>> distributions;
  TournamentSampler race(batch_size);

  random_device rd;
  Fitness bsf_fitness;
  atomic nfes = 0;
  vector<Benchmarks*> benchmarks(100, nullptr);
  for (int i = 0; i < batch_size; i++) {
    benchmarks[i] = generateFuncObj(function_number);
  }
  for (int i = 0; i < batch_size; ++i) {

    benchmarks[i]->nextRun();
    children[i] = (variable*)malloc(sizeof(variable) * problem_size);
    generators.emplace_back(rd());
    distributions.emplace_back(0.0, 1.0);
  }
  //unordered_set<DoubleRegion> popset;
  //int duplicate_count = 0;

  for (int i = 0; i < pop_size; ++i)
  {
    //initialize population
    variable* indiv = makeNewIndividual();
    Fitness fit = benchmarks[0]->compute(indiv);
    int disk_index = addIndividualToDisk(indiv, pop_fd, next_slot, problem_size);
    free(indiv);
    pop.emplace_back(fit, disk_index);
    //children.push_back((variable *)malloc(sizeof(variable) * problem_size));

    /*DoubleRegion dr;
    dr.p = (variable *)malloc(sizeof(variable) * problem_size);
    for (int j = 0; j < problem_size; j++)
      dr.p[j] = pop.back().second[j];
    dr.size = problem_size;
    dr.nfes = nfes;
    dr.parent = 0;
    popset.emplace(dr);*/
    //evaluate the initial population
    // cec14_test_func(pop[i].second, &pop[i].first, problem_size, 1, function_number);

    if ((pop[i].first - optimum) < epsilon)
      pop[i].first = optimum;
    race.insert(pop[i].first, i);
    if (i == 0 || fit < bsf_fitness) {
      bsf_fitness = fit;
      // bsf_solution 仍然需要一次性读回或用 indiv 直接拷
      loadIndividual(disk_index, bsf_solution, pop_fd, problem_size);
    }

    ++nfes;
    // fout << pop[i].first << ',' << nfes << ',' << i << ',' << -1 << ',' << pop_size << ',' << 0.5 << ',' << 0.5 << endl;
    if (nfes >= max_num_evaluations)
      break;
  }

  int num_success_params;
  variable mu_sf, mu_cr;
  int mu_ps;
  vector<variable> success_sf;
  vector<variable> success_cr;
  vector<int> success_ps;
  vector<variable> dif_fitness;

  // the contents of M_f and M_cr are all initialiezed 0.5
  vector<variable> memory_sf(memory_size, 0.5);
  vector<variable> memory_cr(memory_size, 0.5);
  vector<variable> memory_ps(memory_size, pop_size);

  variable temp_sum_sf;
  variable temp_sum_cr;
  variable temp_sum_ps;
  variable sum;
  variable weight;

  mutex mtx;

  //memory index counter
  int memory_pos = 0;

  //for new parameters sampling
  int random_selected_period;
  vector<variable> pop_sf(batch_size);
  vector<variable> pop_cr(batch_size);
  vector<int> pop_ps(batch_size);

  //for current-to-pbest/1
  int p_best_ind, base[batch_size], r1[batch_size], r2[batch_size];
  int p_num = round(pop_size * p_best_rate);
  vector<tuple<Fitness, int, int>> pbests;
  for (int i = 0; i < p_num; ++i)
  {
    pbests.push_back({pop[i].first, pop[i].second,i});
  }
  make_heap(pbests.begin(), pbests.end());
  for (int i = p_num; i < pop_size; ++i)
  {
    if (pop[i].first < get<0>(pbests[0]))
    {
      pop_heap(pbests.begin(), pbests.end());
      get<0>(pbests.back()) = pop[i].first;
      get<1>(pbests.back()) = pop[i].second;
      get<2>(pbests.back()) = i;
      push_heap(pbests.begin(), pbests.end());
    }
  }/**/

  //------------------------------------------
  // for linear population size reduction
  int max_pop_size = pop_size * 20;
  int min_pop_size = 100;



#define CHRONO
#ifdef CHRONO
  using Clock = std::chrono::high_resolution_clock;
  auto t_loop_start = Clock::now();
  auto last_time = Clock::now();
#endif
  bool terminate_flag = false;
  //main loop
  while (nfes < max_num_evaluations) {
#ifdef CHRONO
    auto t_parallel_start = Clock::now();
    chrono::time_point<chrono::system_clock, chrono::system_clock::duration> t_p0_start, t_p1_start, t_p2_start;
#endif
    pop.resize(pop.size() + batch_size, {0, 0});
    int last_nfes = nfes;
    // cout << last_nfes << endl;
    //for (int target = 0; target < pop_size; target++)
#pragma omp parallel for private(mu_sf, mu_cr, mu_ps, random_selected_period, p_best_ind) schedule(dynamic)
    for (int target = 0; target < batch_size; target++) {
      //In each generation, CR_i and F_i used by each individual x_i are generated by first selecting an index r_i randomly from [1, H]
      random_selected_period = generators[target]() % memory_size;
      mu_sf = memory_sf[random_selected_period];
      mu_cr = memory_cr[random_selected_period];
      mu_ps = memory_ps[random_selected_period];
      //generate CR_i and repair its value
      if (mu_cr == -1)
      {
        pop_cr[target] = 0;
      } else {
        pop_cr[target] = gauss(mu_cr, 0.1);
        if (pop_cr[target] > 1) pop_cr[target] = 1;
        else if (pop_cr[target] < 0) pop_cr[target] = 0;
      }

      //generate F_i and repair its value
      do {
        pop_sf[target] = cauchy_g(mu_sf, 0.1);
      } while (pop_sf[target] <= 0);
      if (pop_sf[target] > 1) pop_sf[target] = 1;

      //generate PS_i and repair its value
      //pop_ps[target] = pop_size; //fixed-size
      pop_ps[target] = (int)gauss(mu_ps, 10);
      if (pop_ps[target] < min_pop_size) pop_ps[target] = min_pop_size;
      else if (pop_ps[target] > (int)pop.size()) pop_ps[target] = (int)pop.size();
      else if (pop_ps[target] > max_pop_size) pop_ps[target] = max_pop_size;/**/

      //p-best individual is randomly selected from the top pop_size *  p_i members
      p_best_ind = rand() % p_num;

      int tournament = max(1, (int)(round((int)pop.size() / pop_ps[target])));

      p_best_ind = race.sample_top_y_min(p_best_rate, tournament, generators[target], distributions[target]).id;

            base[target] = race.sample_top_y_min(1, tournament, generators[target], distributions[target]).id;
            do {
              r1[target] = race.sample_top_y_min(1, tournament, generators[target], distributions[target]).id;
            } while (r1[target] == base[target]);
            do {
              r2[target] = race.sample_top_y_min(1, tournament, generators[target], distributions[target]).id;
            } while (r1[target] == base[target] || r1[target] == r2[target]);/* UADE(T)*/
      /*tournament = max(1, tournament);
      int j_base = target;
      int j_r1, j_r2;
      do {
        j_r1 = generators[target]() % batch_size;
      } while (j_base == j_r1);
      do {
        j_r2 = generators[target]() % batch_size;
      } while (j_base == j_r2 || j_r1 == j_r2);
      base[target] = race.sample_min_dpt(tournament, j_base, generators[target], distributions[target]).id;
      r1[target] = race.sample_min_dpt(tournament, j_r1, generators[target], distributions[target]).id;
      r2[target] = race.sample_min_dpt(tournament, j_r2, generators[target], distributions[target]).id;/* UADE(DPT)*/

      operateCurrentToPBest1BinWithArchive(pop, pbests, children[target], base[target], p_best_ind, r1[target], r2[target], pop_sf[target], pop_cr[target], generators[target], pop_fd);

      // evaluate the children's fitness values
      children_fitness[target] = benchmarks[target]->compute(children[target]);
      ++nfes;
      // cec14_test_func(child, &child_fitness, problem_size, 1, function_number);

      if (children_fitness[target] - optimum < epsilon) {
        children_fitness[target] = optimum;
      }
      if (children_fitness[target] - optimum < epsilon)
        children_fitness[target] = optimum;


      //insert child to population
      int new_idx = last_nfes + target;
      saveIndividual(new_idx, children[target], pop_fd, problem_size);
      // addIndividualToDisk(children[target], pop_fd, next_slot, problem_size);
      pop[new_idx] = {children_fitness[target], new_idx};
      race.insert(children_fitness[target], new_idx);

      lock_guard lock(mtx);
      if (children_fitness[target] < bsf_fitness)
      {
        bsf_fitness = children_fitness[target];
        for (int j = 0; j < problem_size; j++)
          bsf_solution[j] = children[target][j];
      }
      //generation alternation
      if (children_fitness[target] < pop[base[target]].first)
      {
        // fout << child_fitness << ',' << nfes << ',' << pop.size() << ',' << base << ',' << pop_ps[target] << ',' << pop_sf[target] << ',' << pop_cr[target] << endl;
        if (nfes < max_num_evaluations) {
          //successful parameters are preserved in S_F and S_CR
          success_sf.push_back(pop_sf[target]);
          success_cr.push_back(pop_cr[target]);
          success_ps.push_back(pop_ps[target]);
          dif_fitness.push_back(fabs(pop[base[target]].first - children_fitness[target]));

        }
      }
    }
#ifdef CHRONO
    auto t_parallel_end = Clock::now();
#endif
    if (bsf_fitness - optimum < epsilon) {
      bsf_fitness = optimum;
      // break;
      terminate_flag = true;
    }
    if (nfes % 100000 == 0) {
      cout  << "Function " << function_number
            << ", Evaluation " << nfes
            <<", Population Size " << pop.size()
            << ", Best " << bsf_fitness
            << ", Time " << (long long)((chrono::duration<double>(t_parallel_end - t_loop_start).count())) << " s"
            << ", Period " << (long long)((chrono::duration<double>(t_parallel_end - last_time).count() * 1000)) << "ms"<< endl;
      fout  << "Function " << function_number
      << ", Evaluation " << nfes
      <<", Population Size " << pop.size()
      << ", Best " << bsf_fitness
      << ", Time " << (long long)((chrono::duration<double>(t_parallel_end - t_loop_start).count())) << " s"
      << ", Period " << (long long)((chrono::duration<double>(t_parallel_end - last_time).count() * 1000)) << "ms"<< endl;
      last_time = Clock::now();
    }

    if (nfes >= max_num_evaluations)
      break;

    if (terminate_flag) break;

    // if numeber of successful parameters > 0, historical memories are updated
    num_success_params = success_sf.size();
    if (num_success_params > 0) {
      memory_sf[memory_pos] = 0;
      memory_cr[memory_pos] = 0;
      memory_ps[memory_pos] = 0;
      temp_sum_sf = 0;
      temp_sum_cr = 0;
      temp_sum_ps = 0;
      sum = 0;

      for (int i = 0; i < num_success_params; i++) sum += dif_fitness[i];

      //weighted lehmer mean
      for (int i = 0; i < num_success_params; i++) {
        weight = dif_fitness[i] / sum;

        memory_sf[memory_pos] += weight * success_sf[i] * success_sf[i];
        temp_sum_sf += weight * success_sf[i];

        memory_cr[memory_pos] += weight * success_cr[i] * success_cr[i];
        temp_sum_cr += weight * success_cr[i];

        memory_ps[memory_pos] += weight * success_ps[i] * success_ps[i];
        temp_sum_ps += weight * success_ps[i];
      }

      memory_sf[memory_pos] /= temp_sum_sf;

      if (temp_sum_cr == 0 || memory_cr[memory_pos] == -1) memory_cr[memory_pos] = -1;
      else memory_cr[memory_pos] /= temp_sum_cr;

      if (temp_sum_ps > 0)
      {
        memory_ps[memory_pos] /= temp_sum_ps;
      }
      //if ((int)pop.size() < max_pop_size) memory_ps[memory_pos] = max_pop_size; //uPADE11

      //pop_size = memory_ps[memory_pos]; // for test

      //increment the counter
      memory_pos++;
      if (memory_pos >= memory_size) memory_pos = 0;

      //clear out the S_F, S_CR and delta fitness
      success_sf.clear();
      success_cr.clear();
      success_ps.clear();
      dif_fitness.clear();
    }

    // linear population reduction
    pop_size = memory_ps[memory_pos];
    //pop_size = round((((min_pop_size - max_pop_size) / (double)max_num_evaluations) * nfes) + max_pop_size); //original
    //pop_size = round((((min_pop_size - max_pop_size) / (double)max_num_evaluations) * nfes * 2) + max_pop_size); //half
    //pop_size = round((((min_pop_size - max_pop_size) / (double)max_num_evaluations) * nfes / 2.) + max_pop_size); //double
    //pop_size = round((((min_pop_size - max_pop_size) / (double)max_num_evaluations) * nfes / 4.) + max_pop_size); //quadruple

    //pop_sf.resize(pop_size);
    //pop_cr.resize(pop_size);
    //pop_ps.resize(pop_size);

    while (p_num > max(2., round(pop_size * p_best_rate))) {
      pop_heap(pbests.begin(), pbests.end());
      pbests.pop_back();
      --p_num;
    }/**/
  }
#ifdef CHRONO
  auto t_loop_end = Clock::now();
  cout <<"Evalution numbers: " << nfes << "\nBest Fitness: " << bsf_fitness << "\nTotal time: " << (int)(chrono::duration<double>(t_loop_end - t_loop_start).count() * 1000) << " ms" << endl;
  fout <<"Evalution numbers: " << nfes << "\nBest Fitness: " << bsf_fitness << "\nTotal time: " << (int)(chrono::duration<double>(t_loop_end - t_loop_start).count() * 1000) << " ms" << endl;
#endif

  for (int i = 0; i < children.size(); i++)
    free(children[i]);

  for (int i = 0; i < p_num; i++) {
    //fout << pbests[0].first << endl;
    pop_heap(pbests.begin(), pbests.end());
    pbests.pop_back(); // removes the largest individual
  }

  /*for (auto itr=popset.begin(); itr!=popset.end(); itr++) {
    free(itr->p);
  }*/
  fout.close();
  if (pop_fd >= 0) close(pop_fd);
  for (int i = 0; i < batch_size; ++i)
    delete benchmarks[i];
  return bsf_fitness - optimum;
}

void UADE_MT::operateCurrentToPBest1BinWithArchive(const vector<pair<Fitness, int>> &pop,
                                      const vector<tuple<Fitness, int, int>> &pbests,
                                      Individual child, int &target, int p_best_individual,
                                      int &r1,
                                      int &r2,
                                      variable &scaling_factor,
                                      variable &cross_rate,
                                      mt19937 &generator,
                                      int pop_fd
  )
{
  vector<variable> buf_target(problem_size);
  vector<variable> buf_pbest(problem_size);
  vector<variable> buf_r1(problem_size);
  vector<variable> buf_r2(problem_size);

    // std::shared_lock<std::shared_mutex> rlock(popfile_mtx);
    loadIndividual(pop[target].second, buf_target.data(), pop_fd, problem_size);
    loadIndividual(pop[p_best_individual].second, buf_pbest.data(), pop_fd, problem_size);
    loadIndividual(pop[r1].second, buf_r1.data(), pop_fd, problem_size);
    loadIndividual(pop[r2].second, buf_r2.data(), pop_fd, problem_size);

  int random_variable = generator() % problem_size;
  for (int i = 0; i < problem_size; i++) {
    if ((((double)generator() / generator.max()) < cross_rate) || (i == random_variable)) {
      //child[i] = pop[target].second[i] + scaling_factor * (get<1>(pbests[p_best_individual])[i] - pop[target].second[i]) + scaling_factor * (pop[r1].second[i] - pop[r2].second[i]);
      child[i] = buf_target[i] + scaling_factor * (buf_pbest[i] - buf_target[i]) + scaling_factor * (buf_r1[i] - buf_r2[i]);
    } else  {
      child[i] = buf_target[i];
    }
  }
  //If the mutant vector violates bounds, the bound handling method is applied
  modifySolutionWithParentMedium(child,  buf_target.data());
}

void UADE_MT::setSHADEParameters() {
  arc_rate = g_arc_rate;
  arc_size = (int)round(pop_size * arc_rate);
  p_best_rate = g_p_best_rate;
  memory_size = g_memory_size;
}

