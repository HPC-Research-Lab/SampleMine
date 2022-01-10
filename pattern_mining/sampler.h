#pragma once


#include <vector>
#include <map>
#include <memory>
#include <random>
#include <chrono>
#include "util.h"
#include "types.h"


namespace euler::pattern_mining {
  std::vector<std::vector<double>> get_table_size(const std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>& subgraph_hist);

  double scale_sampling_param(const SGList& d2, double st);


  class Sampler {
  public:
    virtual std::pair<double, bool> smp_prob(util::span<const int> s, util::span<const int> t, int i, int j, int step, int tid) {
      return {1.0, false};
    }
  };

  class ProportionalSampler : public Sampler {
    std::vector<double> sampling_ratio;
    std::default_random_engine generator;
  public:
    std::vector<std::vector<int64_t>> C;

    ProportionalSampler(const std::vector<double>& sample_ratio) : sampling_ratio(sample_ratio) {
      for (int i=0; i<16; i++) {
        C.push_back(std::vector<int64_t>(sample_ratio.size()+1));
        C[i][0] = (1 << (10 * sample_ratio.size()));
      }

      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      generator = std::default_random_engine(seed);
    }
    inline std::pair<double, bool> smp_prob(util::span<const int> s, util::span<const int> t, int i, int j, int step, int tid) {
      std::binomial_distribution<int64_t> dist(C[tid][step], 1 / sampling_ratio[step] / (1 << 10));
      int64_t n = dist(generator);
      if (n >= 1) {
        C[tid][step + 1] = n;
        return {1/ sampling_ratio[step], false};
      }
      else return {0, true};
    }
  };

  class ProportionalSampler2 : public Sampler {
    std::vector<double> sampling_ratio;
  public:
    ProportionalSampler2(const std::vector<double>& sample_ratio) : sampling_ratio(sample_ratio) {}
    inline std::pair<double, bool> smp_prob(util::span<const int> s, util::span<const int> t, int i, int j, int step, int tid) {
      double pr = 1 / sampling_ratio[step];
      return {pr, util::random_number() >= pr};
    }
  };

  class BudgetSampler : public Sampler {
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>> subgraph_hist;
    std::vector<double> sampling_budget;
  public:
    BudgetSampler(
      std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>> hist, const std::vector<double>& sample_budget) : subgraph_hist(hist), sampling_budget(sample_budget) {
    }
    inline std::pair<double, bool> smp_prob(util::span<const int> s, util::span<const int> t, int i, int j, int step, int tid) {
      double pr = sampling_budget[step] / subgraph_hist->at(step)->at(j).at(t[j + 1]).at(t[0]);
      return {pr, util::random_number() >= pr};
    }
  };


  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>> get_subgraph_hist(const std::vector<SGList>& sgls);

}
