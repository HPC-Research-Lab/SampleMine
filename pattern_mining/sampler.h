#pragma once


#include <vector>
#include <map>
#include <memory>
#include "util.h"
#include "types.h"


namespace euler::pattern_mining {
  std::vector<std::vector<double>> get_table_size(const std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>& subgraph_hist);


  double scale_sampling_param(const SGList& d2, double st);


  class Sampler {
  public:
    virtual bool operator()(util::span<const int> s, util::span<const int> t, int i, int j, int step) {
      return false;
    }
    virtual double smp_prob(util::span<const int> s, util::span<const int> t, int i, int j, int step) {
      return 1.0;
    }
  };

  class ProportionalSampler : public Sampler {
    std::vector<double> sampling_ratio;
  public:
    ProportionalSampler(const std::vector<double>& sample_ratio) : sampling_ratio(sample_ratio) {}
    inline double smp_prob(util::span<const int> s, util::span<const int> t, int i, int j, int step) {
      return 1 / sampling_ratio[step];
    }
  };


  class BudgetSampler : public Sampler {
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>> subgraph_hist;
    std::vector<double> sampling_budget;
  public:
    BudgetSampler(
       std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>> hist,  const std::vector<double> &sample_budget) : subgraph_hist(hist), sampling_budget(sample_budget) {
      }
    inline double smp_prob(util::span<const int> s, util::span<const int> t, int i, int j, int step) {
      return sampling_budget[step] / subgraph_hist->at(step)->at(j).at(t[j+1]).at(t[0]);
    }
  };


  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>> get_subgraph_hist(const std::vector<SGList>& sgls);

}
