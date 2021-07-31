#include "sampler.h"
#include "types.h"

namespace euler::pattern_mining {
double scale_sampling_param(const SGList& d2, double st) {
  double max_s = -1;
  for (int i = 0; i < 2; i++) {
    std::map<int, std::map<int, double>> his;
    for (auto kv1 : d2.sgl->keys) {
      auto buf = d2.sgl->getbuf(kv1.second, d2.sgl->ncols);
      auto it_buf = buf.begin();
      while (true) {
        for (size_t z = 0; z < it_buf.buffer_size / d2.sgl->ncols; z++) {
          const int* it_d = it_buf.buffer + z * d2.sgl->ncols;
          int kk = it_d[i + 1];
          if (his.find(kk) == his.end()) his[kk] = std::map<int, double>();
          int type = it_d[0];
          if (his[kk].find(type) == his[kk].end()) his[kk][type] = 0;
          his[kk][type]++;
        }
        if (!it_buf.has_next) break;
        it_buf.next();
      }
    }
    size_t tot_types = 0;
    for (auto& m : his) {
      tot_types += m.second.size();
    }

    double s = st * his.size() / tot_types;
    if (max_s < s) max_s = s;
  }
  return max_s;
}


std::vector<std::vector<double>> get_table_size(const std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>& subgraph_hist) {
  std::vector<std::vector<double>> res(subgraph_hist.size());
  for (int i = 0; i < subgraph_hist.size(); i++) {
    for (int j = 0; j < subgraph_hist[i]->size(); j++) {
      double sum = 0;
      for (auto& [k1, v1] : subgraph_hist[i]->at(j)) {
        for (auto& [k2, v2] : v1) {
          sum += v2;
        }
      }
      res[i].push_back(sum);
    }
  }
  return res;
}


  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>>  get_subgraph_hist(const std::vector<SGList>& sgls) {
    auto subgraph_hist = std::make_shared<std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>>();

    int res_size = 2;
    std::map<SGList, std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>> processed;
    for (auto& sgl : sgls) {
      int n = sgl.sgl->ncols - 1;
      res_size += n - 1;
      auto it = processed.find(sgl);
      if (it != processed.end()) {
        subgraph_hist->push_back(it->second);
      }
      else {
        subgraph_hist->push_back(std::make_shared<std::vector<std::map<int, std::map<int, double>>>>());
        for (int i = 0; i < n; i++) {
          subgraph_hist->back()->push_back(std::map<int, std::map<int, double>>());

          auto& his = subgraph_hist->back()->back();

          for (auto kv1 : sgl.sgl->keys) {
            auto buf = sgl.sgl->getbuf(kv1.second, sgl.sgl->ncols);
            auto it_buf = buf.begin();
            while (true) {
              for (size_t z = 0; z < it_buf.buffer_size / sgl.sgl->ncols; z++) {
                const int* it_d = it_buf.buffer + z * sgl.sgl->ncols;
                int kk = it_d[i + 1];
                if (his.find(kk) == his.end()) his[kk] = std::map<int, double>();
                int type = it_d[0];
                if (his[kk].find(type) == his[kk].end()) his[kk][type] = 0;
                his[kk][type]++;
              }
              if (!it_buf.has_next) break;
              it_buf.next();
            }
          }
        }
        processed[sgl] = subgraph_hist->back();
      }
    }

    return subgraph_hist;
  }




}