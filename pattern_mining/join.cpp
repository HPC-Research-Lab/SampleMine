#include "join.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>  // generate
#include <algorithm>
#include <functional>  // bind
#include <future>
#include <iterator>  // begin, end, and ostream_iterator
#include <queue>
#include <random>  // mt19937 and uniform_int_distribution
#include <set>
#include <string>
#include <type_traits>
#include <vector>  // vector

#include "db/db.h"
#include "util.h"

using namespace std;

namespace euler::pattern_mining {



#ifdef PROF
  atomic<size_t> memory_load_count = 0;
  atomic<size_t> memory_load_size = 0;
  atomic<size_t> iso_test_count = 0;
  atomic<size_t> auto_test_count = 0;
#endif

  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<> dis(0, 1);

  double random_number() {
    return dis(gen);
  }


  std::pair<std::unordered_set<unsigned long>, std::unordered_set<unsigned long>> join_dummy1 = {};


  bool bsearch(const size_t* x, size_t s, int y) {
    long long low = 0;
    long long high = s - 1;
    while (low <= high) {
      size_t mid = (low + high) / 2;
      if (x[mid] == y) {
        return true;
      }
      else if (x[mid] < y) {
        low = mid + 1;
      }
      else {
        high = mid - 1;
      }
    }
    return false;
  }


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

  std::tuple<vector<vector<shared_ptr<db::MyKV<int>>>>, std::vector<std::shared_ptr<std::vector<std::map<int, map<int, double>>>>>> build_tables(const vector<SGList>& sgls) {
    vector<vector<shared_ptr<db::MyKV<int>>>> H;
    std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>> subgraph_hist;

    int res_size = 2;
    map<SGList, std::pair<vector<shared_ptr<db::MyKV<int>>>, std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>> processed;
    for (auto& sgl : sgls) {
      int n = sgl.sgl->ncols - 1;
      res_size += n - 1;
      auto it = processed.find(sgl);
      if (it != processed.end()) {
        H.push_back(it->second.first);
        subgraph_hist.push_back(it->second.second);
      }
      else {
        vector<shared_ptr<db::MyKV<int>>> ht;
        subgraph_hist.push_back(std::make_shared<std::vector<std::map<int, std::map<int, double>>>>());
        for (int i = 0; i < n; i++) {
          ht.push_back(make_shared<db::MyKV<int>>(sgl.sgl->ncols));
          subgraph_hist.back()->push_back(std::map<int, std::map<int, double>>());

          auto& his = subgraph_hist.back()->back();

          for (auto kv1 : sgl.sgl->keys) {
            auto buf = sgl.sgl->getbuf(kv1.second, sgl.sgl->ncols);
            auto it_buf = buf.begin();
            while (true) {
              for (size_t z = 0; z < it_buf.buffer_size / sgl.sgl->ncols; z++) {
                const int* it_d = it_buf.buffer + z * sgl.sgl->ncols;
                int kk = it_d[i + 1];
                ht.back()->merge(kk, it_d, sgl.sgl->ncols * sizeof(int));
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

        H.push_back(ht);
        processed[sgl] = std::make_pair(ht, subgraph_hist.back());
      }
    }

    return { H, subgraph_hist };
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

  /*
    void update_sampling_weights(double ntot, const SGList& d, const std::vector<std::vector<double>>& original_table_size, std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>>& subgraph_hist) {

      std::vector<double> vecc;
      for (auto it = d.sgl->keys.begin(); it != d.sgl->keys.end(); it++) {
        vecc.push_back(d.sgl->count[it->second]);
      }

      std::sort(vecc.begin(), vecc.end(), std::greater<double>());

      size_t pos = vecc.size() / 10;
      if (pos == 0) pos = 1;
      double threshold = vecc[pos];

      for (auto it = d.sgl->keys.begin(); it != d.sgl->keys.end(); it++) {
        if (d.sgl->count[it->second] < threshold) continue;
        for (auto& [qp, c] : d.sgl->qp_path[it->second]) {
          double s = ntot / (ntot - c);
          size_t length = qp.size();
          for (int j = 0; j < length; j++) {
            for (auto& [k, v] : subgraph_hist[length - j - 1][qp[j][1]]) {
              auto pt = v.find(qp[j][0]);
              pt->second *= s;
            }
          }
        }
      }
      // }

      auto table_size = get_table_size(subgraph_hist);

  #pragma omp parallel for num_threads(_Nthreads)
      for (int i = 0; i < subgraph_hist.size(); i++) {
        for (int j = 0; j < subgraph_hist[i].size(); j++) {
          for (auto& [k1, v1] : subgraph_hist[i][j]) {
            for (auto& [k2, v2] : v1) {
              v2 *= original_table_size[i][j] / table_size[i][j];
            }
          }
        }
      }
    }

  */

  /*
  vector<vector<shared_ptr<db::MyKV<int>>>> build_tables(const vector<SGList> &sgls) {
    vector<vector<shared_ptr<db::MyKV<int>>>> H;
    int res_size = 2;
    for (auto &sgl : sgls) {
      int n = sgl.sgl->ncols - 1;
      res_size += n - 1;
      vector<shared_ptr<db::MyKV<int>>> ht;
      for (int i = 0; i < n; i++) {
        ht.push_back(make_shared<db::MyKV<int>>(sgl.sgl->ncols));
        for (auto kv1 : sgl.sgl->keys) {
          auto buf = sgl.sgl->getbuf(get<0>(kv1.second), sgl.sgl->ncols);
          auto it_buf = buf.begin();
          while (true) {
  #pragma omp parallel for num_threads(_Nthreads)
            for (size_t z = 0; z < it_buf.buffer_size / sgl.sgl->ncols; z++) {
              const int *it_d = it_buf.buffer + z * sgl.sgl->ncols;
              int kk = it_d[i + 1];
              ht.back()->merge(kk, it_d, sgl.sgl->ncols * sizeof(int));
            }
            if (!it_buf.has_next) break;
            it_buf.next();
          }
        }
      }
      H.push_back(ht);
    }
    return H;
  }
  */

}  // namespace euler::pattern_mining
