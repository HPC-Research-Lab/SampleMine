#include "filter.h"

#include <fcntl.h>
#include <limits.h>
#include <omp.h>
#include <unistd.h>

#include <iostream>
#include <set>
#include <string>
#include <vector>
#include "match.h"

using namespace std;
namespace euler::pattern_mining {

  const static int _Nthreads = atoi(getenv("OMP_NUM_THREADS"));

  void filter(SGList& d, double threshold) {
    if (!d.sgl->distinct_vertices.empty()) {
      auto it = d.sgl->keys.begin();
      while (it != d.sgl->keys.end()) {
        size_t file_id = it->second;
        if (d.sgl->mni_met.size() > 0 && d.sgl->mni_met[file_id]) {
          it++;
          continue;
        }
        auto& distincts = d.sgl->distinct_vertices[file_id];

        double mni = (double)INT_MAX;
        for (int j = 0; j < distincts.size(); j++) {
          auto& t = distincts[j];
          std::string mni_fname = d.sgl->DBPath + "/mni_" + std::to_string(file_id) + "_" + std::to_string(j) + ".dat";
          bool file_exist = false;
          if (access(mni_fname.c_str(), F_OK) != -1) {
            std::cout << "wrong in filter" << std::endl;
            ifstream fin(mni_fname);
            int v;
            while (fin >> v) {
              t.insert(v);
            }
            file_exist = true;
          }
          if (mni > t.size()) mni = t.size();
          t.clear();
          if (file_exist) {
            remove(mni_fname.c_str());
          }
        }
        if (mni < threshold) {
          //d.patterns.erase(d.sgl->buf[file_id][0]);
          vector<int>().swap(d.sgl->buf[file_id]);
          if (d.sgl->qp_path.size() > 0) {
            d.sgl->qp_path[file_id].clear();
          }
          if (d.sgl->distinct_vertices.size() > 0) {
            d.sgl->distinct_vertices[file_id].clear();
          }

          std::string fname = d.sgl->DBPath + "/" + std::to_string(file_id) + ".dat";
          if (access(fname.c_str(), F_OK) != -1) {
            remove(fname.c_str());
          }
          it = d.sgl->keys.erase(it);
        }
        else {
          ++it;
        }
      }
    }
  }


  void test_and_filter(const graph::Graph& g, SGList& d, double threshold) {
    if (!d.sgl->distinct_vertices.empty()) {
      auto it = d.sgl->keys.begin();
      while (it != d.sgl->keys.end()) {
        size_t file_id = it->second;
        if (d.sgl->mni_met.size() > 0 && d.sgl->mni_met[file_id]) {
          it++;
          continue;
        }
        auto& distincts = d.sgl->distinct_vertices[file_id];

        double mni = (double)INT_MAX;
        for (int j = 0; j < distincts.size(); j++) {
          auto& t = distincts[j];
          std::string mni_fname = d.sgl->DBPath + "/mni_" + std::to_string(file_id) + "_" + std::to_string(j) + ".dat";
          bool file_exist = false;
          if (access(mni_fname.c_str(), F_OK) != -1) {
            std::cout << "wrong in filter" << std::endl;
            ifstream fin(mni_fname);
            int v;
            while (fin >> v) {
              t.insert(v);
            }
            file_exist = true;
          }
          if (mni > t.size()) mni = t.size();
          t.clear();
          if (file_exist) {
            remove(mni_fname.c_str());
          }
        }
        if (mni < threshold) {

          auto ppt = Pattern::get_pattern(g, d.sgl->buf[file_id], d.sgl->ncols);

          auto dd = match(g, { ppt }, false, false, true, threshold, true, true);
          //cout << "num_keys: " << dd.sgl->keys.size() << endl;
          assert(dd.sgl->keys.size() == 1);
          auto it_dd = dd.sgl->keys.begin();
          if (dd.sgl->mni_met.size() > 0 && dd.sgl->mni_met[it_dd->second]) continue;

          auto& distincts_dd = dd.sgl->distinct_vertices[it_dd->second];

          double mni_dd = (double)INT_MAX;
          for (int j = 0; j < distincts_dd.size(); j++) {
            auto& t = distincts_dd[j];
            std::string mni_fname = dd.sgl->DBPath + "/mni_" + std::to_string(it_dd->second) + "_" + std::to_string(j) + ".dat";
            bool file_exist = false;
            if (access(mni_fname.c_str(), F_OK) != -1) {
              std::cout << "wrong in filter" << std::endl;
              ifstream fin(mni_fname);
              int v;
              while (fin >> v) {
                t.insert(v);
              }
              file_exist = true;
            }
            if (mni_dd > t.size()) mni_dd = t.size();
            t.clear();
            if (file_exist) {
              remove(mni_fname.c_str());
            }
          }

          if (mni_dd >= threshold) continue;


          vector<int>().swap(d.sgl->buf[file_id]);
          if (d.sgl->qp_path.size() > 0) {
            d.sgl->qp_path[file_id].clear();
          }
          if (d.sgl->distinct_vertices.size() > 0) {
            d.sgl->distinct_vertices[file_id].clear();
          }

          std::string fname = d.sgl->DBPath + "/" + std::to_string(file_id) + ".dat";
          if (access(fname.c_str(), F_OK) != -1) {
            remove(fname.c_str());
          }
          it = d.sgl->keys.erase(it);
        }
        else {
          ++it;
        }
      }
    }
  }
}  // namespace euler::pattern_mining