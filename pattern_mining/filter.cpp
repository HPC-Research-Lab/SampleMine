#include "filter.h"

#include <fcntl.h>
#include <limits.h>
#include <omp.h>
#include <unistd.h>

#include <iostream>
#include <set>
#include <string>
#include <vector>

using namespace std;
namespace euler::pattern_mining {

const static int _Nthreads = atoi(getenv("OMP_NUM_THREADS"));

void filter(SGList &d, int threshold) {
  if (!d.sgl->distinct_vertices.empty()) {
    auto it = d.sgl->keys.begin();
    while (it != d.sgl->keys.end()) {
      size_t file_id = it->second;
      if (d.sgl->mni_met.size() > 0 && d.sgl->mni_met[file_id]) {
        it++;
        continue;
      }
      auto &distincts = d.sgl->distinct_vertices[file_id];

      int mni = INT_MAX;
      for (int j=0; j<distincts.size(); j++) {
        auto &t = distincts[j];
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
        std::string fname = d.sgl->DBPath + "/" + std::to_string(file_id) + ".dat";
        if (access(fname.c_str(), F_OK) != -1) {
          remove(fname.c_str());
        }
        it = d.sgl->keys.erase(it);
      } else {
        ++it;
      }
    }
  }
}
/*
void filter(SGList &d, int threshold) {
  auto it = d.sgl->keys.begin();
  while (it != d.sgl->keys.end()) {
    auto buf = d.sgl->getbuf(get<0>(it->second), d.sgl->ncols);
    // int pat_id = buf[0];
    vector<vector<set<int>>> distincts_vec(_Nthreads);
    for (int i=0; i<distincts_vec.size(); i++) distincts_vec[i].resize(d.sgl->ncols - 1);
    auto it_buf = buf.begin();

    while (true) {
      #pragma omp parallel for num_threads(_Nthreads)
      for (int j = 0; j < it_buf.buffer_size; j += d.sgl->ncols) {
        int tid = omp_get_thread_num();
        for (int i = 0; i < d.sgl->ncols - 1; i++) {
          distincts_vec[tid][i].insert(*(it_buf.buffer + j + 1 + i));
        }
      }
      if (!it_buf.has_next) break;
      it_buf.next();
    }


    vector<set<int>> distincts(d.sgl->ncols - 1);

    for (int i=0; i<d.sgl->ncols - 1; i++) {
      for (int j=0; j<_Nthreads; j++) {
        distincts[i].insert(distincts_vec[j][i].begin(), distincts_vec[j][i].end());
      }
    }

    int mni = INT_MAX;
    for (const auto &t : distincts) {
      if (mni > t.size()) mni = t.size();
    }
    if (mni < threshold) {
      size_t idx = get<0>(it->second);
      vector<int>().swap(d.sgl->buf[idx]);
      std::string fname = d.sgl->DBPath + "/" + std::to_string(idx) + ".dat";
      if (access(fname.c_str(), F_OK) != -1) {
        remove(fname.c_str());
      }
      it = d.sgl->keys.erase(it);
    } else {
      ++it;
    }
  }
}*/
}  // namespace euler::pattern_mining