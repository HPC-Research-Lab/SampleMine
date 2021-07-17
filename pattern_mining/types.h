#pragma once

#include "db/db.h"
#include "pattern.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
#include <algorithm>

namespace euler::pattern_mining {
  typedef std::vector<std::shared_ptr<Pattern>> PatList;
  typedef __uint64_t cpp_int;


  enum SamplingMethod { stratified, clustered, none };

  struct SGList {
    std::shared_ptr<db::MyKV<std::string>> sgl = nullptr;
    PatList patterns;

    SGList() {}
    SGList(const SGList& t)
      : sgl(t.sgl), patterns(t.patterns) {}
    SGList(std::shared_ptr<db::MyKV<std::string>> sl, const PatList& up)
      : sgl(move(sl)), patterns(up) {}

    bool operator<(const SGList& sglist) const {
      return sgl < sglist.sgl;
    }

    void combine(SGList& other, bool mni, bool store, bool adaptive_sampling=false) {
      assert(sgl != nullptr);
      sgl->combine(*other.sgl, mni, store, adaptive_sampling);
    }

    void print() {
      sgl->print();
    }

    void get_quick_pattern_path(const std::vector<std::vector<std::array<int, 4>>> &qp_count) {
      sgl->get_pattern_path(qp_count);
    }

    void print_counts() {

      if (sgl != nullptr) {
        std::cout << "total num of patterns: " << sgl->size() << std::endl;
        std::vector<size_t> counts;
        for (int i = 0; i < sgl->count.size(); i++) {
          counts.push_back(sgl->count[i]);
        }
        sort(counts.begin(), counts.end(), std::greater<int>());
        for (int i = 0; i < counts.size(); i++) {
          std::cout << counts[i] << std::endl;
        }
      }
      else {
        std::cout << "total num of patterns: 0" << std::endl;

      }
    }

  };

}  // namespace euler::pattern_mining