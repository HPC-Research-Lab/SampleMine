#pragma once

#include "db/db.h"
#include "pattern.h"
#include <boost/multiprecision/cpp_int.hpp>

namespace euler::pattern_mining {
typedef std::vector<std::shared_ptr<Pattern>> PatList;
typedef __uint64_t cpp_int;


enum SamplingMethod {stratified, clustered, none};

struct SGList {
  std::shared_ptr<db::MyKV<std::string>> sgl = nullptr;
  PatList unlabeled_patterns;

  SGList() {}
  SGList(const SGList &t)
      : sgl(t.sgl), unlabeled_patterns(t.unlabeled_patterns) {}
  SGList(std::shared_ptr<db::MyKV<std::string>> sl, const PatList &up)
      : sgl(move(sl)), unlabeled_patterns(up) {}

  bool operator<(const SGList &sglist) const {
    return sgl < sglist.sgl;
  }

  void combine(SGList& other, bool mni, bool store) {
    assert(sgl != nullptr);
    sgl->combine(*other.sgl, mni, store);
  }

  void print() {
    sgl->print();
  }

};

}  // namespace euler::pattern_mining