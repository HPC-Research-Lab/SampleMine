#pragma once

#include "db/db.h"
#include "pattern.h"
#include <boost/multiprecision/cpp_int.hpp>

namespace euler::pattern_mining {
typedef std::map<int, std::shared_ptr<Pattern>> PatList;
typedef __uint64_t cpp_int;


enum SamplingMethod {stratified, clustered, none};

struct SGList {
  std::shared_ptr<db::MyKV<std::string>> sgl = nullptr;
  PatList patterns;

  SGList() {}
  SGList(const SGList &t)
      : sgl(t.sgl),
        patterns(t.patterns) {}
  SGList(std::shared_ptr<db::MyKV<std::string>> sl, const PatList &pat, size_t ts,
         size_t ns, size_t ninv)
      : sgl(move(sl)),
        patterns(pat) {}

  bool operator<(const SGList &sglist) const {
    return sgl < sglist.sgl;
  }

  void combine(SGList& other, bool mni, bool store) {
    assert(sgl != nullptr);
    sgl->combine(*other.sgl, mni, store);

  }

};

}  // namespace euler::pattern_mining