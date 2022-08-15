#pragma once

#include <list>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include "db/db.h"
#include "graph/graph.h"
#include "pattern.h"
#include "types.h"

namespace euler::pattern_mining {

//static std::vector<size_t> dummy_match = {};


SGList match(
    const graph::Graph &g, const PatList &patterns,
    bool store_data = true, bool edge_induced=false, bool output_labeled=false,  double mni=-1, bool testing = false, bool pattern_labeled = false, double sampling_threshold=0);

std::pair<std::unordered_set<unsigned long>, std::unordered_set<unsigned long>> get_pattern3(const SGList &sgl);

}

namespace euler::pattern_mining::asap {

//static std::vector<size_t> dummy_match = {};

SGList match(
    const graph::Graph &g, const PatList &patterns,
    bool store_data = true, bool edge_induced=false, bool output_labeled=false,  double mni=-1, bool testing = false, bool pattern_labeled = false, double sampling_threshold=0);


}

namespace euler::pattern_mining::old_match {

//static std::vector<size_t> dummy_match = {};

SGList match(
    const graph::Graph &g, const PatList &patterns,
    bool store_data = true, bool edge_induced=false, bool output_labeled=false,  double mni=-1, bool testing = false, bool pattern_labeled = false, double sampling_threshold=0);


}
