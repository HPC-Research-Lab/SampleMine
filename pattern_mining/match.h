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

static std::vector<size_t> dummy = {};


SGList match(
    const graph::Graph &g, const PatList &patterns,
    bool store_data = true, bool edge_induced=false, bool has_label=false, size_t mni=0, double sampling_threshold=0);

//std::pair<std::unordered_set<unsigned long>, std::unordered_set<unsigned long>> get_pattern3(const SGList &sgl);

}
