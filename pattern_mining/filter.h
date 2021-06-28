#pragma once

#include <memory>
#include "db/db.h"
#include "pattern.h"
#include "types.h"

namespace euler::pattern_mining {

void filter(SGList &d, double threshold);
void test_and_filter(const graph::Graph& g, SGList& d, double threshold);

}
