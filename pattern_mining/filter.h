#pragma once

#include <memory>
#include "db/db.h"
#include "pattern.h"
#include "types.h"

namespace euler::pattern_mining {

void filter(SGList &d, double threshold);

}
