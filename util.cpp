#include "util.h"
#include <random>  // mt19937 and uniform_int_distribution


namespace euler::util {

static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<> dis(0, 1);

  double random_number() {
    return dis(gen);
  }


}  // namespace euler::util