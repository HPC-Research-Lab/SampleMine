#include "util.h"

namespace euler::util {

void print_vec(const std::vector<int> &a) {
  for (int i : a) std::cout << i << " ";
  std::cout << std::endl;
}


}  // namespace euler::util