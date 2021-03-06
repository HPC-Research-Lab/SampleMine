#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;

int main(int argc, char *argv[]) {
  graph::Graph_CSR_CPU g;

  g.read_graph(argv[1]);

  auto pat = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(atoi(argv[2])));

  cout << "num of unlabeled size-4 patterns: " << pat.size() << endl;

  cout << "start matchings pat: " << endl;
  auto d3 = match(g, pat, false, true, false);

  return 0;
}
