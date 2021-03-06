#include "pattern_mining/gmine.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
typedef boost::multiprecision::cpp_dec_float_100 value_type;


using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;

int main(int argc, char *argv[]) {
  // system("rm test_temp/*");

  graph::Graph_CSR_CPU g;

  g.read_graph(argv[1]);

  auto pat3 = pattern_mining::PatListing::make_pattern(pattern_mining::PatListing().pattern_listing(3));
  auto pat2 = pattern_mining::PatListing::make_pattern(pattern_mining::PatListing().pattern_listing(2));

  cout << "num of unlabeled size-3 patterns: " << pat3.size() << endl;
  cout << "num of unlabeled size-2 patterns: " << pat2.size() << endl;

  cout << "start matchings pat3: " << endl;
  auto d1 = match(g, pat3, true);

  cout << "start matchings pat2: " << endl;
  auto d3 = match(g, pat2, true);

  auto d2 = join(g, {d3, d3, d3}, none);

  cout << "num of size-4 patterns: " << d2.sgl->size() << endl;
  for (auto c: d2.sgl->count) {
    /*value_type p = value_type(c) / value_type(d2.num_samples);
    value_type z = 1.96 * sqrt(p*(1-p))/sqrt(d2.num_samples);
    value_type p_lower = p - z > 0 ? p - z : value_type(0);
    value_type p_upper = p + z;
    cout << c << "\t" << size_t(p_lower * d2.total_space) << "\t" << size_t(p_upper * d2.total_space) << endl;*/
    cout << c << endl;
  }
  //cout << "computation ratio: " << d2.total_space / d2.num_samples << endl;

  return 0;
}
