#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;

int main(int argc, char *argv[]) {
  // system("rm test_temp/*");

  graph::Graph_CSR_CPU g;

  g.read_graph(argv[1]);

  auto pat3 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(3));
  auto pat2 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(2));

  cout << "num of unlabeled size-3 patterns: " << pat3.size() << endl;
  cout << "num of unlabeled size-2 patterns: " << pat2.size() << endl;

  cout << "start matchings pat3: " << endl;
  util::Timer match_time;
  match_time.start();
  auto d3 = match(g, pat3, true);
  match_time.stop();

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true);

  SamplingMethod sm = none;

  util::Timer t;
  t.start();
  auto d_res = join(g, {d3, d2}, false, false, false, sm, {0.03125, 0.125});
  t.stop();

  cout << "total num of patterns: " << d_res.sgl->size() << endl;
  cout << "num samples: " << d_res.num_samples << endl;
  cout << "estimated total space: " << d_res.total_space << endl;
  cout << "match 3 time: " << match_time.get() << " sec" << endl;
  cout << "join time: " << t.get() << " sec" << endl;
  
  for (int i = 0; i < d_res.sgl->count.size(); i++) {
    if (sm == none) {
      cout << d_res.sgl->count[i] << endl;
    } else {
      value_type p =
          value_type(d_res.sgl->count[i]) / value_type(d_res.num_samples);
      value_type z = 1.96 * sqrt(p * (1 - p)) / sqrt((double)d_res.num_samples);
      value_type p_lower = p - z > 0 ? p - z : value_type(0);
      value_type p_upper = p + z;
      cout << size_t(d_res.sgl->count[i] / value_type(d_res.num_samples) *
                     d_res.total_space)
           << endl;
    }
  }
  return 0;
}
