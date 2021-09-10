#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;


int main(int argc, char* argv[]) {
  // system("rm test_temp/*");

  graph::Graph_CSR_CPU g;

  g.read_graph(argv[1]);

  int pat_size = atoi(argv[2]);

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(pat_size));



  cout << "start sampling: " << endl;

  double timelimit = atof(argv[3]);

  double tot_time = 0.0;

  double rounds = 0;

  SGList res;

  while (tot_time < timelimit) {

    util::Timer t;
    t.start();
    auto d2 = asap::match(g, pat2, false, false, true, -1, false, false, 1);
    t.stop();

    tot_time += t.get();

    if (res.sgl == nullptr) res = d2;
    else res.combine_count(d2);

    rounds += 1;

  }

  cout << "num of patterns: " << res.sgl->size() << endl;
  for (auto& [k, v] : res.sgl->keys) {
    cout << res.sgl->count[v] / rounds << endl;
  }

  return 0;
}
