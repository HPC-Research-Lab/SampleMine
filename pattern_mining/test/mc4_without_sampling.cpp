#include <boost/multiprecision/cpp_dec_float.hpp>
#include <algorithm>

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

  auto pat2 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(2));

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true);

  vector<SGList> sgls = {d2, d2, d2};

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  cout << "build table done" << endl;

  util::Timer t;
  t.start();
  auto [d_res, ess] = join<true, false, false, 3, 3, 3, 3>(g, H, sgls, false, none, {0, 0, 0});
  t.stop();


  vector<int> counts;

  if (d_res.sgl) {
    cout << "total num of patterns: " << d_res.sgl->size() << endl;
    cout << "join time: " << t.get() << " sec" << endl;

    for (int i = 0; i < d_res.sgl->count.size(); i++) {
          counts.push_back(d_res.sgl->count[i]);
    }

    sort(counts.begin(), counts.end(), greater<int>());

    for (int i=0; i < counts.size(); i++) {
      cout << counts[i] << endl;
    }
  } else {
    cout << "num of patterns: 0" << endl;
  }

  return 0;
}
