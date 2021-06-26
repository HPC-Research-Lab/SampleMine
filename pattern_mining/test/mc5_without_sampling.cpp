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
  auto d2 = match(g, pat2, true, true, true);

  vector<SGList> sgls = {d2, d2, d2, d2};

  cout << "building tables..." << endl;
  auto [H, sw] = build_tables(sgls);
  cout << "build table done" << endl;

  util::Timer t;
  t.start();
  auto [d_res, ess]  = join<true, true, false, false, 4, 3, 3, 3, 3>(g, H, sgls, false, none, {0, 0, 0, 0});
  t.stop();


  vector<size_t> counts;

  if (d_res.sgl) {
    cout << "total num of patterns: " << d_res.sgl->size() << endl;
    cout << "join time: " << t.get() << " sec" << endl;

    for (size_t v: d_res.sgl->count) {
          counts.push_back(v);
    }

    sort(counts.begin(), counts.end(), greater<size_t>());

    for (int i=0; i < counts.size(); i++) {
      cout << counts[i] << endl;
    }

   // d_res.print();
  } else {
    cout << "num of patterns: 0" << endl;
  }

  return 0;
}
