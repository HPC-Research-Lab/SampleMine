#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
#include "math.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;

int main(int argc, char* argv[]) {
  // system("rm test_temp/*");

  graph::Graph_CSR_CPU g;


  g.read_graph(argv[1]);

  auto pat3 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(3));

  double st = atof(argv[2]);

  cout << "start matchings pat3: " << endl;
  util::Timer match_time;
  match_time.start();
  auto d3 = match(g, pat3, true, false, false);
  match_time.stop();

  cout << "match 3 time: " << match_time.get() << " sec" << endl;

  vector<SGList> sgls = {d3, d3};

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  cout << "build table done" << endl;

  util::Timer t;
  t.start();
  auto [d_res, ess] = join<false, false, false, 2, 4, 4>(g, H, sgls, false, stratified, {st, st});
  t.stop();

  if (d_res.sgl) {
    cout << "total num of patterns: " << d_res.sgl->size() << endl;
    cout << "join time: " << t.get() << " sec" << endl;
    cpp_int total_count = 0;
    for (int i = 0; i < d_res.sgl->count.size(); i++) {
          total_count += d_res.sgl->count[i];
    }

    for (int i = 0; i < d_res.sgl->count.size(); i++) {
         cout << (size_t)round(d_res.sgl->count[i] * ess / total_count) << endl;
    }
  } else {
    cout << "num of patterns: 0" << endl;
  }
  return 0;
}
