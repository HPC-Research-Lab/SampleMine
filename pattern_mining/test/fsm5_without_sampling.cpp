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

auto pat5 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(5));
  auto pat3 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(3));
  auto pat2 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(2));

  double thh = atof(argv[2]);

  size_t sup = (size_t) round(thh * g.num_nodes());

  cout << "support threshold: " << sup << endl;

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, true, true, sup);


    filter(d2, sup);
  cout << "num of size-2 frequent patterns: " << d2.sgl->size() << endl;

  cout << "start matchings pat3: " << endl;
  util::Timer match_time;
  match_time.start();
  auto d3 = match(g, pat3, true, true, true, sup);
  match_time.stop();

  cout << "match 3 time: " << match_time.get() << " sec" << endl;

  filter(d3, sup);

  cout << "num of size-3 frequent patterns: " << d3.sgl->size() << endl;

  cout << "start matchings pat5: " << endl;
  match_time.start();
  auto d5 = match(g, pat5, true, true, true, sup);
  match_time.stop();

  cout << "match 5 time: " << match_time.get() << " sec" << endl;

  filter(d5, sup);

  cout << "num of size-5 frequent patterns: " << d5.sgl->size() << endl;


  vector<SGList> sgls = {d3, d3};

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  cout << "build table done" << endl;

  util::Timer t;
  t.start();
  auto [d_res, ess] = join<true, true, true, false, 2, 4, 4>(g, H, sgls, false, none, {0, 0}, sup);
  t.stop();

  if (d_res.sgl) {
    cout << "total num of patterns: " << d_res.sgl->size() << endl;
    cout << "join time: " << t.get() << " sec" << endl;

    filter(d_res, sup);

    cout << "num frequent patterns: " << d_res.sgl->keys.size() << endl;
  } else {
    cout << "num frequent patterns: 0" << endl;
  }

  return 0;
}
