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

  auto pat3 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(3));
  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));

  auto pat4 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(4));

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, true, true);

  //d2.print();

  cout << "start matchings pat3: " << endl;
  util::Timer match_time;
  match_time.start();
  auto d3 = match(g, pat3, true, true, true);
  match_time.stop();


  cout << "match 3 time: " << match_time.get() << " sec" << endl;

  cout << "start matchings pat4: " << endl;
  match_time.start();
  auto d4 = match(g, pat4, true, true, true);
  match_time.stop();


  cout << "match 4 time: " << match_time.get() << " sec" << endl;

  d4.print_counts();


  vector<SGList> sgls = { d3, d2 };

  cout << "building tables..." << endl;
  auto [H, sw] = build_tables(sgls);
  cout << "build table done" << endl;

  util::Timer t;
  t.start();
  auto [d_res, ess] = join<true, true, false, false, 2, 4, 3>(g, H, sgls, false, none, { 0, 0 });
  t.stop();


  cout << "join time: " << t.get() << " sec" << endl;

  d_res.print_counts();
  return 0;
}
