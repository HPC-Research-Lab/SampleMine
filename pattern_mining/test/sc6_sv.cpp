#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;


// citeseer 4,5,6,7  // sampling ratio: 2,3,4
// mico 4
int main(int argc, char* argv[]) {
  // system("rm test_temp/*");

  graph::Graph_CSR_CPU g;


  g.read_graph(argv[1]);

  cout << "max degree: " << g.max_degree() << endl;

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));

  double st1 = atof(argv[2]);

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true);

  cout << "start join for pat6: " << endl;
  vector<SGList> sgls2 = { d2, d2, d2, d2, d2};

  Sampler *sm1;
  if (st1 != 1)
    sm1 = new ProportionalSampler2({ st1, st1, st1, st1, st1 });
  else sm1 = &default_sampler;

  util::Timer match_time;
  match_time.start();
  auto H2 = build_tables(sgls2);

  auto [d, ess] = join<true, true, false, false, 5, 3, 3, 3, 3, 3>(g, H2, sgls2, true, *sm1, -1, true);

  match_time.stop();

  cout << "join for pat6 time: " << match_time.get() << " sec" << endl << flush;

  size_t npat = d.sgl->size();

  cout << "num of  patterns: " << npat << endl;


 
  return 0;
}
