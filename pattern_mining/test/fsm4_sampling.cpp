#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;


/* Example: ./fsm4_sampling.exe /data/not_backed_up/shared/gpm_data/citeseer.lg 0.001 1 10


 0.001 is the support threshold
 1 is the samping threshold, set it to 1 or 2 for experiment
 10 is the sampling rounds, let's use 10

*/

int main(int argc, char* argv[]) {
  // system("rm test_temp/*");

  graph::Graph_CSR_CPU g;


  g.read_graph(argv[1]);

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));
  auto pat3 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(3));

  double thh = atof(argv[2]);

  double st = atof(argv[3]);

  SamplingMethod sm;

  if (st > 0) {
    sm = clustered;
  }
  else {
    sm = none;
  }

  cout << "sampling method: " << sm << endl;

  int sampling_rounds = atoi(argv[4]);

  double sup = (size_t)round(thh * g.num_nodes());

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


  vector<SGList> sgls = { d2, d3 };

  cout << "building tables..." << endl;
  auto [H, sw] = build_tables(sgls);
  cout << "build table done" << endl;


  double tot_time = 0;
  SGList tot_res;

  for (int i = 0; i < sampling_rounds; i++) {

    util::Timer t;
    t.start();
    auto [d_res, ess] = join<true, true, true, false, 2, 3, 4>(g, H, sgls, false, sm, { st, st }, sw, sup);
    t.stop();

    filter(d_res, sup);

    if (tot_res.sgl == nullptr) tot_res = d_res;
    else
      tot_res.combine(d_res, true, false);

    tot_time += t.get();

    cout << "Time: " << tot_time << " sec, ";
    if (tot_res.sgl) {
      cout << "Num patterns: " << tot_res.sgl->keys.size() << endl;
    }
    else {
      cout << "Num patterns: 0" << endl;
    }

  }


  return 0;
}

