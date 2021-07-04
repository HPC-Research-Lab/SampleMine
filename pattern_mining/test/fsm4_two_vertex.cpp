#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;


/* Example: ./fsm4_sampling.exe /data/not_backed_up/shared/gpm_data/citeseer.lg 0.001 0 1


 0.001 is the support threshold
 1 is the samping threshold, 0 means no sampling
 1 is the experiments rounds

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

  double st1 = atof(argv[3]);

  double st2 = atof(argv[4]);

  SamplingMethod sm1, sm2;

  if (st1 > 0) {
    sm1 = clustered;
  }
  else {
    sm1 = none;
  }

  if (st2 > 0) {
    sm2 = clustered;
  }
  else {
    sm2 = none;
  }


  cout << "sampling method: " << sm1 << " " << sm2 << endl;

  int sampling_rounds = atoi(argv[5]);

  double sup = (size_t)round(thh * g.num_nodes());

  cout << "support threshold: " << sup << endl;

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, true, true, sup);


  filter(d2, sup);
  cout << "num of size-2 frequent patterns: " << d2.sgl->size() << endl;

  double scaled_st1 = scale_sampling_param(d2, st1);
  double scaled_st2 = scale_sampling_param(d2, st2);

  cout << "start join for pat3: " << endl;
  vector<SGList> sgls2 = { d2, d2 };

  util::Timer match_time;
  match_time.start();
  auto [H2, subgraph_hist2] = build_tables(sgls2);

  auto [d3, ess3] = join<true, true, true, false, 2, 3, 3>(g, H2, sgls2, true, sm1, { scaled_st1, scaled_st1 }, subgraph_hist2, sup, true);

  match_time.stop();

  cout << "join for pat3 time: " << match_time.get() << " sec" << endl;

  filter(d3, sup);

  size_t npat3 = d3.sgl->size();

  cout << "num of size-3 frequent patterns: " << npat3 << endl;


  double st2_scaled = scale_sampling_param(d2, st2);

  cout << "scaled sampling param: " << st2_scaled << endl;


  //  auto pat4 = pattern_mining::PatListing::make_pattern(
    //  pattern_mining::PatListing().pattern_listing(4));

  /*
    cout << "start matchings pat4: " << endl;
    match_time.start();
    auto d4 = match(g, pat4, true, true, true, sup);
    match_time.stop();

    cout << "match 4 time: " << match_time.get() << " sec" << endl;

    filter(d4, sup);

    cout << "num of size-4 frequent patterns: " << d4.sgl->size() << endl;

  */

  vector<SGList> sgls = { d3, d2 };

  cout << "building tables..." << endl;
  auto [H, subgraph_hist] = build_tables(sgls);
  cout << "build table done" << endl;


  double tot_time = 0;
  SGList tot_res;

  for (int i = 0; i < sampling_rounds; i++) {

    util::Timer t;
    t.start();
    auto [d_res, ess] = join<true, true, true, false, 2, 4, 3>(g, H, sgls, false, sm2, { st2_scaled * st2_scaled, st2_scaled }, subgraph_hist, sup);
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
