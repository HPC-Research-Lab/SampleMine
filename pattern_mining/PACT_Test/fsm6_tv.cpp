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

  double thh = atof(argv[2]);

  double st2 = atof(argv[3]);

  double sup = (size_t)round(thh * g.num_nodes());


  cout << "support threshold: " << sup << endl;


  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true, sup);


  filter(d2, sup);
  cout << "num of size-2 frequent patterns: " << d2.sgl->size() << endl;

  double scaled_st2 = scale_sampling_param(d2, st2);

  cout << "start join for pat3: " << endl;
  vector<SGList> sgls2 = { d2, d2 };

  util::Timer match_time;
  match_time.start();
  auto H2 = build_tables(sgls2);

  auto [d3, ess3] = join<true, true, true, true, 2, 3, 3>(g, H2, sgls2, true, default_sampler, sup, true);

  match_time.stop();

  cout << "join for pat3 time: " << match_time.get() << " sec" << endl;

  filter(d3, sup);

  size_t npat3 = d3.sgl->size();

  cout << "num of size-3 frequent patterns: " << npat3 << endl;


  double st2_scaled = scale_sampling_param(d2, st2);

  cout << "scaled sampling param: " << st2_scaled << endl;


  vector<SGList> sgls = { d3, d3, d2 };

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  auto  subgraph_hist = get_subgraph_hist(sgls);
  cout << "build table done" << endl;

  Sampler *sm2;
  if (st2 > 0)
    sm2 = new BudgetSampler(subgraph_hist, {st2_scaled*st2_scaled, st2_scaled*st2_scaled, st2_scaled});
  else sm2 = &default_sampler;

  util::Timer t;
  t.start();
  auto [d_res, ess] = join<true, true, true, true, 3, 4, 4, 3>(g, H, sgls, false, *sm2, sup, false, false);
  t.stop();

  filter(d_res, sup);

  cout << "Time: " << t.get() << " sec, ";
  if (d_res.sgl) {
    cout << "Num patterns: " << d_res.sgl->keys.size() << endl;
  }
  else {
    cout << "Num patterns: 0" << endl;
  }

  return 0;
}
