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

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));

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

  double st_scaled = scale_sampling_param(d2, st);

  cout << "scaled sampling param: " << st_scaled << endl;



  vector<SGList> sgls = { d2, d2, d2 };

  cout << "building tables..." << endl;
  auto [H, subgraph_hist] = build_tables(sgls);
  cout << "build table done" << endl;


  double tot_time = 0;
  SGList tot_res;

  for (int i = 0; i < sampling_rounds; i++) {

    util::Timer t;
    t.start();
    auto [d_res, ess] = join<true, true, true, false, 3, 3, 3, 3>(g, H, sgls, false, sm, { st_scaled, st_scaled, st_scaled }, subgraph_hist, sup);
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
