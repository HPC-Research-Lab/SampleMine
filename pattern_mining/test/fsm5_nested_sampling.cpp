#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;


/*

Nested sampling is used when thre are too many size-3 subgraphs and cannot be store in memory

 Example: ./fsm5_nested_sampling.exe /data/not_backed_up/shared/gpm_data/citeseer.lg 0.001 8 1 10

 0.001 is the support threshold
 8 is the sampling threshold for obtaining size-3 subgraphs
 1 is the sampling threshold for obtaining size-5 subgraphs
 10 is the sampling rounds


*/


int main(int argc, char* argv[]) {
  // system("rm test_temp/*");

  graph::Graph_CSR_CPU g;


  g.read_graph(argv[1]);

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));

  double thh = atof(argv[2]);

  double st1 = atof(argv[3]);
  double st2 = atof(argv[4]);

  int sampling_rounds = atoi(argv[5]);

  double sup = (size_t)round(thh * g.num_nodes());

  cout << "support threshold: " << sup << endl;

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, true, true, sup);

  filter(d2, sup);
  cout << "num of size-2 frequent patterns: " << d2.sgl->size() << endl;


  double tot_time = 0;
  SGList tot_res;

  cout << "start join for pat3: " << endl;
  vector<SGList> sgls2 = { d2, d2 };

  for (int i = 0; i < sampling_rounds; i++) {
    util::Timer match_time;
    match_time.start();
    auto [H2, sw2] = build_tables(sgls2);

    auto [d3, ess3] = join<true, true, true, false, 2, 3, 3>(g, H2, sgls2, true, clustered, { st1, st1 }, sw2, sup, true);

    match_time.stop();

    cout << "join for pat3 time: " << match_time.get() << " sec" << endl;

    tot_time += match_time.get();

    filter(d3, sup);

    cout << "num of size-3 frequent patterns: " << d3.sgl->size() << endl;

    vector<SGList> sgls = { d3, d3 };

    cout << "building tables..." << endl;
    auto [H, sw] = build_tables(sgls);
    cout << "build table done" << endl;

    util::Timer t;
    t.start();
    auto [d_res, ess] = join<true, true, true, false, 2, 4, 4>(g, H, sgls, false, clustered, { st2, st2 }, sw, sup);
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
