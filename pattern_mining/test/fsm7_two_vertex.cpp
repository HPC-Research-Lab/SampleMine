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

  double sup = (size_t)round(thh * g.num_nodes());

  cout << "support threshold: " << sup << endl;

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, true, true, sup);


  filter(d2, sup);
  cout << "num of size-2 frequent patterns: " << d2.sgl->size() << endl;


  cout << "start join for pat3: " << endl;
  vector<SGList> sgls2 = { d2, d2 };

  util::Timer match_time;
  match_time.start();
  auto [H2, subgraph_hist2] = build_tables(sgls2);

  auto [d3, ess3] = join<true, true, true, false, 2, 3, 3>(g, H2, sgls2, true, none, { 0, 0 }, subgraph_hist2, sup, true);

  match_time.stop();

  cout << "join for pat3 time: " << match_time.get() << " sec" << endl;

  filter(d3, sup);

  size_t npat3 = d3.sgl->size();

  cout << "num of size-3 frequent patterns: " << npat3 << endl;

  vector<SGList> sgls = { d3, d3 };

  cout << "building tables..." << endl;
  auto [H3, subgraph_hist3] = build_tables(sgls);
  cout << "build table done" << endl;


  double st_scaled = scale_sampling_param(d2, atof(argv[3]));

  cout << "st_scaled: " << st_scaled << endl;


  util::Timer t;
  t.start();
  auto [d5, ess5] = join<true, true, true, false, 2, 4, 4>(g, H3, sgls, true, none, { 0,0 }, subgraph_hist3, sup, true);


  filter(d5, sup);

  cout << "num of size-5 frequent patterns: " << d5.sgl->size() << endl;

  auto [H5, subgraph_hist5] = build_tables({ d5, d3 });

  cout << "build d5 table done" << endl;
  


  auto [d7, ess7] = join<true, true, true, false, 2, 6, 4>(g, H5, { d5, d3 }, false, none, { st_scaled*st_scaled*st_scaled*st_scaled, st_scaled*st_scaled }, subgraph_hist5, sup, false);

  filter(d7, sup);

  t.stop();

  cout << "Time: " << t.get() << " sec, ";
  if (d7.sgl) {
    cout << "Num patterns: " << d7.sgl->keys.size() << endl;
  }
  else {
    cout << "Num patterns: 0" << endl;
  }



  return 0;
}
