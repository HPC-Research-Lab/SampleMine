#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;

//     find the frequent subgraphs with a label 1 OR a label 2
//    ./q3.exe /data/not_backed_up/shared/gpm_data/citeseer.lg 0.001 0
//    ./q3.exe /data/not_backed_up/shared/gpm_data/citeseer.lg 0.005 0
// sampling threshold 4,6,8,10


// citeseer 4,5,6,7
// mico 4
// size-5
// 3 + 3 (step=1)
// 3 + 3 + 3 (step=2)
int my_query(const graph::Graph& g, const int *s, std::shared_ptr<Pattern> pat2, int step) {
  if (step == 1) {
    int n1 = 0;
    int n2 = 0;
    for (int i=1; i<6; i++) {
      int l = g.get_vertex_label(s[i]);
      if (l == 1) n1++;
      if (l == 2) n2++;
    }
    if (n1 == 0 && n2 == 0) return -1;
  }
  return 0;
}

int main(int argc, char* argv[]) {
  // system("rm test_temp/*");

  graph::Graph_CSR_CPU g;


  g.read_graph(argv[1]);

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));

  double thh = atof(argv[2]);

  double st2 = atof(argv[3]);

  SamplingMethod sm2;

  if (st2 > 0) {
    sm2 = clustered;
  }
  else {
    sm2 = none;
  }

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
  auto [H2, subgraph_hist2] = build_tables(sgls2);

  auto [d3, ess3] = join<true, true, true, true, 2, 3, 3>(g, H2, sgls2, true, none, { 0, 0 }, subgraph_hist2, sup, true);

  match_time.stop();

  cout << "join for pat3 time: " << match_time.get() << " sec" << endl;

  filter(d3, sup);

  size_t npat3 = d3.sgl->size();

  cout << "num of size-3 frequent patterns: " << npat3 << endl;


  double st2_scaled = scale_sampling_param(d2, st2);

  cout << "scaled sampling param: " << st2_scaled << endl;


  vector<SGList> sgls = { d3, d3 };

  cout << "building tables..." << endl;
  auto [H, subgraph_hist] = build_tables(sgls);
  cout << "build table done" << endl;

  util::Timer t;
  t.start();
  auto [d_res, ess] = join<true, true, true, true, 2, 4, 4>(g, H, sgls, false, sm2, { st2_scaled * st2_scaled, st2_scaled * st2_scaled }, subgraph_hist, sup, false, false, false, join_dummy1, my_query);
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
