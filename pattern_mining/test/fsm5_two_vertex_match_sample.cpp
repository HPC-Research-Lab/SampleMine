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

  auto pat3 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(3));

  double thh = atof(argv[2]);

  double st1 = atof(argv[3]);

  double sup = (size_t)round(thh * g.num_nodes());


  cout << "support threshold: " << sup << endl;


  cout << "start matchings pat3: " << endl;
  auto d3 = match(g, pat3, true, false, true, sup, false, false, st1);

  filter(d3, sup);

  size_t npat3 = d3.sgl->size();

  cout << "num of size-3 frequent patterns: " << npat3 << endl;

  vector<SGList> sgls = { d3, d3 };

  cout << "building tables..." << endl;
  auto [H, subgraph_hist] = build_tables(sgls);
  cout << "build table done" << endl;

  d3.sgl->clear();

  double tot_time = 0;
  SGList tot_res;

  // auto original_table_size = get_table_size(subgraph_hist);

  vector<double> vst = {2,4,6,8 };

  for (double st2 : vst) {

    cout << "st2: " << st2 << endl;

    SamplingMethod sm2;


    if (st2 > 0) {
      sm2 = clustered;
    }
    else {
      sm2 = none;
    }

    util::Timer t;
    t.start();
    auto [d_res, ess] = join<true, true, true, false, 2, 4, 4>(g, H, sgls, false, sm2, { st2, st2 }, subgraph_hist, sup, false);
    t.stop();

    double ntot = d_res.sgl->tot_count();

    filter(d_res, sup);


    cout << "Time: " << t.get() << " sec, ";
    if (d_res.sgl) {
      cout << "Num patterns: " << d_res.sgl->keys.size() << endl;
    }
    else {
      cout << "Num patterns: 0" << endl;
    }

  }


  return 0;
}
