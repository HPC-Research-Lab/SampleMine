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

  double st = atof(argv[2]);

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true);

  vector<SGList> sgls = { d2, d2, d2, d2 };

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  cout << "build table done" << endl;

  std::map<std::string, double> ess;
  double tot_time = 0;
  int num_experiment = atoi(argv[3]);

  for (int iter = 0; iter < num_experiment; iter++) {

    util::Timer t;
    t.start();
    auto [d_res, est] = join<true, false, false, true, 4, 3, 3, 3, 3>(g, H, sgls, true, clustered, { 20, 10, 10, 5});
    t.stop();
    tot_time += t.get();
    for (auto &[k,v]: est) {
      if (ess.find(k) == ess.end()) ess[k] = 0;
      ess[k] += v;
    }
  }


  vector<size_t> counts;

  cout << "total num of patterns: " << ess.size() << endl;
  cout << "join time: " << tot_time << " sec" << endl;

  for (auto& [k, v] : ess) {
    size_t t = (size_t)round(v / num_experiment);
    if (t == 0) t = 1;
    counts.push_back(t);
  }

  sort(counts.begin(), counts.end(), greater<size_t>());

  for (int i = 0; i < counts.size(); i++) {
    cout << counts[i] << endl;
  }

  return 0;
}
