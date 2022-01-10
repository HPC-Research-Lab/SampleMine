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
  double st2 = atof(argv[3]);

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true);

  cout << "start join for pat3: " << endl;
  vector<SGList> sgls2 = { d2, d2 };

   Sampler *sm1;
  if (st1 != 1)
    sm1 = new ProportionalSampler({ st1, st1 });
  else sm1 = &default_sampler;

  util::Timer match_time;
  match_time.start();
  auto H2 = build_tables(sgls2);

  auto [d3, ess3] = join<true, true, false, false, 2, 3, 3>(g, H2, sgls2, true, *sm1, -1, true);

  match_time.stop();

  cout << "join for pat3 time: " << match_time.get() << " sec" << endl << flush;

  size_t npat3 = d3.sgl->size();

  cout << "num of size-3 patterns: " << npat3 << endl;


  vector<SGList> sgls = { d3, d3, d2 };

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  cout << "build table done" << endl;

  Sampler *sm2;
  if (st2 != 1)
    sm2 = new ProportionalSampler({ st2*st2, st2*st2, st2 });
  else sm2 = &default_sampler;

  util::Timer t;
  t.start();
  auto [d_res, ess] = join<true, true, false, false, 3, 4, 4, 3>(g, H, sgls, false, *sm2, -1, false, st2 > 1);
  t.stop();

  cout << "Time: " << t.get() << " sec, ";
  if (d_res.sgl) {
    if (st2 == 1) {
      cout << "Num patterns: " << d_res.sgl->keys.size() << endl;

      vector<double> counts(d_res.sgl->count.begin(), d_res.sgl->count.end());
      sort(counts.begin(), counts.end(), std::greater<double>());

      for (int i = 0; i < (50 > counts.size() ? counts.size() : 50); i++) {
        cout << counts[i] * st1 * st1* st1 * st1 << endl;
      }
    }
    else {
      cout << "Num patterns: " << ess.size() << endl;

      vector<double> counts;

      for (auto& [k, v] : ess) counts.push_back(v);

      sort(counts.begin(), counts.end(), std::greater<double>());

      for (int i = 0; i < (50 > counts.size() ? counts.size() : 50); i++) {
        cout << counts[i] * st1 * st1* st1 * st1 << endl;
      }
    }
  }
  else {
    cout << "Num patterns: 0" << endl;
  }

  return 0;
}
