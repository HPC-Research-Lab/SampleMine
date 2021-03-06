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

  vector<double> th = {0.005, 0.01, 0.05};
  vector<double> sp = {20, 40, 80, 160, 320, 640};

  SamplingMethod sm = clustered;

  graph::Graph_CSR_CPU g;

  g.read_graph(argv[1]);

  auto pat3 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(3));
  auto pat2 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(2));

  cout << "num of unlabeled size-3 patterns: " << pat3.size() << endl;
  cout << "num of unlabeled size-2 patterns: " << pat2.size() << endl;

  cout << "start matchings pat3: " << endl;
  util::Timer match_time;
  match_time.start();
  auto d3 = match(g, pat3, true, false, true, true);
  match_time.stop();


  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true, true);

  cout << "match 3 time: " << match_time.get() << " sec" << endl;

  for (auto thh : th) {
    cout << "\n=========threshold: " << thh << "  =============="<< endl;
    int mni_threshold = (int)round(g.num_nodes() * thh);
    filter(d3, mni_threshold);

    auto q3 = get_pattern3(d3);

    cout << "num of size-3 patterns: " << d3.sgl->size() << endl;
    filter(d2, mni_threshold);
    cout << "num of size-2 patterns: " << d2.sgl->size() << endl;

    vector<SGList> sgls = {d3, d3};

    cout << "building tables..." << endl;
    auto H = build_tables(sgls);
    cout << "build table done" << endl;

    for (auto spp : sp) {
      cout << "\nsample param: " << spp << endl;
      {
        cout << "---------no pruning----------" << endl;
        util::Timer t;
        t.start();
        auto d_res = join<true, true, true, 2, 4, 4>(g, H, sgls, false, sm, q3, {spp, spp});
        t.stop();

        cout << "total num of patterns: " << d_res.sgl->size() << endl;
        cout << "num samples: " << d_res.num_samples << endl;
        cout << "estimated total space: " << d_res.total_space << endl;
        cout << "join time: " << t.get() << " sec" << endl;

        filter(d_res, mni_threshold);

        cout << "num frequent patterns: " << d_res.sgl->keys.size() << endl;
      }
    }
  }

  return 0;
}
