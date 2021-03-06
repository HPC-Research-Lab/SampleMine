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

  cout << "support threshold: " << thh << endl;

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true, true);


  int mni_threshold = (int)round(g.num_nodes() * thh);

  filter(d2, mni_threshold);
  cout << "num of size-2 patterns: " << d2.sgl->size() << endl;

  vector<SGList> sgls = {d2, d2, d2, d2};

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  cout << "build table done" << endl;

  util::Timer t;
  t.start();
  auto d_res = join<true, true, true, 4, 3, 3, 3, 3>(g, H, sgls, false, none);
  t.stop();

  if (d_res.sgl) {
    cout << "total num of patterns: " << d_res.sgl->size() << endl;
    cout << "join time: " << t.get() << " sec" << endl;

    filter(d_res, mni_threshold);

    cout << "num frequent patterns: " << d_res.sgl->keys.size() << endl;
  } else {
    cout << "num frequent patterns: 0" << endl;
  }

  return 0;
}
