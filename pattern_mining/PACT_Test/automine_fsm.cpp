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

  int pat_size = atoi(argv[2]);

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(pat_size));

  cout << "start matching: " << endl;
  util::Timer t;
  t.start();
  auto d2 = match(g, pat2, false, true, true, -1);
  t.stop();

  cout << "Time: " << t.get() << " sec, ";
  cout << "num of frequent patterns: " << d2.sgl->size() << endl;


  return 0;
}
