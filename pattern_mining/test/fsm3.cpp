#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;


int main(int argc, char *argv[]) {

 graph::Graph_CSR_CPU g;

  g.read_graph(argv[1]);

  auto pat3 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(atoi(argv[2])));

  cout << "num of unlabeled size-3 patterns: " << pat3.size() << endl;

  cout << "start matchings pat3: " << endl;
  auto d1 =
      match(g, pat3, true,
            false, true); 


  size_t sum = 0;
  for (size_t t: d1.sgl->count) {
    sum += t;
  }

  cout << "num patterns: " << d1.sgl->keys.size() << endl;

  filter(d1, g.num_nodes() / 100);

  cout << "num frequent patterns: " << d1.sgl->keys.size() << endl; 

  vector<size_t> ss;
  for (auto &s: d1.sgl->keys) {
    ss.push_back(d1.sgl->count[get<0>(s.second)]);
  }

  sort(ss.begin(), ss.end());

  size_t partial = 0;
  for (size_t t: ss) {
    partial += t;
  }




 cout << float(partial) / sum << endl;
  
  return 0;

}
