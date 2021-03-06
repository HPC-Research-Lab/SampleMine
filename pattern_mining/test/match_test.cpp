#include "../gmine.h"

using namespace std;
using namespace euler;
using namespace euler::graph;
using namespace euler::pattern_mining;

int main(int argc, char *argv[]) {
  graph::Graph_CSR_CPU g;

  g.read_graph(argv[1]);

  auto pat3 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(3));
  auto pat4 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(4));
  auto pat5 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(5));
  auto pat2 = pattern_mining::PatListing::make_pattern(
      pattern_mining::PatListing().pattern_listing(2));

  cout << "num of unlabeled size-3 patterns: " << pat3.size() << endl;
  cout << "num of unlabeled size-2 patterns: " << pat2.size() << endl;
  cout << "num of unlabeled size-4 patterns: " << pat4.size() << endl;

  cout << "start matchings pat3: " << endl;
  //auto d1 = match(g, pat3, true);

  cout << "start matchings pat2: " << endl;
  //auto d3 = match(g, pat2, true);

  cout << "start matchings pat4: " << endl;
  //auto d4 = match(g, pat4, true);

  cout << "start matchings pat5: " << endl;
  // auto d4 = match(g, pat5, true);

  auto d5 = match(g, pattern_mining::PatListing::make_pattern({{{0,1}, {1,2}, {2,3}, {3,4}}}), true);
}
