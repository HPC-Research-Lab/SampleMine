#include "../gmine.h"

using namespace std;

typedef vector<pair<int, int>> pat_t;

int main(int argc, char *argv[]) {
  // system("rm test_temp/*");

  Graph g;

  g.read_file(argv[1]);

  auto pat5 = PatListing::make_pattern(PatListing().pattern_listing(4));

  cout << "num of unlabeled size-4 patterns: " << pat5.size() << endl;

  cout << "start matchings pat4: " << endl;
  auto d1 = match(g, pat5, false, true, false);


  return 0;
}
