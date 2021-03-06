#include "../gmine.h"

using namespace std;

typedef vector<pair<int, int>> pat_t;

int main(int argc, char *argv[]) {
  // system("rm test_temp/*");

  Graph g;

  g.read_file(argv[1]);

  auto pat7 = PatListing::make_pattern(PatListing().pattern_listing(7));

  cout << "num of unlabeled size-7 patterns: " << pat7.size() << endl;

  cout << "start matchings pat7: " << endl;
  auto d1 = match(g, pat7, false, true, false);
 // cout << "num of matching size-7 patterns: " << d1.sgl->size() << endl;
 // for (auto c: d1.sgl->count) {
  //  cout << c << endl;
 // }


  return 0;
}
