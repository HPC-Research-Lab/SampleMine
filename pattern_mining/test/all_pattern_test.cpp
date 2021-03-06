#include <iostream>

#include "../gmine.h"

using namespace std;

int main(int argc, char *argv[]) {
  auto t =
      PatListing::make_pattern(PatListing().pattern_listing(atoi(argv[1])));

  int sep = atoi(argv[2]);

  cout << t.size() << endl;

  while (!t.empty()) {
    auto pats = t;
    t.clear();
    int c = 0;
    for (auto &p : pats) {
      if (!p.second->is_separable(sep)) {
        c++;
        t.insert(p);
      }
    }
    cout << sep << " " << c << endl;
  }

  return 0;
}