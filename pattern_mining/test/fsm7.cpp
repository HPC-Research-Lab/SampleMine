#include "gmine.h"

using namespace std;

typedef vector<pair<int, int>> pat_t;

int main(int argc, char *argv[]) {
  // system("rm test_temp/*");

  Graph g;

  g.read_file(argv[1]);

  auto pat4 = PatListing::make_pattern(PatListing().pattern_listing(4));

  cout << "num of unlabeled size-4 patterns: " << pat4.size() << endl;

  cout << "start matchings pat4: " << endl;
  auto d1 =
      match(g, pat4, true,
            false, true);  // for MNI support measure, symmetric need to be stored


  filter(d1, atoi(argv[2]));

  cout << "num of matching size-4 frequent patterns: " << d1->size() << endl;

  auto pat4_fq = Pattern::get_pattern(g, d1, 4);

  auto d2 = match(g, pat4_fq, true, true, false);


  cout << "num of matching size-4 frequent patterns: " << d2->size() << endl;

  auto d3 = join(g, pat4_fq, d2, pat4_fq, d2);
  cout << "num of size-7 frequent patterns: " << d3->size() << endl;

  /*
    vector<shared_ptr<Pattern>> new_pat;
    for (auto &p : d3) {
      p.second.first->print();
      new_pat.push_back(p.second.first);
      print_vec(p.second.second);
    }

    auto d4 = match(g, new_pat, false);*/

  return 0;
}
