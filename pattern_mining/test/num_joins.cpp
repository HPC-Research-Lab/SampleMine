#include "../gmine.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
typedef boost::multiprecision::cpp_dec_float_100 value_type;


using namespace std;

typedef vector<pair<int, int>> pat_t;

int main(int argc, char *argv[]) {
  // system("rm test_temp/*");

  Graph g;

  g.read_file(argv[1]);

  auto pat3 = PatListing::make_pattern(PatListing().pattern_listing(2));

  cout << "num of unlabeled size-3 patterns: " << pat3.size() << endl;

  cout << "start matchings pat3: " << endl;
  auto d1 = match(g, pat3, true, true);

  auto d2 = join(g, {d1, d1, d1, d1, d1, d1}, stratified, atoi(argv[2]));

  cout << "num of size-5 patterns: " << d2.sgl->size() << endl;
  for (auto c: d2.sgl->count) {
    value_type p = value_type(c) / value_type(d2.num_samples);
    value_type z = 1.96 * sqrt(p*(1-p))/sqrt(d2.num_samples);
    value_type p_lower = p - z > 0 ? p - z : value_type(0);
    value_type p_upper = p + z;
    cout << c << "\t" << size_t(p_lower * d2.total_space) << "\t" << size_t(p_upper * d2.total_space) << endl;
  }
  cout << "computation ratio: " << d2.total_space / d2.num_samples << endl;

  return 0;
}
