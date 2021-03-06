#include "../gmine.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
typedef boost::multiprecision::cpp_dec_float_100 value_type;


using namespace std;

typedef vector<pair<int, int>> pat_t;

int main(int argc, char *argv[]) {
  // system("rm test_temp/*");

  Graph g;

  g.read_file(argv[1]);

  auto pat3 = PatListing::make_pattern(PatListing().pattern_listing(3));

  cout << "num of unlabeled size-3 patterns: " << pat3.size() << endl;

  cout << "start matchings pat3: " << endl;
  auto d1 = match(g, pat3, true, true);

  auto pat2 = PatListing::make_pattern(PatListing().pattern_listing(2));

  cout << "num of unlabeled size-2 patterns: " << pat2.size() << endl;

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, true);


  auto d4 = join(g, {d2, d1, d1}, stratified, atoi(argv[2]));

  cout << "num of size-6 patterns: " << d4.sgl->size() << endl;
  for (auto c: d4.sgl->count) {
    value_type p = value_type(c) / value_type(d4.num_samples);
    value_type z = 1.96 * sqrt(p*(1-p))/sqrt(d4.num_samples);
    value_type p_lower = p - z > 0 ? p - z : value_type(0);
    value_type p_upper = p + z;
    cout << c << "\t" << size_t(p_lower * d4.total_space) << "\t" << size_t(p_upper * d4.total_space) << endl;
  }
  cout << "computation ratio: " << d4.total_space / d4.num_samples << endl;

  return 0;
}
