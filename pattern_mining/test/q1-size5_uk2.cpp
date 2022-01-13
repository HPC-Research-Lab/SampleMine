#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;


// p x x x x x


// 3 + 3 = 5
// 3 + 2 = 4
// 3 + 3 + 2 = 6
// 3 + 3 + 3 = 7

// citeseer-4,5,6,7
// mico-4
// find the size-5 subgraphs with a label 1 AND a label 2
class MyQuery: public Query {
int operator()(const graph::Graph& g, util::span<const int> s, std::shared_ptr<Pattern> pat, int step) {
  if (step == 1) {
    int n1 = 0;
    int n2 = 0;
    for (int i=1; i<6; i++) {
      int l = g.get_vertex_label(s[i]);
      if (l == 1) n1++;
      if (l == 2) n2++;
    }
    if (n1 == 0 || n2 == 0) return -1;
  }
  return 0;
}
};


int main(int argc, char* argv[]) {
  // system("rm test_temp/*");

  graph::Graph_CSR_CPU g;


  g.read_graph(argv[1]);

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));


  double st2 = atof(argv[2]);

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true);

  cout << "start join for pat3: " << endl;
  vector<SGList> sgls2 = { d2, d2 };

  util::Timer match_time;
  match_time.start();
  auto H2 = build_tables(sgls2);

//    Sampler *sm_tmp = new ProportionalSampler({ 16, 16 });
    Sampler *sm_tmp = new ProportionalSampler({ 32, 32 });

  auto [d3, ess3] = join<true, true, false, false, 2, 3, 3>(g, H2, sgls2, true, *sm_tmp, -1, true);

  match_time.stop();

  cout << "join for pat3 time: " << match_time.get() << " sec" << endl;

  size_t npat3 = d3.sgl->size();

  cout << "num of size-3 patterns: " << npat3 << endl;
    


  vector<SGList> sgls = { d3, d3 };

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  cout << "build table done" << endl;
    
    for(int i=0; i < H[0].size(); i++)
    {
        double sum = 0, max = H[0][i]->count[0];
        int Len = H[0][i]->count.size();
        for(int j=0; j < Len; j++)
        {
            sum += H[0][i]->count[j];
            if (max < H[0][i]->count[j])
                max = H[0][i]->count[j];
        }
        sum /= Len;
        cout << "H[0][" << i << "]: number=" << Len << ", avg=" << sum << ", max=" << max << endl;
    }

//  Sampler *sm2;
//  if (st2 > 0)
//    sm2 = new ProportionalSampler({ st2, st2 });
//  else sm2 = &default_sampler;
//
//  auto query = MyQuery();
//
//  util::Timer t;
//  t.start();
//  auto [d_res, ess] = join<false, true, false, false, 2, 4, 4>(g, H, sgls, false, *sm2, -1, false, st2 > 0, false, join_dummy1, query);
//  t.stop();
//
//  cout << "Time: " << t.get() << " sec, ";
//  if (d_res.sgl) {
//    if (st2 == 0) {
//      cout << "Num patterns: " << d_res.sgl->keys.size() << endl;
//
//      vector<double> counts(d_res.sgl->count.begin(), d_res.sgl->count.end());
//      sort(counts.begin(), counts.end(), std::greater<double>());
//
//      for (int i = 0; i < (50 > counts.size() ? counts.size() : 50); i++) {
//        cout << counts[i] << endl;
//      }
//    }
//    else {
//      cout << "Num patterns: " << ess.size() << endl;
//
//      vector<double> counts;
//
//      for (auto& [k, v] : ess) {
//        counts.push_back(v);
//      }
//
//      sort(counts.begin(), counts.end(), std::greater<double>());
//
//      for (int i = 0; i < (50 > counts.size() ? counts.size() : 50); i++) {
//        cout << counts[i]*256*256 << endl;
//        cout << counts[i]*32*32*32*32 << endl;
//      }
//    }
//
//  }
//  else {
//    cout << "Num patterns: 0" << endl;
//  }

  return 0;
}