#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;


// citeseer 4,5,6,7  // sampling ratio: 2,3,4
// mico 4
int main(int argc, char* argv[]) {
  // system("rm test_temp/*");

  graph::Graph_CSR_CPU g;


  g.read_graph(argv[1]);

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));

   auto pat5 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(5));


  double st2 = atof(argv[2]);

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true);

  cout << "start join for pat3: " << endl;
  vector<SGList> sgls2 = { d2, d2 };

  util::Timer match_time;
  match_time.start();
  auto H2 = build_tables(sgls2);

  auto [d3, ess3] = join<true, true, false, false, 2, 3, 3>(g, H2, sgls2, true, default_sampler, -1, true);

  match_time.stop();

  cout << "join for pat3 time: " << match_time.get() << " sec" << endl;

  size_t npat3 = d3.sgl->size();

  cout << "num of size-3 patterns: " << npat3 << endl;


  vector<SGList> sgls = { d3, d3 };

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  cout << "build table done" << endl;

  Sampler* sm2;
  if (st2 > 0)
    sm2 = new ProportionalSampler2({ st2, st2 });
  else sm2 = &default_sampler;

  util::Timer t;
  t.start();
  auto [d_res, ess] = join<true, true, false, false, 2, 4, 4>(g, H, sgls, false, *sm2, -1, false, st2 > 0);
  t.stop();

  double timelimit = t.get();

  cout << "Time: " << timelimit << " sec" << endl;



  SGList asap_res;
  double asap_rounds = 0;

  {
    cout << "start asap sampling: " << endl;

    double tot_time = 0.0;

    while (tot_time < timelimit) {

      util::Timer t;
      t.start();
      auto d2 = asap::match(g, pat5, false, false, true, -1, false, false, 1);


      if (asap_res.sgl == nullptr) asap_res = d2;
      else asap_res.combine_count(d2);

      asap_rounds += 1;

      t.stop();
      tot_time += t.get();
    }
  }


  if (d_res.sgl) {
    if (st2 == 0) {
      cout << "Num patterns: " << d_res.sgl->keys.size() << endl;

      vector<pair<string, size_t>> result;
      for (auto& kv : d_res.sgl->keys) {
        result.push_back(kv);
      }

      sort(result.begin(), result.end(), [&](const pair<string, size_t>& a, const pair<string, size_t>& b) {
        return d_res.sgl->count[a.second] > d_res.sgl->count[b.second];
        });

      int found_by_asap = 0;
      for (int i = 0; i < result.size(); i++) {
        if (asap_res.sgl->keys.find(result[i].first) != asap_res.sgl->keys.end()) {
          if (i<50) found_by_asap++;
          cout << d_res.sgl->count[result[i].second] << "\t" << asap_res.sgl->count[asap_res.sgl->keys[result[i].first]] / asap_rounds << endl;
        }
      }
      cout << "num of first 50 patterns found by asap: " << found_by_asap << endl;
    }
    else {
      cout << "Num patterns: " << ess.size() << endl;

      vector<pair<string, size_t>> result;
      for (auto& kv : ess) {
        result.push_back(kv);
      }


            cout << "num patterns returned by asap: " << size_t(asap_res.sgl->size() * 18.7) << endl;


      sort(result.begin(), result.end(), [&](const pair<string, size_t>& a, const pair<string, size_t>& b) {
        return a.second > b.second;
        });


      size_t asap_50 = 0;
      for (int i = 0; i < (50 > result.size() ? result.size() : 50); i++) {
        //cout << result[i].second << "\t";
        if (asap_res.sgl->keys.find(result[i].first) != asap_res.sgl->keys.end()) {
          //cout << asap_res.sgl->count[asap_res.sgl->keys[result[i].first]] / asap_rounds << endl;
          asap_50+=1;
        }
        else {
          //cout << 0 << endl;
        }
      }
      cout<< "num of top 50 patterns returned by asap: "<<asap_50<<endl;
    }
  }
  else {
    cout << "Num patterns: 0" << endl;
  }

  return 0;
}
