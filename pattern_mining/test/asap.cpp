#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;

typedef vector<pair<int, int>> pat_t;



std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::map<int, double>>>>>  get_edge_hist(const std::vector<SGList>& sgls) {
  auto subgraph_hist = std::make_shared<std::vector<std::shared_ptr<std::vector<std::map<int, double>>>>>();

  auto degrees = std::make_shared<std::vector<std::map<int, double>>>();
  for (int i = 0; i < 2; i++) {
    degrees->push_back(std::map<int, double>());

    auto& his = degrees->back();

    for (auto kv1 : sgls[0].sgl->keys) {
      auto buf = sgls[0].sgl->getbuf(kv1.second, sgls[0].sgl->ncols);
      auto it_buf = buf.begin();
      while (true) {
        for (size_t z = 0; z < it_buf.buffer_size / sgls[0].sgl->ncols; z++) {
          const int* it_d = it_buf.buffer + z * sgls[0].sgl->ncols;
          int kk = it_d[i + 1];
          if (his.find(kk) == his.end()) his[kk] = 0;
          his[kk]++;
        }
        if (!it_buf.has_next) break;
        it_buf.next();
      }
    }
  }
  for (auto &sg: sgls) {
    subgraph_hist->push_back(degrees);
  }
  return subgraph_hist;
}




class ASAPSampler : public Sampler {
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::map<int, double>>>>> edge_hist;
  std::vector<int> num_samples;
public:
  ASAPSampler(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::map<int, double>>>>> hist) : edge_hist(hist) {
    num_samples.resize(hist->size(), 0);
  }
  inline SampleRes smp_prob(util::span<const int> s, util::span<const int> t, int i, int j, int step) {
    if (num_samples[step] == 0) {
      num_samples[step]++;
      return {BREAK, 1 / edge_hist->at(step)->at(j).at(t[j + 1])};
    }
    else {
      return 0;
    }
  }
};


int main(int argc, char* argv[]) {
  // system("rm test_temp/*");

  graph::Graph_CSR_CPU g;


  g.read_graph(argv[1]);

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true);


  vector<SGList> sgls = { d2, d2, d2 };

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  auto edge_hist = get_edge_hist(sgls); 
  cout << "build table done" << endl;

  Sampler *sm = new ASAPSampler(edge_hist);

  double tot_time = 0;
  SGList tot_res;

  for (int i = 0; i < 10; i++) {

    util::Timer t;
    t.start();
    auto [d_res, ess] = join<true, true, false, false, 3, 3, 3, 3>(g, H, sgls, false, *sm, -1, false, true);
    t.stop();

    if (tot_res.sgl == nullptr) tot_res = d_res;
    else
      tot_res.combine(d_res, true, false);

    tot_time += t.get();

    cout << "Time: " << tot_time << " sec, ";
    if (tot_res.sgl) {
      cout << "Num patterns: " << tot_res.sgl->keys.size() << endl;
    }
    else {
      cout << "Num patterns: 0" << endl;
    }

  }


  return 0;
}
