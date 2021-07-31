#include "join.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>  // generate
#include <algorithm>
#include <functional>  // bind
#include <future>
#include <iterator>  // begin, end, and ostream_iterator
#include <queue>
#include <set>
#include <string>
#include <type_traits>
#include <vector>  // vector

#include "db/db.h"
#include "util.h"

using namespace std;

namespace euler::pattern_mining {



#ifdef PROF
  atomic<size_t> memory_load_count = 0;
  atomic<size_t> memory_load_size = 0;
  atomic<size_t> iso_test_count = 0;
  atomic<size_t> auto_test_count = 0;
#endif

  


  std::pair<std::unordered_set<unsigned long>, std::unordered_set<unsigned long>> join_dummy1 = {};

  Sampler default_sampler = Sampler();
  Query default_query = Query();


  bool bsearch(const size_t* x, size_t s, int y) {
    long long low = 0;
    long long high = s - 1;
    while (low <= high) {
      size_t mid = (low + high) / 2;
      if (x[mid] == y) {
        return true;
      }
      else if (x[mid] < y) {
        low = mid + 1;
      }
      else {
        high = mid - 1;
      }
    }
    return false;
  }


  

vector<vector<shared_ptr<db::MyKV<int>>>> build_tables(const vector<SGList>& sgls) {
    vector<vector<shared_ptr<db::MyKV<int>>>> H;
    std::vector<std::shared_ptr<std::vector<std::map<int, std::map<int, double>>>>> subgraph_hist;

    int res_size = 2;
    map<SGList, vector<shared_ptr<db::MyKV<int>>>> processed;
    for (auto& sgl : sgls) {
      int n = sgl.sgl->ncols - 1;
      res_size += n - 1;
      auto it = processed.find(sgl);
      if (it != processed.end()) {
        H.push_back(it->second);
      }
      else {
        vector<shared_ptr<db::MyKV<int>>> ht;
        for (int i = 0; i < n; i++) {
          ht.push_back(make_shared<db::MyKV<int>>(sgl.sgl->ncols));
          for (auto kv1 : sgl.sgl->keys) {
            auto buf = sgl.sgl->getbuf(kv1.second, sgl.sgl->ncols);
            auto it_buf = buf.begin();
            while (true) {
              for (size_t z = 0; z < it_buf.buffer_size / sgl.sgl->ncols; z++) {
                const int* it_d = it_buf.buffer + z * sgl.sgl->ncols;
                int kk = it_d[i + 1];
                ht.back()->merge(kk, it_d, sgl.sgl->ncols * sizeof(int));
              }
              if (!it_buf.has_next) break;
              it_buf.next();
            }
          }
        }
        H.push_back(ht);
        processed[sgl] = ht;
      }
    }

    return H;
  }

 

}  // namespace euler::pattern_mining
