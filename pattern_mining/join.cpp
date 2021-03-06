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
#include <random>  // mt19937 and uniform_int_distribution
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

static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<> dis(0, 1);

double random_number() {
  return dis(gen);
}

bool bsearch(const int *x, size_t s, int y) {
  int low = 0;
  int high = s - 1;
  while (low <= high) {
    int mid = (low + high) / 2;
    if (x[mid] == y) {
      return true;
    } else if (x[mid] < y) {
      low = mid + 1;
    } else {
      high = mid - 1;
    }
  }
  return false;
}

/*
vector<pair<connection_t, shared_ptr<Pattern>>> get_connectivity_edge_induced(
    const vector<pair<const int *, size_t>> &nbv, const int *d1, const int *d2,
    int ncols1, int ncols2, int c1, int c2,
    shared_ptr<Pattern> pat1, shared_ptr<Pattern> pat2, vector<int> &value,
    bool has_labels) {
  for (int xi = 1; xi < ncols1; xi++) {
    int n1 = d1[xi];
    for (int yi = 1; yi < ncols2; yi++) {
      int n2 = d2[yi];
      if ((xi != c1 || yi != c2) && n1 == n2) {
        return {make_pair(
            1 << 31, nullptr)};  // support join of maximum size-5 and size-6
      }
    }
  }

  value.reserve(ncols1 + ncols2 - 2);
  value.push_back(0);
  vector<int> idx2;
  idx2.reserve(ncols2 - 1);

  auto ptt = std::make_shared<Pattern>();
  ptt->enable_label();
  int jp;
  for (int xi = 1; xi < ncols1; xi++) {
    value.push_back(d1[xi]);
    int lab = has_labels ? pat1->vertex_label[xi - 1] : 0;
    ptt->add_node(lab);
  }  // left < ncols1
  int txi = 0;
  for (int xi = 1; xi < ncols2; xi++) {
    if (xi == c2) {
      idx2.push_back(c1 - 1);
    } else {
      idx2.push_back((txi++) + ncols1 - 1);
      value.push_back(d2[xi]);
      int lab = has_labels ? pat2->vertex_label[xi - 1] : 0;
      ptt->add_node(lab);
    }
  }
  jp = c1;

  for (int ti = 0; ti < pat1->adj_list.size(); ti++) {
    for (int tmp : pat1->adj_list[ti]) {
      if (tmp > ti) ptt->add_edge(ti, tmp);
    }
  }

  for (int ti = 0; ti < pat2->adj_list.size(); ti++) {
    for (int tmp : pat2->adj_list[ti]) {
      if (tmp > ti) ptt->add_edge(idx2[ti], idx2[tmp]);
    }
  }

  connection_t key_full = 0;

  for (int yi = 0; yi < nbv.size(); yi++) {
    for (int xi = 1; xi < ncols1; xi++) {
      if (xi == jp) continue;
      if (bsearch(nbv[yi].first, nbv[yi].second, value[xi])) {
        ptt->add_edge(xi - 1, yi + ncols1 - 1);
        key_full |= (1 << ((xi - 1) * (ncols2 - 1) + (yi)));
      }
    }
  }

  vector<pair<connection_t, shared_ptr<Pattern>>> res;

  unsigned int nbits = 0;
  vector<int> bits_idx;

  for (int i = 0; i < 32; i++) {
    if (key_full & (1 << i)) {
      nbits++;
      bits_idx.push_back(i);
    }
  }

  for (unsigned int i = 0; i < (1 << nbits); i++) {
    // for each (key, ptt) in all possible connections
    connection_t key = 0;
    shared_ptr<Pattern> pttt = make_shared<Pattern>(*ptt);

    for (int j = 0; j < 32; j++) {
      if (i & (1 << j)) {
        int c = bits_idx[j];
        int x = c / (ncols2 - 1);
        int y = c % (ncols2 - 1);
        pttt->add_edge(x, y + ncols1 - 1);
        key |= (1 << c);
      }
    }

    // insertion sort, sort the vertices by ID
    vector<int> start_vec;
    start_vec.reserve(value.size());
    for (int xi = 1; xi < value.size(); xi++) {
      int pos = 0;
      while (pos < start_vec.size() && value[start_vec[pos]] < value[xi]) pos++;
      start_vec.insert(start_vec.begin() + pos, xi);
    }
    

    vector<int> visited(value.size(), 0);

#ifdef PROF
    auto_test_count++;
#endif
    for (int start : start_vec) {
      for (int z = 0; z < visited.size(); z++) visited[z] = 0;
      int ret = list_left(pttt, start, ncols1 - 2, visited, value, jp);

      if (ret == -1) {
        return {make_pair(1 << 31, nullptr)};
      }
      if (ret == 1) break;
    }

    unsigned int left_mask = 0;
    for (int xi = 1; xi < value.size(); xi++)
      if (visited[xi]) left_mask |= (1 << xi);
    unsigned int left_mask1 = 0;
    for (int xi = 1; xi < ncols1; xi++) {
      left_mask1 |= (1 << xi);
    }

    if (left_mask1 != left_mask) {
      return {make_pair(1 << 31, nullptr)};
    }
    

    res.emplace_back(key, pttt);
  }
  return res;
}
*/

/*
void for_loop2_edge_induced(const vector<SGList> &sgls, vector<int> &s,
                            shared_ptr<Pattern> pat,
                            vector<vector<shared_ptr<db::MyKV<int>>>> &H,
                            int level, vector<int> &iterates, SGList &res,
                            vector<map<int, string>> &qp2cp,
                            vector<vector<int>> &qp_count,
                            vector<vector<map<vector<int>, int>>> &qp_idx,
                            const graph::Graph &g, bool has_labels,
                            SamplingMethod sm, vector<double> sampling_param,
                            cpp_int &num_samples, cpp_int &total_population) {
  int tid = omp_get_thread_num();

  if (H.size() == level) {
    auto it = qp2cp[tid].find(s[0]);
    if (it != qp2cp[tid].end()) {
      res.sgl->merge(it->second, s.data(), s.size() * sizeof(int), false);
    } else {
#ifdef PROF
      iso_test_count++;
#endif
      string coding = pat->dfs_coding();
      qp2cp[tid][s[0]] = coding;
      res.sgl->merge(coding, s.data(), s.size() * sizeof(int), false);
    }

    total_population = 1;
    num_samples++;

  } else {
    auto &pats1 = sgls[level].patterns;

    for (int i = 1; i < s.size(); i++) {
      if (H[level][iterates[level]]->keys.find(s[i]) ==
          H[level][iterates[level]]->keys.end())
        continue;
#ifdef PROF
      memory_load_count++;
#endif

      //vector<int> buf;

      int ncols = H[level][iterates[level]]->ncols;

      auto [buf, size] = H[level][iterates[level]]->getbuf(get<0>(H[level][iterates[level]]->keys[s[i]]));

#ifdef PROF
      memory_load_size += size * sizeof(int);
#endif

      size_t len = size / ncols;
      //size_t lena;

      int j = iterates[level];

      cpp_int tp_estimate = 0;

      //  vector<size_t> indices;
      

      int type1 = s[0];
      vector<pair<const int *, size_t>> nbv;
      for (int li = 1; li < s.size(); li++) {
        if (li == i) continue;
        nbv.push_back(g.get_neighbors(s[li]));
      }

      for (size_t z = 0; z < len; z++) {
        // cout << z << " " << z1 << endl;

        int type2 = buf[z * ncols];

        const int *it_buf = buf + z * ncols;

        vector<int> value;

        auto key_pat_vec = get_connectivity_edge_induced(
            nbv, it_buf, s.data(), ncols, s.size(), j + 1, i, pats1.at(type2), pat,
            value, has_labels);

        for (int kpi = 0; kpi < key_pat_vec.size(); kpi++) {
          auto key2 = key_pat_vec[kpi].first;
          auto &ptt = key_pat_vec[kpi].second;

          if (key2 == (1 << 31)) {
            continue;
          }

          vector<int> key;

          unsigned int ta = ((i - 1) * (ncols - 1) + j);
          key.push_back(type1);
          key.push_back(type2);
          key.push_back(ta);
          key.push_back(key2);

          if (has_labels) {
            for (int xi = 1; xi < s.size(); xi++) {
              key.push_back(g.get_vertex_label(s[xi]));
            }
            for (int xi = 1; xi < ncols; xi++) {
              if (xi == j + 1) continue;
              key.push_back(g.get_vertex_label(it_buf[xi]));
            }
          }

          {
            auto itt = qp_idx[tid][level].find(key);
            if (itt == qp_idx[tid][level].end()) {
              qp_idx[tid][level][key] = qp_count[tid][level];
              value[0] = qp_count[tid][level];
              qp_count[tid][level]++;
            } else {
              value[0] = itt->second;
            }
          }

          cpp_int tp = 0;
          for_loop2_edge_induced(sgls, value, ptt, H, level + 1, iterates, res, qp2cp,
                                 qp_count, qp_idx, g, has_labels, sm, sampling_param,
                                 num_samples, tp);

          tp_estimate = tp + tp_estimate;
        }
      }
    }
  }
}
*/

/*
void for_loop1_edge_induced(const vector<SGList> &sgls,
                            vector<vector<shared_ptr<db::MyKV<int>>>> &H,
                            vector<int> &iterates, int level, SGList &res,
                            vector<map<int, string>> &qp2cp,
                            vector<vector<int>> &qp_count,
                            vector<vector<map<vector<int>, int>>> &qp_idx,
                            const graph::Graph &g, bool has_labels,
                            SamplingMethod sm, vector<double> sampling_param, bool store) {
  if (level < H.size()) {
    for (int i = 0; i < H[level].size(); i++) {
      iterates.push_back(i);
      for_loop1_edge_induced(sgls, H, iterates, level + 1, res, qp2cp, qp_count, qp_idx, g,
                             has_labels, sm, sampling_param, store);
      iterates.pop_back();
    }
  } else {
    euler::util::print_vec(iterates);

    int i = iterates[0];
    int j = iterates[1];

    auto &tab1 = H[0][i];
    auto &tab2 = H[1][j];

    auto &pats1 = sgls[0].patterns;
    auto &pats2 = sgls[1].patterns;

    int ncols1 = tab1->ncols;
    int ncols2 = tab2->ncols;

    // cout << "tab size: " << tab1->keys.size() << "\t" << tab2->keys.size() <<
    // endl;
    atomic<cpp_int> num_samples(0);
    cpp_int total_population = 0;

    for (auto bi = tab1->keys.begin(); bi != tab1->keys.end(); bi++) {
#ifdef PROF
      memory_load_count++;
#endif
      //vector<int> reordered_d1;
      //tab1->getbuf(reordered_d1, get<0>(bi->second));
      auto [reordered_d1, size1] = tab1->getbuf(get<0>(bi->second));
      auto it_tab2 = tab2->keys.find(bi->first);

#ifdef PROF
      memory_load_size += size1 * sizeof(int);
#endif

      if (it_tab2 != tab2->keys.end()) {
#ifdef PROF
        memory_load_count++;
#endif

        //vector<int> reordered_d2;
        //tab2->getbuf(reordered_d2, get<0>(it_tab2->second));
        auto [reordered_d2, size2] = tab2->getbuf(get<0>(it_tab2->second));

#ifdef PROF
        memory_load_size += size2 * sizeof(int);
#endif

        atomic<cpp_int> tp_estimate(0);

        // now we have reordered_d1 and reordered_d2 ready
        // the first element in each row is the pattern_idx
        size_t len1 = size1 / ncols1;
        size_t len2 = size2 / ncols2;

        

#pragma omp parallel for num_threads(_Nthreads)
        for (size_t z = 0; z < len1; z++) {
          int tid = omp_get_thread_num();
          int type1 = reordered_d1[z * ncols1];
          const int *it_d1 = reordered_d1 + z * ncols1;

          vector<pair<const int *, size_t>> nbv;
          for (int li = 1; li < ncols1; li++) {
            if (li == i + 1) continue;
            nbv.push_back(g.get_neighbors(it_d1[li]));
          }

          for (size_t z1 = 0; z1 < len2; z1++) {
            int tid = omp_get_thread_num();

            int type2 = reordered_d2[z1 * ncols2];

            vector<int> value;

            const int *it_d2 = reordered_d2 + z1 * ncols2;
            const int *it_d1 = reordered_d1 + z * ncols1;

            auto key_pat_vec = get_connectivity_edge_induced(
                nbv, it_d2, it_d1, ncols2, ncols1, j + 1,
                i + 1, pats2.at(type2), pats1.at(type1), value, has_labels);

            for (auto [key2, ptt] : key_pat_vec) {
              if (key2 == (1 << 31)) {
                continue;
              }

              vector<int> key;

              unsigned int ta = (i * (ncols2 - 1) + j);
              key.push_back(type1);
              key.push_back(type2);
              key.push_back(ta);
              key.push_back(key2);

              if (has_labels) {
                for (int xi = 1; xi < ncols1; xi++) {
                  key.push_back(
                      g.get_vertex_label(it_d1[xi]));
                }
                for (int xi = 1; xi < ncols2; xi++) {
                  if (xi == j + 1) continue;
                  key.push_back(
                      g.get_vertex_label(it_d2[xi]));
                }
              }

              // now we get the key for quick pattern: if no label is
              // considered, the key is simply the topological pattern,
              // otherwise, the topological pattern appended with the vertex
              // labels

              // quick pattern: qp_tmp only stores one value

              auto itt = qp_idx[tid][0].find(key);
              if (itt == qp_idx[tid][0].end()) {
                // the first time we found a quick pattern, compute its
                // topological canonical form

                qp_idx[tid][0][key] = qp_count[tid][0];
                value[0] = qp_count[tid][0];
                qp_count[tid][0]++;

              } else {
                value[0] = itt->second;
              }
              cpp_int ns = 0, tp = 0;
              for_loop2_edge_induced(sgls, value, ptt, H, 2, iterates, res, qp2cp, qp_count,
                                     qp_idx, g, has_labels, sm, sampling_param, ns, tp);

              tp_estimate = tp_estimate + tp;
              num_samples = ns + num_samples;
            }
          }
        }
      }
    }
    res.num_samples += num_samples;
    res.total_space += total_population;
  }
}
*/

vector<vector<shared_ptr<db::MyKV<int>>>> build_tables(const vector<SGList> &sgls) {
  vector<vector<shared_ptr<db::MyKV<int>>>> H;
  int res_size = 2;
  map<SGList, vector<shared_ptr<db::MyKV<int>>>> processed;
  for (auto &sgl : sgls) {
    int n = sgl.sgl->ncols - 1;
    res_size += n - 1;
    auto it = processed.find(sgl);
    if (it != processed.end()) {
      H.push_back(it->second);
    } else {
      vector<shared_ptr<db::MyKV<int>>> ht;
      for (int i = 0; i < n; i++) {
        ht.push_back(make_shared<db::MyKV<int>>(sgl.sgl->ncols));
        for (auto kv1 : sgl.sgl->keys) {
          auto buf = sgl.sgl->getbuf(kv1.second, sgl.sgl->ncols);
          auto it_buf = buf.begin();
          while (true) {
            for (size_t z = 0; z < it_buf.buffer_size / sgl.sgl->ncols; z++) {
              const int *it_d = it_buf.buffer + z * sgl.sgl->ncols;
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
/*
vector<vector<shared_ptr<db::MyKV<int>>>> build_tables(const vector<SGList> &sgls) {
  vector<vector<shared_ptr<db::MyKV<int>>>> H;
  int res_size = 2;
  for (auto &sgl : sgls) {
    int n = sgl.sgl->ncols - 1;
    res_size += n - 1;
    vector<shared_ptr<db::MyKV<int>>> ht;
    for (int i = 0; i < n; i++) {
      ht.push_back(make_shared<db::MyKV<int>>(sgl.sgl->ncols));
      for (auto kv1 : sgl.sgl->keys) {
        auto buf = sgl.sgl->getbuf(get<0>(kv1.second), sgl.sgl->ncols);
        auto it_buf = buf.begin();
        while (true) {
#pragma omp parallel for num_threads(_Nthreads)
          for (size_t z = 0; z < it_buf.buffer_size / sgl.sgl->ncols; z++) {
            const int *it_d = it_buf.buffer + z * sgl.sgl->ncols;
            int kk = it_d[i + 1];
            ht.back()->merge(kk, it_d, sgl.sgl->ncols * sizeof(int));
          }
          if (!it_buf.has_next) break;
          it_buf.next();
        }
      }
    }
    H.push_back(ht);
  }
  return H;
}
*/

}  // namespace euler::pattern_mining
