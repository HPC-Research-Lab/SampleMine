#pragma once

#include <limits.h>
#include <omp.h>

#include <array>
#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <type_traits>


#include "graph/graph.h"
#include "pattern.h"
#include "types.h"
#include "sampler.h"

namespace euler::pattern_mining {

  const static int _Nthreads = atoi(getenv("OMP_NUM_THREADS"));
  const static size_t MAX_NUM_QPATTERN_PER_THREAD = 1024l * 1024l;

  extern std::pair<std::unordered_set<unsigned long>, std::unordered_set<unsigned long>> join_dummy1;

  std::vector<std::vector<std::shared_ptr<db::MyKV<int>>>> build_tables(const std::vector<SGList>& sgls);



  extern Sampler default_sampler;
  extern Query default_query;




  bool bsearch(const size_t* x, size_t s, int y);


  template <size_t ncols>
  bool is_connected(std::shared_ptr<Pattern> pat, std::array<int, ncols> visited) {
    int source;

    /* if (omp_get_thread_num() == 0) {
       t_is_connected.start();
     }*/
    for (int i = 1; i < visited.size(); i++)
      if (!visited[i]) {
        source = i;
        break;
      }
    std::queue<int> q;
    q.push(source);
    while (!q.empty()) {
      int cur = q.front();
      q.pop();
      visited[cur] = 1;
      for (int nb : pat->adj_list[cur - 1]) {
        if (!visited[nb + 1]) q.push(nb + 1);
      }
    }
    /*if (omp_get_thread_num() == 0) {
      t_is_connected.stop();
    }*/
    for (int i = 1; i < visited.size(); i++)
      if (!visited[i]) return false;
    return true;
  }

  template <size_t ncols>
  int list_left(std::shared_ptr<Pattern> pat, int pos, int n, std::array<int, ncols>& visited,
    std::array<int, ncols>& value, int jp) {
    visited[pos] = 1;
    if (n == 0) {
      std::array<int, visited.size()> visited_tmp = visited;
      bool flag = false;
      for (int xi = 1; xi < visited_tmp.size(); xi++) {
        if (visited_tmp[xi] == 1) {
          visited_tmp[xi] = 0;
          if (is_connected(pat, visited_tmp)) {
            flag = true;
            if (value[xi] < value[jp]) return -1;
          }
          visited_tmp[xi] = 1;
        }
      }
      if (flag)
        return 1;
      else {
        visited[pos] = 0;
        return 0;
      }
    }
    else {
      std::array<int, visited.size()> nbs{};

      unsigned nelem = 0;
      for (int xi = 1; xi < visited.size(); xi++) {
        if (visited[xi]) {
          for (int nb : pat->adj_list[xi - 1]) {
            if (!visited[nb + 1]) {
              int yi;
              for (yi = 0; yi < nelem; yi++) {
                if (value[nbs[yi]] >= value[nb + 1]) break;
              }
              // if it is not a duplicate
              if (nbs[yi] != nb + 1) {
                unsigned yyi = nelem;
                while (yyi > yi) {
                  nbs[yyi] = nbs[yyi - 1];
                  yyi--;
                }
                nbs[yyi] = nb + 1;
                nelem++;
              }
            }
          }
        }
      }

      for (int nb : nbs) {
        if (nb > 0) {
          int ret = list_left(pat, nb, n - 1, visited, value, jp);
          if (ret == 1)
            return 1;
          else if (ret == -1)
            return -1;
          visited[nb] = 0;
        }
      }
      return 0;
    }
  }

  typedef unsigned int connection_t;

  template <bool has_labels, bool edge_induced, size_t ncols1, size_t ncols2, size_t ncols>
  std::vector<std::pair<connection_t, std::shared_ptr<Pattern>>> get_connectivity(
    const std::array<std::pair<const size_t*, size_t>, ncols2 - 2>& nbv, const int* d1, const int* d2,
    int c1, int c2,
    std::shared_ptr<Pattern> pat1, std::shared_ptr<Pattern> pat2, std::array<int, ncols>& value,
    const std::pair<std::unordered_set<unsigned long>, std::unordered_set<unsigned long>>& sgl3) {


    for (int xi = 1; xi < ncols1; xi++) {
      int n1 = d1[xi];
      for (int yi = 1; yi < ncols2; yi++) {
        int n2 = d2[yi];
        if ((xi != c1 || yi != c2) && n1 == n2) {
          return { std::make_pair(1 << 31,
                                 nullptr) };  // support join of maximum size-5 and size-6
        }
      }
    }

    value[0] = 0;
    std::array<int, ncols2 - 1> idx2;

    auto ptt = std::make_shared<Pattern>();
    ptt->enable_label();
    int jp;
    int t = 1;
    for (int xi = 1; xi < ncols1; xi++) {
      value[t++] = d1[xi];
      int lab = has_labels ? pat1->vertex_label[xi - 1] : 0;
      ptt->add_node(lab);
    }  // left < ncols1
    int txi = 0;
    int idx2_i = 0;
    for (int xi = 1; xi < ncols2; xi++) {
      if (xi == c2) {
        idx2[idx2_i++] = (c1 - 1);
      }
      else {
        idx2[idx2_i++] = ((txi++) + ncols1 - 1);
        value[t++] = d2[xi];
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
    unsigned int nbits = 0;
    std::vector<int> bits_idx;

    for (int yi = 0; yi < nbv.size(); yi++) {
      for (int xi = 1; xi < ncols1; xi++) {
        if (xi == jp) continue;
        if (bsearch(nbv[yi].first, nbv[yi].second, value[xi])) {
          if constexpr (!edge_induced) ptt->add_edge(xi - 1, yi + ncols1 - 1);
          int idxx = (xi - 1) * (ncols2 - 1) + (yi);
          key_full |= (1 << idxx);
          nbits++;
          bits_idx.push_back(idxx);
        }
      }
    }

    std::vector<std::pair<connection_t, std::shared_ptr<Pattern>>> res;

    unsigned int num_combinations;
    if constexpr (edge_induced)
      num_combinations = (1 << nbits);
    else
      num_combinations = 1;

    /*std::cout << "((((((((" << std::endl;
    pat1->print();
    pat2->print();
    ptt->print();
    std::cout << "))))))))" << std::endl;*/

    for (unsigned int i = 0; i < num_combinations; i++) {
      connection_t key = 0;
      std::shared_ptr<Pattern> pttt = std::make_shared<Pattern>(*ptt);

      if constexpr (edge_induced) {
        for (int j = 0; j < 32; j++) {
          if (i & (1 << j)) {
            int c = bits_idx[j];
            int x = c / (ncols2 - 1);
            int y = c % (ncols2 - 1);
            pttt->add_edge(x, y + ncols1 - 1);
            key |= (1 << c);
          }
        }
      }
      else {
        key = key_full;
      }

      if constexpr (has_labels) {
        if (sgl3.first.size() > 0 || sgl3.second.size() > 0) {
          bool skip = false;
          for (int n1 : pat1->adj_list[c1 - 1]) {
            for (int n2 : pat2->adj_list[c2 - 1]) {
              int idx2 = n2;
              if (n2 > c2 - 1) idx2--;
              unsigned long lab = 0;
              int labs[3];
              labs[0] = pat1->vertex_label[c1 - 1];
              labs[1] = pat1->vertex_label[n1];
              labs[2] = pat2->vertex_label[n2];
              lab |= (1ULL << labs[0]);
              lab |= (1ULL << labs[1]);
              lab |= (1ULL << labs[2]);
              if (key & (1 << (n1 * (ncols2 - 1) + n2))) {
                if ((labs[0] == labs[1]) && (labs[0] < labs[2])) lab |= (1ULL << 63);
                if ((labs[0] == labs[2]) && (labs[0] < labs[1])) lab |= (1ULL << 63);
                if ((labs[1] == labs[2]) && (labs[1] < labs[0])) lab |= (1ULL << 63);
                if (sgl3.first.find(lab) == sgl3.first.end()) {
                  skip = true;
                  break;
                }
              }
              else {
                if ((labs[1] > labs[0]) && (labs[2] > labs[0]))
                  lab |= (1ULL << 62);
                else if ((labs[1] > labs[0]) || (labs[2] > labs[0]))
                  lab |= (1ULL << 63);
                else if (((labs[1] == labs[0]) && (labs[2] != labs[0])) || ((labs[1] != labs[0]) && (labs[2] == labs[0])))
                  lab |= (1ULL << 61);
                if (sgl3.second.find(lab) == sgl3.second.end()) {
                  skip = true;
                  break;
                }
              }
            }
            if (skip) break;
          }
          if (skip) continue;
        }
      }

      std::array<int, value.size() - 1> start_vec;
      unsigned nelem = 0;
      for (int xi = 1; xi < value.size(); xi++) {
        int j = nelem;
        while (j > 0 && value[start_vec[j - 1]] > value[xi]) {
          start_vec[j] = start_vec[j - 1];
          j--;
        }
        start_vec[j] = xi;
        nelem++;
      }

      std::array<int, ncols> visited{};

      //  if (omp_get_thread_num() == 0) t_list_left.start();
#ifdef PROF
      auto_test_count++;
#endif
      bool valid = true;
      for (int start : start_vec) {
        for (int z = 0; z < visited.size(); z++) visited[z] = 0;
        int ret = list_left<visited.size()>(pttt, start, ncols1 - 2, visited, value, jp);

        if (ret == -1) {
          //    if (omp_get_thread_num() == 0) t_list_left.stop();
          valid = false;
          break;
        }
        if (ret == 1) break;
      }
      //  if (omp_get_thread_num() == 0) t_list_left.stop();

      if (!valid) continue;

      unsigned int left_mask = 0;
      for (int xi = 1; xi < value.size(); xi++)
        if (visited[xi]) left_mask |= (1 << xi);
      unsigned int left_mask1 = 0;
      for (int xi = 1; xi < ncols1; xi++) {
        left_mask1 |= (1 << xi);
      }

      if (left_mask1 != left_mask) {
        continue;
      }

      res.emplace_back(key, pttt);
    }

    // assert(key != (1 << 31));
    return res;
  }

  template <bool pat_agg, bool mni, size_t ncols_left>
  void for_loop2_end(const graph::Graph& g, std::array<int, ncols_left>& s,
    std::shared_ptr<Pattern> pat, int level,
    std::vector<SGList>& res, std::vector<std::map<std::string, double>>& estimate_counts, std::vector<std::vector<double>>& sampling_probs, bool est,
    std::vector<std::map<int, typename std::conditional<mni, std::tuple<std::string, std::vector<std::vector<unsigned>>, std::vector<unsigned>>, std::string>::type>>& qp2cp, bool store, double mni_threshold, bool need_actual_pattern, std::map<std::shared_ptr<Pattern>, size_t, cmpByPattern>& actual_patterns, bool adaptive_sampling, Query& query, Sampler& sampler) {


    int tid = omp_get_thread_num();

    int q = query(g, { s.data(), s.size() }, pat, level - 1);

    if (q == -1) return;

    if constexpr (pat_agg) {

      auto it = qp2cp[tid].find(s[0]);
      if (it != qp2cp[tid].end()) {

        if (need_actual_pattern) {
          int pat_id = -1;
#pragma omp critical
          {
            if (actual_patterns.find(pat) == actual_patterns.end()) {
              pat_id = actual_patterns.size();
              actual_patterns[pat] = pat_id;
            }
            else {
              pat_id = actual_patterns[pat];
            }
          }
          s[0] = pat_id;
        }


        if constexpr (mni) {
          res[tid].sgl->merge(std::get<0>(it->second), s.data(), s.size() * sizeof(int), store, mni_threshold, std::get<1>(it->second), std::get<2>(it->second), adaptive_sampling);

          if (est) {
            double est_ct = ((ProportionalSampler*)&sampler)->C[tid].back();
            //std::cout << est_ct << std::endl;
            for (double& pr : sampling_probs[tid]) est_ct /= pr;
            auto itt = estimate_counts[tid].find(std::get<0>(it->second));
            if (itt != estimate_counts[tid].end()) {
              itt->second += est_ct;
            }
            else {
              estimate_counts[tid].insert({ std::get<0>(it->second), est_ct });
            }
          }

        }
        else {
          res[tid].sgl->merge(it->second, s.data(), s.size() * sizeof(int), store);
          if (est) {
            double est_ct = ((ProportionalSampler*)&sampler)->C[tid].back();
            //double est_ct = 1;
            for (double& pr : sampling_probs[tid]) {
              est_ct /= pr;
            }
            auto itt = estimate_counts[tid].find(it->second);
            if (itt != estimate_counts[tid].end()) {
              itt->second += est_ct;
            }
            else {
              estimate_counts[tid].insert({ it->second, est_ct });
            }
          }
        }
      }
      else {
#ifdef PROF
        iso_test_count++;
#endif
        // if (qp2cp[tid].size() > MAX_NUM_QPATTERN_PER_THREAD) {
         //  std::cerr << "qp2cp flushed" << std::endl;
         // qp2cp[tid].clear();
        // }

        int s0 = s[0];

        if (need_actual_pattern) {
          int pat_id = -1;
#pragma omp critical
          {
            if (actual_patterns.find(pat) == actual_patterns.end()) {
              pat_id = actual_patterns.size();
              actual_patterns[pat] = pat_id;
            }
            else {
              pat_id = actual_patterns[pat];
            }
          }
          s[0] = pat_id;
        }


        if constexpr (mni) {
          auto coding = pat->canonical_form();
          qp2cp[tid][s0] = coding;
          res[tid].sgl->merge(std::get<0>(coding), s.data(), s.size() * sizeof(int), store, mni_threshold, std::get<1>(coding), std::get<2>(coding), adaptive_sampling);

          if (est) {
            double est_ct = 1;
            for (double& pr : sampling_probs[tid]) est_ct /= pr;
            auto itt = estimate_counts[tid].find(std::get<0>(coding));
            if (itt != estimate_counts[tid].end()) {
              itt->second += est_ct;
            }
            else {
              estimate_counts[tid].insert({ std::get<0>(coding), est_ct });
            }
          }
        }
        else {
          std::string coding = pat->dfs_coding();
          qp2cp[tid][s0] = coding;
          res[tid].sgl->merge(coding, s.data(), s.size() * sizeof(int), store);


          if (est) {
            double est_ct = 1;
            for (double& pr : sampling_probs[tid]) est_ct /= pr;
            auto itt = estimate_counts[tid].find(coding);
            if (itt != estimate_counts[tid].end()) {
              itt->second += est_ct;
            }
            else {
              estimate_counts[tid].insert({ coding, est_ct });
            }
          }
        }
      }
    }
    else {
      std::string kk = std::to_string(q);
      res[tid].sgl->merge(kk, s.data(), s.size() * sizeof(int), store);

      if (est) {
        double est_ct = 1;
        for (double& pr : sampling_probs[tid]) est_ct /= pr;
        auto itt = estimate_counts[tid].find(kk);
        if (itt != estimate_counts[tid].end()) {
          itt->second += est_ct;
        }
        else {
          estimate_counts[tid].insert({ kk, est_ct });
        }
      }
    }
  }

  template <bool pat_agg, bool has_labels, bool edge_induced, bool mni, int K, typename key_type, size_t ncols_left, size_t ncols, size_t... ncols_right>
  void for_loop2(const std::vector<SGList>& sgls, std::array<int, ncols_left>& s,
    std::shared_ptr<Pattern> pat,
    const std::vector<std::vector<std::shared_ptr<db::MyKV<int>>>>& H, int level,
    std::vector<int>& iterates, std::vector<SGList>& res, std::vector<std::map<std::string, double>>& estimate_counts, std::vector<std::vector<double>>& sampling_probs,
    std::vector<std::map<int, typename std::conditional<mni, std::tuple<std::string, std::vector<std::vector<unsigned>>, std::vector<unsigned>>, std::string>::type>>& qp2cp, std::vector<std::vector<std::vector<std::array<int, 4>>>>& qp_count,
    std::vector<std::vector<std::map<key_type, int>>>& qp_idx,
    const graph::Graph& g, Sampler& sampler,
    bool store, const std::pair<std::unordered_set<unsigned long>, std::unordered_set<unsigned long>>& sgl3, double mni_threshold, bool need_actual_pattern, std::map<std::shared_ptr<Pattern>, size_t, cmpByPattern>& actual_patterns, bool est, bool adaptive_sampling, Query& query) {


    int qret = query(g, { s.data(), s.size() }, pat, level - 1);
    if (qret < 0) return;


    int tid = omp_get_thread_num();

    auto& pats1 = sgls[level].patterns;

    for (int i = 1; i < s.size(); i++) {
      if (H[level][iterates[level]]->keys.find(s[i]) ==
        H[level][iterates[level]]->keys.end())
        continue;

#ifdef PROF
      memory_load_count++;
#endif
      //vector<int> buf;
      auto buf = H[level][iterates[level]]->getbuf(H[level][iterates[level]]->keys[s[i]], ncols);

#ifdef PROF
      memory_load_size += size * sizeof(int);
#endif

      int j = iterates[level];

      int type1 = s[0];
      int t = 0;
      std::array<std::pair<const size_t*, size_t>, ncols_left - 2> nbv;
      for (int li = 1; li < s.size(); li++) {
        if (li == i) continue;
        nbv[t++] = (g.get_neighbors(s[li]));
      }

      auto it1 = buf.begin();

      while (true) {
        size_t length = it1.buffer_size / ncols;


        for (size_t z = 0; z < length; z++) {
          //size_t z = indices[zi];
          // cout << z << " " << z1 << endl;
          int type2 = it1.buffer[z * ncols];

          //std::cout << (length * (*sampling_weights[level])[type2] / (tot_weight * sampling_param[level])) << std::endl;


          const int* it_buf = it1.buffer + z * ncols;

          auto [pr, skip] = sampler.smp_prob({ s.data(), s.size() }, { it_buf, ncols }, i - 1, j, level, tid);
          if (skip) continue;

          sampling_probs[tid][level] = pr;

          std::array<int, ncols_left + ncols - 2> value;

          //auto pat_buf = Pattern::get_labels(g, it_buf, ncols, pats1[type2]);

          auto key_pat_vec =
            get_connectivity<has_labels, edge_induced, ncols, s.size(), value.size()>(nbv, it_buf, s.data(), j + 1, i,
              pats1[type2], pat, value, sgl3);

          for (auto& kp : key_pat_vec) {
            auto& key2 = kp.first;
            auto& ptt = kp.second;
            if (key2 == (1 << 31)) {
              continue;
            }

            key_type key{};

            unsigned int ta = ((i - 1) * (ncols - 1) + j);
            key[0] = (type1);
            key[1] = (type2);
            key[2] = (ta);
            key[3] = (key2);

            if constexpr (has_labels) {
              int t = 4;
              for (int xi = 0; xi < ptt->nn; xi++) {
                key[t++] = ptt->vertex_label[xi];
              }
            }
            {
              auto itt = qp_idx[tid][level].find(key);
              if (itt == qp_idx[tid][level].end()) {
                //  if (qp_idx[tid][level].size() > MAX_NUM_QPATTERN_PER_THREAD) {
                 // std::cerr << "qp_idx flushed" << std::endl;
                  // qp_idx[tid][level].clear();
                   //}
                qp_idx[tid][level][key] = qp_count[tid][level - 1].size();
                value[0] = qp_count[tid][level - 1].size();
                qp_count[tid][level - 1].push_back({ key[0], key[1], i - 1, j });
              }
              else {
                value[0] = itt->second;
              }
            }


            if constexpr (K > 1) {
              for_loop2<pat_agg, has_labels, edge_induced, mni, K - 1, key_type, value.size(), ncols_right...>(sgls, value, ptt, H, level + 1, iterates, res, estimate_counts, sampling_probs, qp2cp,
                qp_count, qp_idx, g, sampler,
                store, sgl3, mni_threshold, need_actual_pattern, actual_patterns, est, adaptive_sampling, query);
            }
            else {
              for_loop2_end<pat_agg, mni, value.size()>(g, value, ptt, level + 1, res, estimate_counts, sampling_probs, est, qp2cp,
                store, mni_threshold, need_actual_pattern, actual_patterns, adaptive_sampling, query, sampler);
            }
          }
        }
        if (!it1.has_next) break;
        it1.next();
      }
    }
  }

  template <bool pat_agg, bool has_labels, bool edge_induced, bool mni, int K, typename key_type, size_t ncols1, size_t ncols2, size_t... ncols>
  void for_loop1(const std::vector<SGList>& sgls,
    const std::vector<std::vector<std::shared_ptr<db::MyKV<int>>>>& H,
    std::vector<int>& iterates, int level, std::vector<SGList>& res,
    std::vector<std::map<int, typename std::conditional<mni, std::tuple<std::string, std::vector<std::vector<unsigned>>, std::vector<unsigned>>, std::string>::type>>& qp2cp, std::vector<std::vector<std::vector<std::array<int, 4>>>>& qp_count,
    std::vector<std::vector<std::map<key_type, int>>>& qp_idx,
    const graph::Graph& g, Sampler& sampler, bool store, const std::pair<std::unordered_set<unsigned long>, std::unordered_set<unsigned long>>& sgl3, double mni_threshold, std::vector<std::map<std::string, double>>& estimate_counts, std::vector<std::vector<double>>& sampling_probs, bool need_actual_pattern, std::map<std::shared_ptr<Pattern>, size_t, cmpByPattern>& actual_patterns, bool adaptive_sampling, bool est, Query& query) {
    if (level < H.size()) {
      for (int i = 0; i < H[level].size(); i++) {
        iterates.push_back(i);
        for_loop1<pat_agg, has_labels, edge_induced, mni, K, key_type, ncols1, ncols2, ncols...>(sgls, H, iterates, level + 1, res, qp2cp, qp_count, qp_idx, g,
          sampler, store, sgl3, mni_threshold, estimate_counts, sampling_probs, need_actual_pattern, actual_patterns, adaptive_sampling, est, query);
        iterates.pop_back();
      }
    }
    else {
      //euler::util::print_vec(iterates);

      int i = iterates[0];
      int j = iterates[1];

      auto& tab1 = H[0][i];
      auto& tab2 = H[1][j];

      auto& pats1 = sgls[0].patterns;
      auto& pats2 = sgls[1].patterns;

      // cout << "tab size: " << tab1->keys.size() << "\t" << tab2->keys.size() <<
      // endl;

      //cout << tab1->keys.size() << endl;
      //size_t key_idx = 0;

      for (auto bi = tab1->keys.begin(); bi != tab1->keys.end(); bi++) {
        //cout << key_idx++ << endl;
#ifdef PROF
        memory_load_count++;
#endif
        //vector<int> reordered_d1;
        //tab1->getbuf(reordered_d1, get<0>(bi->second));
        auto reordered_d1 = tab1->getbuf(bi->second, ncols1);
        auto it_tab2 = tab2->keys.find(bi->first);

#ifdef PROF
        memory_load_size += size1 * sizeof(int);
#endif

        if (it_tab2 != tab2->keys.end()) {
#ifdef PROF
          memory_load_count++;
#endif
          // vector<int> reordered_d2;
          //tab2->getbuf(reordered_d2, get<0>(it_tab2->second));
          auto reordered_d2 = tab2->getbuf(it_tab2->second, ncols2);
#ifdef PROF
          memory_load_size += size2 * sizeof(int);
#endif


          auto it1 = reordered_d1.begin();


          while (true) {
            size_t length1 = it1.buffer_size / ncols1;

#pragma omp parallel for num_threads(_Nthreads) 
            for (size_t z = 0; z < length1; z++) {
              int type1 = it1.buffer[z * ncols1];
              const int* it_d1 = it1.buffer + z * ncols1;




              int qret = query(g, { it_d1, ncols1 }, pats1[type1], 0);
              if (qret < 0) continue;

              int tid = omp_get_thread_num();


              auto [pr1, skip1] = sampler.smp_prob({ nullptr,0 }, { it_d1, ncols1 }, -1, i, 0, tid);

              if (skip1) continue;


              sampling_probs[tid][0] = pr1;


              std::array<std::pair<const size_t*, size_t>, ncols1 - 2> nbv;
              int t = 0;
              for (int li = 1; li < ncols1; li++) {
                if (li == i + 1) continue;
                nbv[t++] = g.get_neighbors(it_d1[li]);
              }

              auto it2 = reordered_d2.begin();


              while (true) {
                size_t length2 = it2.buffer_size / ncols2;

                for (size_t z1 = 0; z1 < length2; z1++) {
                  int type2 = it2.buffer[z1 * ncols2];

                  std::array<int, ncols1 + ncols2 - 2> value;

                  const int* it_d2 = it2.buffer + z1 * ncols2;

                  auto [pr2, skip2] = sampler.smp_prob({ it_d1, ncols1 }, { it_d2, ncols2 }, i, j, 1, tid);
                  if (skip2) continue;

                  sampling_probs[tid][1] = pr2;

                  auto key_pat_vec = get_connectivity<has_labels, edge_induced, ncols2, ncols1, value.size()>(
                    nbv, it_d2, it_d1, j + 1,
                    i + 1, pats2[type2], pats1[type1], value, sgl3);

                  for (auto& kp : key_pat_vec) {
                    auto& key2 = kp.first;
                    auto& ptt = kp.second;

                    if (key2 == (1 << 31)) {
                      continue;
                    }

                    key_type key{};
                    // key.reserve(4);

                    unsigned int ta = (i * (ncols2 - 1) + j);
                    key[0] = type1;
                    key[1] = type2;
                    key[2] = ta;
                    key[3] = key2;

                    if constexpr (has_labels) {
                      int t = 4;
                      for (int xi = 0; xi < ptt->nn; xi++) {
                        key[t++] = ptt->vertex_label[xi];
                      }
                    }

                    // now we get the key for quick pattern: if no label is
                    // considered, the key is simply the topological pattern,
                    // otherwise, the topological pattern appended with the vertex
                    // labels

                    // quick pattern: qp_tmp only stores one value

                    auto itt = qp_idx[tid][0].find(key);
                    if (itt == qp_idx[tid][0].end()) {

                      qp_idx[tid][0][key] = qp_count[tid][0].size();
                      value[0] = qp_count[tid][0].size();
                      qp_count[tid][0].push_back({ key[0], key[1], i, j });

                    }
                    else {
                      value[0] = itt->second;
                    }

                    // if(omp_get_thread_num() == 0) t_for_loop.start();
                    if constexpr (K > 2) {

                      for_loop2<pat_agg, has_labels, edge_induced, mni, K - 2, key_type, value.size(), ncols...>(sgls, value, ptt, H, 2, iterates, res, estimate_counts, sampling_probs, qp2cp,
                        qp_count, qp_idx, g, sampler, store, sgl3, mni_threshold, need_actual_pattern, actual_patterns, est, adaptive_sampling, query);
                    }
                    else {
                      for_loop2_end<pat_agg, mni, value.size()>(g, value, ptt, 2, res, estimate_counts, sampling_probs, est, qp2cp,
                        store, mni_threshold, need_actual_pattern, actual_patterns, adaptive_sampling, query, sampler);
                    }
                  }
                }

                //std::cout << it2.has_next << std::endl;
                if (!it2.has_next) break;
                it2.next();
              }
            }
            //std::cout << it1.has_next << std::endl;
            if (!it1.has_next) break;
            it1.next();
          }
        }
      }
    }
  }

  template<int...>
  struct sum;

  template<int s>
  struct sum<s> {
    enum { value = s - 2 };
  };

  template<int s, int... others>
  struct sum<s, others...> {
    enum { value = s - 2 + sum<others...>::value };
  };

  template<int s1, int s2, int... others>
  struct sum<s1, s2, others...> {
    enum { value = s1 - 2 + sum<s2, others...>::value };
  };


  template <bool pat_agg, bool has_labels, bool edge_induced, bool mni, int K, size_t ncols1, size_t ncols2, size_t... ncols>
  std::tuple<SGList, std::map<std::string, double>> join(const graph::Graph& g, const std::vector<std::vector<std::shared_ptr<db::MyKV<int>>>>& H, const std::vector<SGList>& sgls, bool store, Sampler& sampler = default_sampler, double mni_threshold = -1, bool need_actual_pattern = false, bool est = false, bool adaptive_sampling = false, const std::pair<std::unordered_set<unsigned long>, std::unordered_set<unsigned long>>& sgl3 = join_dummy1, Query& query = default_query) {
    assert((!mni && mni_threshold == -1) || (mni && mni_threshold >= 0));
    int res_size = 2;
    for (auto& d : sgls) {
      if (d.patterns.empty()) return { SGList(), std::map<std::string, double>() };
      int n = d.sgl->ncols - 1;
      res_size += n - 1;
    }

    /*
        for (auto& p : sampling_weights) {
          util::print_vec(*p);
        }*/

    std::vector<SGList> res(_Nthreads);
    for (int i = 0; i < _Nthreads; i++)
      res[i].sgl = std::make_shared<db::MyKV<std::string>>(res_size);
    std::vector<int> iterates;

    //std::cout << sum<ncols1, ncols2, ncols...>::value << std::endl;
    typedef typename std::conditional<has_labels, std::array<int, sum<ncols1, ncols2, ncols...>::value + 5>, std::array<int, 4>>::type key_type;

    std::vector<std::map<int, typename std::conditional<mni, std::tuple<std::string, std::vector<std::vector<unsigned>>, std::vector<unsigned>>, std::string>::type>> qp2cp(_Nthreads);
    std::vector<std::vector<std::vector<std::array<int, 4>>>> qp_count(_Nthreads);
    for (auto& q : qp_count) q.resize(sgls.size() - 1);
    std::vector<std::vector<std::map<key_type, int>>> qp_idx(_Nthreads);
    for (auto& q : qp_idx) q.resize(sgls.size());

    std::vector<std::map<std::string, double>> estimate_counts(_Nthreads);
    std::vector<std::vector<double>> sampling_probs;
    for (int i = 0; i < _Nthreads; i++) {
      sampling_probs.push_back(std::vector<double>(K, 1));
    }

    std::map<std::shared_ptr<Pattern>, size_t, cmpByPattern> actual_patterns;

    for_loop1<pat_agg, has_labels, edge_induced, mni, K, key_type, ncols1, ncols2, ncols...>(sgls, H, iterates, 0, res, qp2cp, qp_count, qp_idx, g,
      sampler, store, sgl3, mni_threshold, estimate_counts, sampling_probs, need_actual_pattern, actual_patterns, adaptive_sampling, est, query);

    //std::cout << "exploration space size: " << exploration_space_size << std::endl;

#ifdef PROF
    //cout << "auto test time: " << t_list_left.get() << endl;
    //cout << "get pattern time: " << t_get_pattern.get() << endl;
    //cout << "iso test time: " << t_iso_check.get() << endl;
    //cout << "for loop time: " << t_for_loop.get() << endl;
    cout << "memory load count: " << memory_load_count << endl;
    cout << "memory load size: " << memory_load_size << endl;
    cout << "iso test count: " << iso_test_count << endl;
    cout << "auto test count: " << auto_test_count << endl;

#endif
    if (adaptive_sampling) {
      for (int i = 0; i < res.size(); i++) {
        res[i].get_quick_pattern_path(qp_count[i]);
      }
    }

    for (int i = 1; i < res.size(); i++) {
      res[0].combine(res[i], mni, store, adaptive_sampling);
    }

    if (need_actual_pattern) {
      res[0].patterns.resize(actual_patterns.size());
      for (auto& [key, value] : actual_patterns) {
        res[0].patterns[value] = key;
      }
    }

    if (est) {
      for (int i = 1; i < _Nthreads; i++) {
        for (auto& [k, v] : estimate_counts[i]) {
          auto it = estimate_counts[0].find(k);
          if (it != estimate_counts[0].end()) {
            it->second += v;
          }
          else {
            estimate_counts[0].insert({ k, v });
          }
        }
      }
    }

    return { res[0], estimate_counts[0] };
  }
}  // namespace euler::pattern_mining
