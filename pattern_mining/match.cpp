#include "match.h"

#include <dlfcn.h>
#include <fcntl.h>
#include <omp.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <random>  // mt19937 and uniform_int_distribution
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

#include "util.h"

namespace fs = std::experimental::filesystem;

namespace euler::pattern_mining {

  const static int _Nthreads = atoi(getenv("OMP_NUM_THREADS"));

  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<> dis(0, 1);

  static double random_number() {
    return dis(gen);
  }




  bool nested_for_loop(
    const graph::Graph& g, const std::set<std::pair<int, int>>& L,
    std::vector<std::shared_ptr<db::MyKV<std::string>>> data, int level, int nn,
    std::vector<std::vector<int>>& v,
    std::vector<std::vector<std::vector<int>>>& vertex_set,
    const std::vector<int>& start_set, std::atomic<size_t>& count,
    const std::shared_ptr<Pattern> pat, const std::vector<int>& node_order,
    std::map<std::shared_ptr<Pattern>, std::tuple<size_t, std::string, std::vector<std::vector<unsigned int>>, std::vector<unsigned int>>, cmpByPattern>& actual_patterns,
    std::string& key, bool store_data, bool output_labeled,
    std::vector<std::vector<std::pair<bool, std::vector<int>>>>& z,
    bool edge_induced, std::vector<size_t>& count_per_vertex, double sampling_threshold, double mni, bool testing, bool pattern_labeled) {
    if (level == 0) {
#pragma omp parallel for num_threads(_Nthreads)
      for (int i = 0; i < start_set.size(); i++) {
        int tid = omp_get_thread_num();
        v[tid][level + 1] = start_set[i];

        if (pattern_labeled && g.get_vertex_label(v[tid][level + 1]) != pat->vertex_label[node_order[level]]) continue;

        auto& vertex_set_next = vertex_set[tid][level + 1];
        vertex_set_next.clear();
        int set_size = 0;

        int first = -1;
        for (int j = level; j >= 0; j--) {
          if (pat->connected(node_order[j], node_order[level + 1])) {
            g.get_neighbors(v[tid][j + 1], vertex_set_next);
            first = j;
            break;
          }
        }
        int pivot = -1;

        for (int j = level; j >= 0; j--) {
          if (L.find(std::make_pair(j, level + 1)) != L.end()) {
            pivot = j + 1;
            break;
          }
        }
        if (!testing && pivot != -1) {
          int pos = std::lower_bound(vertex_set_next.begin(),
            vertex_set_next.end(), v[tid][pivot]) -
            vertex_set_next.begin();
          vertex_set_next.resize(pos);
        }
        set_size = vertex_set_next.size();
        for (int j = level; j >= 0; j--) {
          if (j == first)
            continue;
          std::vector<int> vertex_set_vj;
          g.get_neighbors(v[tid][j + 1], vertex_set_vj);
          int set_size_j = vertex_set_vj.size();
          if (pat->connected(node_order[j], node_order[level + 1])) {
            int pi = 0, pj = 0, pos = 0;
            while (pi != set_size && pj != set_size_j) {
              if (vertex_set_next[pi] < vertex_set_vj[pj])
                pi++;
              else if (vertex_set_next[pi] > vertex_set_vj[pj])
                pj++;
              else {
                vertex_set_next[pos++] = vertex_set_next[pi];
                pi++;
                pj++;
              }
            }
            set_size = pos;
          }
          else if (!testing) {
            int pi = 0, pj = 0, pos = 0, pos2 = 0;
            while (pi != set_size && pj != set_size_j) {
              if (vertex_set_next[pi] < vertex_set_vj[pj])
                pi++;
              else if (vertex_set_next[pi] > vertex_set_vj[pj])
                pj++;
              else {
                vertex_set_vj[pos++] = pi;
                pi++;
                pj++;
              }
            }
            set_size_j = pos;
            pos = 0;
            for (pi = 0; pi < set_size; pi++) {
              if (pos >= set_size_j) {
                vertex_set_next[pos2++] = vertex_set_next[pi];
              }
              else if (pi == vertex_set_vj[pos]) {
                pos++;
                continue;
              }
              else {
                vertex_set_next[pos2++] = vertex_set_next[pi];
              }
            }
            set_size = pos2;
          }
        }

        // print_vec(vertex_set_next);

        std::vector<int> vertex_set_vj;
        vertex_set_vj.resize(level + 1);
        int set_size_j = 0;
        for (int tz = 0; tz < level + 1; tz++) {
          int tzz = 0;
          while (tzz < set_size_j&& vertex_set_vj[tzz] <= v[tid][tz + 1])
            tzz++;
          for (int tx = set_size_j; tx > tzz; tx--)
            vertex_set_vj[tx] = vertex_set_vj[tx - 1];
          vertex_set_vj[tzz] = v[tid][tz + 1];
          set_size_j++;
        }

        if (set_size_j > 0) {
          int pi = 0, pj = 0, pos = 0, pos2 = 0;
          while (pi != set_size && pj != set_size_j) {
            if (vertex_set_next[pi] < vertex_set_vj[pj])
              pi++;
            else if (vertex_set_next[pi] > vertex_set_vj[pj])
              pj++;
            else {
              vertex_set_vj[pos++] = pi;
              pi++;
              pj++;
            }
          }
          set_size_j = pos;
          pos = 0;
          for (pi = 0; pi < set_size; pi++) {
            if (pos >= set_size_j) {
              vertex_set_next[pos2++] = vertex_set_next[pi];
            }
            else if (pi == vertex_set_vj[pos]) {
              pos++;
              continue;
            }
            else {
              vertex_set_next[pos2++] = vertex_set_next[pi];
            }
          }
          set_size = pos2;
        }

        vertex_set_next.resize(set_size);

        // print_vec(vertex_set_next);

        auto& z_next = z[tid][level + 1].second;
        z_next.clear();

        first = -1;
        if (level + 2 < nn) {
          for (int j = level; j >= 0; j--) {
            if (pat->connected(node_order[j], node_order[level + 2])) {
              g.get_neighbors(v[tid][j + 1], z_next);
              first = j;
              break;
            }
          }
        }
        if (first != -1) {
          z[tid][level + 1].first = true;
          int set_size_z = z_next.size();
          for (int j = level; j >= 0; j--) {
            if (j == first)
              continue;
            std::vector<int> vertex_set_vj;
            g.get_neighbors(v[tid][j + 1], vertex_set_vj);
            int set_size_j = vertex_set_vj.size();
            if (pat->connected(node_order[j], node_order[level + 2])) {
              int pi = 0, pj = 0, pos = 0;
              while (pi != set_size_z && pj != set_size_j) {
                if (z_next[pi] < vertex_set_vj[pj])
                  pi++;
                else if (z_next[pi] > vertex_set_vj[pj])
                  pj++;
                else {
                  z_next[pos++] = z_next[pi];
                  pi++;
                  pj++;
                }
              }
              set_size_z = pos;
            }
            else if (!testing) {
              int pi = 0, pj = 0, pos = 0, pos2 = 0;
              while (pi != set_size_z && pj != set_size_j) {
                if (z_next[pi] < vertex_set_vj[pj])
                  pi++;
                else if (z_next[pi] > vertex_set_vj[pj])
                  pj++;
                else {
                  vertex_set_vj[pos++] = pi;
                  pi++;
                  pj++;
                }
              }
              set_size_j = pos;
              pos = 0;
              for (pi = 0; pi < set_size_z; pi++) {
                if (pos >= set_size_j) {
                  z_next[pos2++] = z_next[pi];
                }
                else if (pi == vertex_set_vj[pos]) {
                  pos++;
                  continue;
                }
                else {
                  z_next[pos2++] = z_next[pi];
                }
              }
              set_size_z = pos2;
            }
          }
          z_next.resize(set_size_z);
        }
        else {
          z[tid][level + 1].first = false;
        }
        //  print_vec(vertex_set_next);
        // std::cout << "vs: " << vertex_set[tid][level+1].size() << std::endl;

        nested_for_loop(g, L, data, level + 1, nn, v, vertex_set, start_set,
          count, pat, node_order, actual_patterns, key, store_data, output_labeled, z, edge_induced, count_per_vertex, sampling_threshold, mni, testing, pattern_labeled);
      }
      return true;
    }
    else if (level < nn - 1) {
      int tid = omp_get_thread_num();
      if (sampling_threshold > 0) {
        std::random_shuffle(vertex_set[tid][level].begin(), vertex_set[tid][level].end());
      }
      for (int i = 0; i < vertex_set[tid][level].size(); i++) {
        v[tid][level + 1] = vertex_set[tid][level][i];

        if (pattern_labeled && g.get_vertex_label(v[tid][level + 1]) != pat->vertex_label[node_order[level]]) continue;


        auto& vertex_set_next = vertex_set[tid][level + 1];
        vertex_set_next.clear();
        int set_size = 0;

        if (z[tid][level].first == false) {
          int first = -1;
          for (int j = level; j >= 0; j--) {
            if (pat->connected(node_order[j], node_order[level + 1])) {
              g.get_neighbors(v[tid][j + 1], vertex_set_next);
              first = j;
              break;
            }
          }
          int pivot = -1;

          for (int j = level; j >= 0; j--) {
            if (L.find(std::make_pair(j, level + 1)) != L.end()) {
              pivot = j + 1;
              break;
            }
          }
          if (!testing && pivot != -1) {
            int pos = std::lower_bound(vertex_set_next.begin(),
              vertex_set_next.end(), v[tid][pivot]) -
              vertex_set_next.begin();
            vertex_set_next.resize(pos);
          }
          set_size = vertex_set_next.size();
          for (int j = level; j >= 0; j--) {
            if (j == first)
              continue;
            std::vector<int> vertex_set_vj;
            g.get_neighbors(v[tid][j + 1], vertex_set_vj);
            int set_size_j = vertex_set_vj.size();
            if (pat->connected(node_order[j], node_order[level + 1])) {
              int pi = 0, pj = 0, pos = 0;
              while (pi != set_size && pj != set_size_j) {
                if (vertex_set_next[pi] < vertex_set_vj[pj])
                  pi++;
                else if (vertex_set_next[pi] > vertex_set_vj[pj])
                  pj++;
                else {
                  vertex_set_next[pos++] = vertex_set_next[pi];
                  pi++;
                  pj++;
                }
              }
              set_size = pos;
            }
            else if (!testing) {
              int pi = 0, pj = 0, pos = 0, pos2 = 0;
              while (pi != set_size && pj != set_size_j) {
                if (vertex_set_next[pi] < vertex_set_vj[pj])
                  pi++;
                else if (vertex_set_next[pi] > vertex_set_vj[pj])
                  pj++;
                else {
                  vertex_set_vj[pos++] = pi;
                  pi++;
                  pj++;
                }
              }
              set_size_j = pos;
              pos = 0;
              for (pi = 0; pi < set_size; pi++) {
                if (pos >= set_size_j) {
                  vertex_set_next[pos2++] = vertex_set_next[pi];
                }
                else if (pi == vertex_set_vj[pos]) {
                  pos++;
                  continue;
                }
                else {
                  vertex_set_next[pos2++] = vertex_set_next[pi];
                }
              }
              set_size = pos2;
            }
          }
        }
        else {
          vertex_set_next = z[tid][level].second;
          int pivot = -1;

          for (int j = level; j >= 0; j--) {
            if (L.find(std::make_pair(j, level + 1)) != L.end()) {
              pivot = j + 1;
              break;
            }
          }
          if (!testing && pivot != -1) {
            int pos = std::lower_bound(vertex_set_next.begin(),
              vertex_set_next.end(), v[tid][pivot]) -
              vertex_set_next.begin();
            vertex_set_next.resize(pos);
          }

          set_size = vertex_set_next.size();

          std::vector<int> vertex_set_vj;
          g.get_neighbors(v[tid][level + 1], vertex_set_vj);
          int set_size_j = vertex_set_vj.size();
          if (pat->connected(node_order[level + 1], node_order[level])) {
            int pi = 0, pj = 0, pos = 0;
            while (pi != set_size && pj != set_size_j) {
              if (vertex_set_next[pi] < vertex_set_vj[pj])
                pi++;
              else if (vertex_set_next[pi] > vertex_set_vj[pj])
                pj++;
              else {
                vertex_set_next[pos++] = vertex_set_next[pi];
                pi++;
                pj++;
              }
            }
            set_size = pos;
          }
          else if (!testing) {
            int pi = 0, pj = 0, pos = 0, pos2 = 0;
            while (pi != set_size && pj != set_size_j) {
              if (vertex_set_next[pi] < vertex_set_vj[pj])
                pi++;
              else if (vertex_set_next[pi] > vertex_set_vj[pj])
                pj++;
              else {
                vertex_set_vj[pos++] = pi;
                pi++;
                pj++;
              }
            }
            set_size_j = pos;
            pos = 0;
            for (pi = 0; pi < set_size; pi++) {
              if (pos >= set_size_j) {
                vertex_set_next[pos2++] = vertex_set_next[pi];
              }
              else if (pi == vertex_set_vj[pos]) {
                pos++;
                continue;
              }
              else {
                vertex_set_next[pos2++] = vertex_set_next[pi];
              }
            }
            set_size = pos2;
          }
        }

        // print_vec(vertex_set_next);

        std::vector<int> vertex_set_vj;
        vertex_set_vj.resize(level + 1);
        int set_size_j = 0;
        for (int tz = 0; tz < level + 1; tz++) {
          int tzz = 0;
          while (tzz < set_size_j&& vertex_set_vj[tzz] <= v[tid][tz + 1])
            tzz++;
          for (int tx = set_size_j; tx > tzz; tx--)
            vertex_set_vj[tx] = vertex_set_vj[tx - 1];
          vertex_set_vj[tzz] = v[tid][tz + 1];
          set_size_j++;
        }

        if (set_size_j > 0) {
          int pi = 0, pj = 0, pos = 0, pos2 = 0;
          while (pi != set_size && pj != set_size_j) {
            if (vertex_set_next[pi] < vertex_set_vj[pj])
              pi++;
            else if (vertex_set_next[pi] > vertex_set_vj[pj])
              pj++;
            else {
              vertex_set_vj[pos++] = pi;
              pi++;
              pj++;
            }
          }
          set_size_j = pos;
          pos = 0;
          for (pi = 0; pi < set_size; pi++) {
            if (pos >= set_size_j) {
              vertex_set_next[pos2++] = vertex_set_next[pi];
            }
            else if (pi == vertex_set_vj[pos]) {
              pos++;
              continue;
            }
            else {
              vertex_set_next[pos2++] = vertex_set_next[pi];
            }
          }
          set_size = pos2;
        }

        vertex_set_next.resize(set_size);

        // print_vec(vertex_set_next);

        auto& z_next = z[tid][level + 1].second;
        z_next.clear();

        int first = -1;
        if (level + 2 < nn) {
          for (int j = level; j >= 0; j--) {
            if (pat->connected(node_order[j], node_order[level + 2])) {
              g.get_neighbors(v[tid][j + 1], z_next);
              first = j;
              break;
            }
          }
        }
        if (first != -1) {
          z[tid][level + 1].first = true;
          int set_size_z = z_next.size();
          for (int j = level; j >= 0; j--) {
            if (j == first)
              continue;
            std::vector<int> vertex_set_vj;
            g.get_neighbors(v[tid][j + 1], vertex_set_vj);
            int set_size_j = vertex_set_vj.size();
            if (pat->connected(node_order[j], node_order[level + 2])) {
              int pi = 0, pj = 0, pos = 0;
              while (pi != set_size_z && pj != set_size_j) {
                if (z_next[pi] < vertex_set_vj[pj])
                  pi++;
                else if (z_next[pi] > vertex_set_vj[pj])
                  pj++;
                else {
                  z_next[pos++] = z_next[pi];
                  pi++;
                  pj++;
                }
              }
              set_size_z = pos;
            }
            else if (!testing) {
              int pi = 0, pj = 0, pos = 0, pos2 = 0;
              while (pi != set_size_z && pj != set_size_j) {
                if (z_next[pi] < vertex_set_vj[pj])
                  pi++;
                else if (z_next[pi] > vertex_set_vj[pj])
                  pj++;
                else {
                  vertex_set_vj[pos++] = pi;
                  pi++;
                  pj++;
                }
              }
              set_size_j = pos;
              pos = 0;
              for (pi = 0; pi < set_size_z; pi++) {
                if (pos >= set_size_j) {
                  z_next[pos2++] = z_next[pi];
                }
                else if (pi == vertex_set_vj[pos]) {
                  pos++;
                  continue;
                }
                else {
                  z_next[pos2++] = z_next[pi];
                }
              }
              set_size_z = pos2;
            }
          }
          z_next.resize(set_size_z);
        }
        else {
          z[tid][level + 1].first = false;
        }
        //  print_vec(vertex_set_next);

        bool ret = nested_for_loop(g, L, data, level + 1, nn, v, vertex_set, start_set,
          count, pat, node_order, actual_patterns, key, store_data, output_labeled, z, edge_induced, count_per_vertex, sampling_threshold, mni, testing, pattern_labeled);

        if (!ret) return false;
      }
      return true;
    }
    else {
      int tid = omp_get_thread_num();
      if (sampling_threshold > 0) {
        std::random_shuffle(vertex_set[tid][level].begin(), vertex_set[tid][level].end());
      }
      for (int i = 0; i < vertex_set[tid][level].size(); i++) {
        v[tid][level + 1] = vertex_set[tid][level][i];

        if (pattern_labeled && g.get_vertex_label(v[tid][level + 1]) != pat->vertex_label[node_order[level]]) continue;

        if (level > 0) {
          bool to_break = false;

          for (int j = level - 1; j >= 0; j--) {
            if (L.find(std::make_pair(j, level)) != L.end()) {
              if (v[tid][j + 1] < v[tid][level + 1]) {
                to_break = true;
              }
              break;
            }
          }
          if (!testing && to_break)
            break;
        }
        // std::cout << tid << " " << level << std::endl;
        //util::print_vec(v[tid]);

        std::vector<std::shared_ptr<Pattern>> possible_patterns;
        //util::print_vec(v[tid]);


        if (edge_induced) {
          auto patt = Pattern::get_pattern(g, v[tid], nn + 1, false);


          unsigned num_combinations = (1 << patt->ne);

          for (unsigned ii = 0; ii < num_combinations; ii++) {
            std::shared_ptr<Pattern> pttt = std::make_shared<Pattern>(*patt);

            bool connected = true;

            int cc = 0;

            for (int j = 0; j < patt->ne; j++) {
              if (ii & (1 << j)) {
                pttt->remove_edge(j - cc);
                cc++;
                if (!pttt->connected()) {
                  connected = false;
                  break;
                }
              }
            }
            if (connected) {
              //pttt->print();
              possible_patterns.push_back(pttt);
            }
          }
          //std::cout << "!!!!!!!!!!!" << std::endl;
        }
        else {
          possible_patterns.push_back(pat);
        }

        for (auto ptr : possible_patterns) {

          count++;
          count_per_vertex[v[tid][1]]++;
          if (sampling_threshold > 0 && count_per_vertex[v[tid][1]] >= sampling_threshold) {
            return false;
          }

          if (output_labeled) {
            auto ptt = std::make_shared<Pattern>(*ptr);

            ptt->enable_label();
            for (int k = 1; k < nn + 1; k++) {
              int l = g.get_vertex_label(v[tid][k]);
              ptt->add_label(l);
            }

            // if the matching is labeled, we aggregagate the subgraph with their actual pattern. 
            // The first entry of each subgraph store the pattern id

            //ptt->print();

            //std::cout << std::get<0>(cp) << std::endl;
            //util::print_vec(cp.second);

            std::tuple<std::string, std::vector<std::vector<unsigned int>>, std::vector<unsigned int>> cp;
#pragma omp critical
            {
              auto it_pat = actual_patterns.find(ptt);
              if (it_pat == actual_patterns.end()) {
                cp = ptt->canonical_form();
                actual_patterns[ptt] = std::make_tuple(actual_patterns.size(), std::get<0>(cp), std::get<1>(cp), std::get<2>(cp));
                v[tid][0] = actual_patterns.size() - 1;
              }
              else {
                v[tid][0] = std::get<0>(it_pat->second);
                cp = std::make_tuple(std::get<1>(it_pat->second), std::get<2>(it_pat->second), std::get<3>(it_pat->second));
              }
            }
            data[tid]->merge(std::get<0>(cp), v[tid].data(), v[tid].size() * sizeof(int), store_data, mni, std::get<1>(cp), std::get<2>(cp));
          }
          else {
            // if the matching is unlabeled, we simply aggregate them based on the unlabeled patterns
            // The first entry of each subgaph is the unlabeled pattern id. 
            if (store_data) {
              {
                data[tid]->merge(key, v[tid].data(), v[tid].size() * sizeof(int));
              }
            }
          }
        }
      }
      return true;
    }  // namespace euler::pattern_mining
  }  // namespace euler::pattern_mining

  void permutation(std::vector<std::vector<int>>& all, std::vector<int>& a, int l,
    int r) {
    // Base case
    if (l == r)
      all.push_back(a);
    else {
      // Permutations made
      for (int i = l; i <= r; i++) {
        // Swapping done
        std::swap(a[l], a[i]);
        // Recursion called
        permutation(all, a, l + 1, r);
        // backtrack
        std::swap(a[l], a[i]);
      }
    }
  }



  SGList match(const graph::Graph& g, const PatList& patterns, bool store_data, bool edge_induced, bool output_labeled, double mni, bool testing, bool pattern_labeled, double sampling_threshold) {
    // if mni, the output subgraph has to be labeled
    assert(mni == 0 || output_labeled);

    if (patterns.size() == 0)
      return SGList();

    assert(sampling_threshold >= 0);

    static int tag = 0;
    const int s1 = patterns[0]->nn + 1;

    auto data = std::vector<std::shared_ptr<db::MyKV<std::string>>>();

    for (int i = 0; i < _Nthreads; i++) {
      data.push_back(std::make_shared<db::MyKV<std::string>>(s1));
    }

    std::vector<int> p1;
    for (int i = 0; i < s1 - 1; i++)
      p1.push_back(i);
    std::vector<std::vector<int>> permute;
    permutation(permute, p1, 0, s1 - 2);

    size_t pi = 0;

    std::map<std::shared_ptr<Pattern>, std::tuple<size_t, std::string, std::vector<std::vector<unsigned int>>, std::vector<unsigned int>>, cmpByPattern> actual_patterns;

    for (auto& pat : patterns) {

      // pat->print();

      int root = 0;
      size_t max_degree = pat->adj_list[0].size();
      for (int i = 1; i < pat->nn; i++) {
        if (pat->adj_list[i].size() > max_degree) {
          root = i;
          max_degree = pat->adj_list[i].size();
        }
      }

      std::queue<int> q;
      q.push(root);
      std::vector<int> node_order;
      std::vector<int> visited(pat->nn, 0);
      while (!q.empty()) {
        int a = q.front();
        q.pop();
        if (!visited[a])
          node_order.push_back(a);
        visited[a] = 1;
        for (int b : pat->adj_list[a])
          if (!visited[b])
            q.push(b);
      }
      std::vector<int> order_map(pat->nn);
      for (int i = 0; i < pat->nn; i++)
        order_map[node_order[i]] = i;
      // cout << permute.size() << endl;

      std::vector<std::vector<int>> valid_permute;

      for (auto& pp : permute) {
        std::vector<std::set<int>> adj_tmp(pat->nn);
        for (int i = 0; i < pat->nn; i++) {
          std::set<int> tp;
          for (int j : pat->adj_list[i])
            tp.insert(pp[j]);
          adj_tmp[pp[i]] = tp;
        }
        bool valid = true;
        for (int i = 0; i < pat->nn; i++) {
          if (adj_tmp[i] != pat->adj_list[i]) {
            valid = false;
            break;
          }
        }
        if (valid)
          valid_permute.push_back(pp);
      }

      std::set<std::pair<int, int>> L;
      for (int v : node_order) {
        std::vector<std::vector<int>> stabilized_aut;
        for (auto& x : valid_permute) {
          if (x[v] == v) {
            stabilized_aut.push_back(x);
          }
          else {
            L.insert(std::make_pair(order_map[v], order_map[x[v]]));
          }
        }
        valid_permute = stabilized_aut;
      }

      auto vv = std::vector<std::vector<int>>(_Nthreads);

      for (auto& v : vv) {
        v.resize(pat->nn + 1);
        v[0] = pi;
      }
      auto vertex_set = std::vector<std::vector<std::vector<int>>>(_Nthreads);
      for (auto& v : vertex_set) {
        v.resize(pat->nn);
        for (auto& x : v)
          x.reserve(g.max_degree());
      }
      std::vector<int> start_set;
      for (int i = 0; i < g.num_nodes(); i++)
        start_set.push_back(i);
      // std::cout << start_set.size() << std::endl;
      std::atomic<size_t> count = 0;
      std::string key = std::to_string(pi++);
      auto z = std::vector<std::vector<std::pair<bool, std::vector<int>>>>(_Nthreads);
      for (auto& v : z) {
        v.resize(pat->nn);
        for (auto& x : v)
          x.second.reserve(g.max_degree());
      }

      std::vector<size_t> count_per_vertex(g.num_nodes(), 0);

      nested_for_loop(g, L, data, 0, pat->nn, vv, vertex_set, start_set, count,
        pat, node_order, actual_patterns, key, store_data, output_labeled, z, edge_induced, count_per_vertex, sampling_threshold, mni, testing, pattern_labeled);

      if (!testing)
        std::cout << "count: " << count << std::endl;
    }


    for (int i = 1; i < data.size(); i++) {
      data[0]->combine(*data[i], mni >= 0, store_data, false);
    }

    if (!output_labeled) {
      return SGList(data[0], patterns);
    }
    else {
      PatList labeled_patterns(actual_patterns.size());
      for (auto& [key, value] : actual_patterns) {
        labeled_patterns[std::get<0>(value)] = key;
      }
      return SGList(data[0], labeled_patterns);
    }
  }


  std::pair<std::unordered_set<unsigned long>, std::unordered_set<unsigned long>> get_pattern3(const SGList& sgl) {
    assert(sizeof(unsigned long) == 8);
    std::unordered_set<unsigned long> triangles;
    std::unordered_set<unsigned long> wedges;
    for (auto& s : sgl.patterns) {
      //std::cout << s->ne << std::endl;
      if (s->ne == 3) {
        unsigned long lab = 0;
        int labs[3];
        for (int i = 0; i < 3; i++) {
          int l = s->vertex_label[i];
          assert(l < 63);
          labs[i] = l;
          lab |= (1ULL << l);
        }
        if ((labs[0] == labs[1]) && (labs[0] < labs[2])) lab |= (1ULL << 63);
        if ((labs[0] == labs[2]) && (labs[0] < labs[1])) lab |= (1ULL << 63);
        if ((labs[1] == labs[2]) && (labs[1] < labs[0])) lab |= (1ULL << 63);
        triangles.insert(lab);
      }
      else if (s->ne == 2) {
        unsigned long lab = 0;
        int labs[3];
        int key;
        for (int i = 0; i < 3; i++) {
          int l = s->vertex_label[i];
          assert(l < 61);
          labs[i] = l;
          // std::cout << l << " ";
          lab |= (1ULL << l);
          if (s->adj_list[i].size() == 2) {
            key = l;
          }
        }
        //  std::cout << key << " ";
        int c = 0;
        int e = 0;
        for (int i = 0; i < 3; i++) {
          if (labs[i] > key) c++;
          if (labs[i] == key) e++;
        }
        if (c == 1) lab |= (1ULL << 63);
        if (c == 2) lab |= (1ULL << 62);
        if (c == 0 && e == 2) lab |= (1ULL << 61);

        // std::cout << " : " << std::bitset<64>(lab) << std::endl;

        wedges.insert(lab);
      }
      else {
        exit(-1);
      }
    }
    //std::cout << triangles.size() << " " << wedges.size() << " " << sgl.patterns.size() << std::endl;
    //assert(triangles.size() + wedges.size() == sgl.patterns.size());
    return std::make_pair(triangles, wedges);
  }
}  // namespace euler::pattern_mining