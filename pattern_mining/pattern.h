#pragma once

#include <dlfcn.h>
#include <sys/time.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "db/db.h"
#include "graph/graph.h"
#include "util.h"

namespace euler::pattern_mining {
  struct Pattern {
    size_t nn;
    size_t ne;
    std::vector<std::set<int>> adj_list;
    std::vector<int> vertex_label;
    bool use_vertex_label;

    Pattern() : nn(0), ne(0), use_vertex_label(false) {}

    Pattern(const Pattern& pat) : nn(pat.nn), ne(pat.ne), vertex_label(pat.vertex_label), use_vertex_label(pat.use_vertex_label) {
      for (auto& s : pat.adj_list) {
        adj_list.push_back(s);
      }
    }

    Pattern(const std::vector<std::pair<int, int>>& pat) {
      use_vertex_label = false;
      std::set<int> nodes;
      ne = pat.size();  // assume there is no duplicate edge in pat
      for (auto& p : pat) {
        nodes.insert(p.first);
        nodes.insert(p.second);
      }
      nn = nodes.size();
      adj_list.resize(nn);
      vertex_label.resize(nn, 0);
      for (auto& p : pat) {
        adj_list[p.first].insert(p.second);
        adj_list[p.second].insert(p.first);
      }
    }

    bool operator<(const Pattern& pp) const {
      return nn < pp.nn || ne < pp.ne || (use_vertex_label && pp.use_vertex_label && vertex_label < pp.vertex_label);
    }

    void clear() {
      nn = 0; ne = 0;
      adj_list.clear();
      vertex_label.clear();
      use_vertex_label = false;
    }

    void set(const Pattern& pat) {
      nn = pat.nn; ne = pat.ne; use_vertex_label = pat.use_vertex_label;
      vertex_label.clear();
      for (auto& s : pat.vertex_label) {
        vertex_label.push_back(s);
      }
      adj_list.clear();
      for (auto& s : pat.adj_list) {
        adj_list.push_back(s);
      }
    }

    void enable_label() {
      use_vertex_label = true;
      vertex_label.clear();
    }

    void combine(const Pattern& p, int j1, int j2) {
      if (p.use_vertex_label != use_vertex_label) {
        std::cerr << "not compatible patterns" << std::endl;
      }
      for (int i = 0; i < p.vertex_label.size(); i++) {
        if (i != j2) vertex_label.push_back(p.vertex_label[i]);
      }
      int ec = 0;
      for (int i = 0; i < p.adj_list.size(); i++) {
        if (i != j2) {
          std::set<int> ts;
          for (int t : p.adj_list[i]) {
            ec++;
            if (t < j2) {
              ts.insert(t + nn);
            }
            else if (t > j2) {
              ts.insert(t - 1 + nn);
            }
            else {
              ts.insert(j1);
            }
          }
          adj_list.push_back(ts);
        }
        else {
          for (int t : p.adj_list[i]) {
            ec++;
            if (t < j2)
              adj_list[j1].insert(t + nn);
            else if (t > j2)
              adj_list[j1].insert(t - 1 + nn);
          }
        }
      }

      nn += p.nn - 1;
      ne += ec / 2;
      assert(adj_list.size() == nn);
    }

    void add_label(int lbl) { vertex_label.push_back(lbl); }

    void add_node(int lbl) {
      if (use_vertex_label == false) {
        std::cerr << "no label allowed" << std::endl;
        exit(-1);
      }
      vertex_label.push_back(lbl);
      nn = vertex_label.size();
    }

    void add_edge(int i, int j) {
      if (adj_list.size() != nn) adj_list.resize(nn);
      adj_list[i].insert(j);
      adj_list[j].insert(i);
      ne++;
    }


    std::shared_ptr<Pattern> permute(const std::vector<int>& perm) {
      auto pat = std::make_shared<Pattern>();
      pat->enable_label();
      std::vector<int> rperm(perm.size());

      for (int i = 0; i < perm.size(); i++) {
        rperm[perm[i]] = i;
      }

      for (int i = 0; i < nn; i++) {
        pat->add_node(vertex_label[rperm[i]]);
      }
      for (int i = 0; i < adj_list.size(); i++) {
        for (int j : adj_list[i]) {
          pat->add_edge(perm[i], perm[j]);
        }
      }
      return pat;
    }


    std::string dfs_coding() const;

    std::tuple<std::string, std::vector<std::vector<unsigned>>, std::vector<unsigned>> canonical_form() const;

    //bool is_separable(int sep) const;

    //bool sep_test(std::vector<int> &visited, int idx, int n) const;

    bool is_connected(std::vector<int> visited) const;

    static std::shared_ptr<Pattern> get_pattern(const graph::Graph& g,
      const std::vector<int>& buf,
      int ncols,
      bool with_label = true);

    static std::shared_ptr<Pattern> get_pattern(const graph::Graph& g,
      const int* buf,
      int ncols,
      bool with_label = true);

    static std::shared_ptr<Pattern> get_labels(const graph::Graph& g,
      const int* buf,
      int ncols, const std::shared_ptr<Pattern> unlabeled_pat);

    // static std::map<int, std::shared_ptr<Pattern>> set_pattern(const Graph &g,
    // std::shared_ptr<MyKV<std::string>> &d);

    bool connected(int x, int y) const {
      return (adj_list[x].find(y) != adj_list[x].end());
    }

    //std::string backtracking(size_t idx, bool store_data, bool remove_sym,
      //                       bool labeled) const;

    void print() const {
      std::cout << " ======== pat start =========" << std::endl;
      std::cout << "has label: " << use_vertex_label << std::endl;
      std::cout << "nn: " << nn << "\t ne: " << ne << std::endl;
      for (int i = 0; i < nn; i++) {
        std::cout << "v " << i << " " << vertex_label[i] << std::endl;
      }
      for (int i = 0; i < nn; i++) {
        for (int j : adj_list[i]) {
          if (i < j) std::cout << "e " << i << " " << j << std::endl;
        }
      }
      std::cout << "======== pat end ============" << std::endl;
    }

    std::string serialize() const {
      std::string res = std::to_string(nn) + " " + std::to_string(ne) + "\n";
      for (int i = 0; i < nn; i++) {
        res += std::to_string(vertex_label[i]) + "\n";
      }
      int ec = 0;
      for (int i = 0; i < nn; i++) {
        for (int j : adj_list[i]) {
          if (i < j) {
            res += std::to_string(i) + " " + std::to_string(j) + "\n";
            ec++;
          }
        }
      }
      assert(ec == ne);
      res.pop_back();
      return res;
    }
  };

}  // namespace euler::mining
