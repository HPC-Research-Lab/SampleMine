#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/isomorphism.hpp>
#include <fstream>
#include <iostream>
#include <list>
#include <string>

#include "pattern.h"
namespace euler::pattern_mining {

class PatListing {
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
      graph_t;

  std::vector<std::pair<int, int>> idx2edge;
  std::list<graph_t> pat_list;

  void dfs(std::vector<int> &visited, int x, graph_t &g) {
    if (!visited[x]) {
      visited[x] = 1;
      typename boost::graph_traits<graph_t>::out_edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = out_edges(x, g); ei != ei_end; ++ei) {
        int t = target(*ei, g);
        dfs(visited, t, g);
      }
    }
  }

  bool is_connected(graph_t &g) {
    std::vector<int> visited(num_vertices(g));
    dfs(visited, 0, g);
    for (int i = 0; i < visited.size(); i++) {
      if (visited[i] == 0) return false;
    }
    return true;
  }

  void enumerate(graph_t &g, int t) {
    if (t >= 0) {
      graph_t g1(g);
      add_edge(idx2edge[t].first, idx2edge[t].second, g1);
      enumerate(g, t - 1);
      enumerate(g1, t - 1);
    } else {
      if (is_connected(g)) pat_list.push_back(g);
    }
  }

 public:
  std::vector<std::vector<std::pair<int, int>>> pattern_listing(int n) {
    for (int i = 0; i < n - 1; i++) {
      for (int j = i + 1; j < n; j++) {
        idx2edge.push_back(std::make_pair(i, j));
      }
    }

    graph_t g(n);

    enumerate(g, idx2edge.size() - 1);

    for (std::list<graph_t>::iterator p1 = pat_list.begin();
         p1 != pat_list.end(); ++p1) {
      std::list<graph_t>::iterator p2 = p1;
      if (p2 == pat_list.end()) break;
      ++p2;
      while (p2 != pat_list.end()) {
        std::vector<graph_t::vertex_descriptor> f(n);
        if (boost::isomorphism(*p1, *p2, boost::isomorphism_map(f.data()))) {
          pat_list.erase(p2++);
        } else {
          p2++;
        }
      }
    }

    std::vector<std::vector<std::pair<int, int>>> res;
    for (auto &g : pat_list) {
      res.push_back(std::vector<std::pair<int, int>>());
      typename boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        res.back().push_back(
            std::make_pair(boost::source(*ei, g), boost::target(*ei, g)));
      }
    }
    return res;
  }

  static std::vector<std::vector<std::pair<int, int>>> spanning_tree_listing(int n);

  static std::vector<std::shared_ptr<Pattern>> make_pattern(
      const std::vector<std::vector<std::pair<int, int>>> &pat) {
    std::vector<std::shared_ptr<Pattern>> res(pat.size());
    for (int i = 0; i < pat.size(); i++) {
      res[i] = std::make_shared<Pattern>(pat[i]);
    }
    return res;
  }
};
}