#include "pattern.h"

#include <omp.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <sstream>

#include "bliss/defs.hh"
#include "bliss/graph.hh"
#include "bliss/timer.hh"
#include "bliss/utils.hh"

namespace euler::pattern_mining {
bool Pattern::is_connected(std::vector<int> visited) const {
  int source;
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
    for (int nb : adj_list[cur]) {
      if (!visited[nb]) q.push(nb);
    }
  }
  for (int i = 1; i < visited.size(); i++)
    if (!visited[i]) return false;
  return true;
}
/*
bool Pattern::sep_test(std::vector<int> &visited, int idx, int n) const {
  visited[idx] = 1;
  if (n == 0) {
    std::vector<int> visited_tmp(visited.begin(), visited.end());
    for (int xi = 0; xi < visited_tmp.size(); xi++) {
      if (visited_tmp[xi] == 1) {
        visited_tmp[xi] = 0;
        if (is_connected(visited_tmp)) {
          return true;
        }
        visited_tmp[xi] = 1;
      }
    }
    visited[idx] = 0;
    return false;
  } else {
    for (int xi = 0; xi < visited.size(); xi++) {
      if (visited[xi]) {
        for (int nb : adj_list[xi]) {
          if (!visited[nb]) {
            if (sep_test(visited, nb, n - 1)) return true;
            visited[nb] = 0;
          }
        }
      }
    }
  }
}

bool Pattern::is_separable(int sep) const {
  for (int i = 0; i < nn; i++) {
    std::vector<int> visited(nn, 0);
    if (sep_test(visited, i, sep - 1)) return true;
  }
  return false;
}*/

std::shared_ptr<Pattern> Pattern::get_pattern(const graph::Graph &g,
                                              const std::vector<int> &buf,
                                              int ncols, bool with_label) {
  auto p = std::make_shared<Pattern>();
  p->enable_label();
  for (int j = 1; j < ncols; j++) {
    if (with_label)
      p->add_node(g.get_vertex_label(buf[j]));
    else
      p->add_node(0);
  }

  for (int xi = 1; xi < ncols - 1; xi++) {
    for (int yi = xi + 1; yi < ncols; yi++) {
      if (g.is_neighbor(buf[xi], buf[yi])) {
        p->add_edge(xi - 1, yi - 1);
      } 
    }
  }
  // std::cout << p->dfs_coding() << std::endl;
  return p;
}

std::shared_ptr<Pattern> Pattern::get_pattern(const graph::Graph &g,
                                              const int *buf,
                                              int ncols, bool with_label) {
  auto p = std::make_shared<Pattern>();
  p->enable_label();
  for (int j = 1; j < ncols; j++) {
    if (with_label)
      p->add_node(g.get_vertex_label(buf[j]));
    else
      p->add_node(0);
  }

  for (int xi = 1; xi < ncols - 1; xi++) {
    for (int yi = xi + 1; yi < ncols; yi++) {
      if (g.is_neighbor(buf[xi], buf[yi])) {
        p->add_edge(xi - 1, yi - 1);
      } 
    }
  }
  // std::cout << p->dfs_coding() << std::endl;
  return p;
}


std::shared_ptr<Pattern> Pattern::get_labels(const graph::Graph &g,
                                              const int *buf,
                                              int ncols, const std::shared_ptr<Pattern> pat) {
  auto p = std::make_shared<Pattern>(*pat);
  p->enable_label();
  for (int j = 1; j < ncols; j++) {
      p->add_label(g.get_vertex_label(buf[j]));
  }
  // std::cout << p->dfs_coding() << std::endl;
  return p;
}


std::string Pattern::dfs_coding() const {
  if (!use_vertex_label) {
    std::cerr << "no label" << std::endl;
    exit(-1);
  }
  bliss::Graph g;
  for (int v : vertex_label) {
    g.add_vertex(v);
  }
  for (int i = 0; i < adj_list.size(); i++) {
    for (int j : adj_list[i]) {
      if (i < j) {
        g.add_edge(i, j);
      }
    }
  }
  bliss::Stats stats;
  const unsigned int *cl = g.canonical_form(stats, NULL, stdout);
  bliss::Graph *cf = g.permute(cl);
  std::ostringstream sout;
  cf->write_dimacs(sout);
  delete cf;

  return sout.str();
}

static void
report_aut(void* param, const unsigned int n, const unsigned int* aut)
{
  assert(param);
  fprintf((FILE*)param, "Generator: ");
  bliss::print_permutation((FILE*)param, n, aut, 1);
  fprintf((FILE*)param, "\n");
}


std::pair<std::string, std::vector<int>> Pattern::canonical_form() const {
  if (!use_vertex_label) {
    std::cerr << "no label" << std::endl;
    exit(-1);
  }
  bliss::Graph g;
  for (int v : vertex_label) {
    g.add_vertex(v);
  }
  for (int i = 0; i < adj_list.size(); i++) {
    for (int j : adj_list[i]) {
      if (i < j) {
        g.add_edge(i, j);
      }
    }
  }
  bliss::Stats stats;
  const unsigned int *cl = g.canonical_form(stats, report_aut, stdout);
  std::vector<int> perm(cl, cl+nn);
  bliss::Graph *cf = g.permute(cl);
  std::ostringstream sout;
  cf->write_dimacs(sout);
  delete cf;

  return {sout.str(), perm};
}



}  // namespace euler::pattern_mining
