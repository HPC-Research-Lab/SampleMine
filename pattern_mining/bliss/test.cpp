#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "defs.hh"
#include "graph.hh"
#include "timer.hh"
#include "utils.hh"

int main() {
  bliss::Graph g;

  bliss::Stats stats;

  g.add_vertex(0);
  g.add_vertex(0);
  g.add_vertex(0);

  g.add_edge(0, 1);
  g.add_edge(1, 2);

  const unsigned int* cl = g.canonical_form(stats, NULL, stdout);
  bliss::Graph *cf = g.permute(cl);
  cf->write_dimacs(stdout);
  return 0;
}