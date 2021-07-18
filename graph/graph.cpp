#include "graph.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <functional>

namespace euler::graph {

  bool Graph_CSR_CPU::is_neighbor(int n1, int n2, const std::function<void(int, int)>& f) const {
    int low = rowptr[n1];
    int high = rowptr[n1 + 1] - 1;
    while (low <= high) {
      int mid = (low + high) / 2;
      if (colidx[mid] == n2) {
        if (f) f(n1, n2);
        return true;
      }
      else if (colidx[mid] < n2) {
        low = mid + 1;
      }
      else {
        high = mid - 1;
      }
    }
    return false;
  }

  std::pair<const int*, size_t> Graph_CSR_CPU::get_neighbors(int vid) const {
    return { &colidx[rowptr[vid]], rowptr[vid + 1] - rowptr[vid] };
  }

  void Graph_CSR_CPU::get_neighbors(int vid, std::vector<int>& neibs) const {
    neibs.insert(neibs.begin(), colidx.begin() + rowptr[vid], colidx.begin() + rowptr[vid + 1]);
  }

  int Graph_CSR_CPU::get_vertex_label(int vid) const {
    return vertex_labels[vid];
  }


  size_t Graph_CSR_CPU::max_degree() const {
    size_t l = 0;
    for (int i = 0; i < rowptr.size() - 1; i++) {
      size_t d = rowptr[i + 1] - rowptr[i];
      if (d > l) l = d;
    }
    return l;
  }

  std::vector<size_t> Graph_CSR_CPU::degree() const {
    std::vector<size_t> l;
    for (int i = 0; i < rowptr.size() - 1; i++) {
      size_t d = rowptr[i + 1] - rowptr[i];
      l.push_back(d);
    }
    return l;
  }

  void Graph_CSR_CPU::read_graph(const char* filename) {
    std::ifstream fin(filename);
    std::string line;
    while (std::getline(fin, line) && (line[0] == '#'))
      ;
    nnodes = 0;
    do {
      std::istringstream sin(line);
      std::string tmp;
      int v;
      int label;
      sin >> tmp >> v >> label;
      // cout << tmp << " " << v << " " << nnodes << endl;
      assert(tmp == "v" && v == nnodes);
      vertex_labels.push_back(label);
      nnodes++;
    } while (std::getline(fin, line) && (line[0] == 'v'));
    std::vector<std::vector<int>> adj_list(nnodes);
    do {
      std::istringstream sin(line);
      std::string tmp;
      int v1, v2, label;
      sin >> tmp >> v1 >> v2 >> label;
      assert(tmp == "e");
      // edge_set.insert(std::make_pair(v1, v2));
      // edge_set.insert(std::make_pair(v2, v1));
      adj_list[v1].push_back(v2);
      adj_list[v2].push_back(v1);
    } while (getline(fin, line));
    rowptr.push_back(0);
    for (int i = 0; i < nnodes; i++) {
      sort(adj_list[i].begin(), adj_list[i].end());
      int pos = 0;
      for (int j = 1; j < adj_list[i].size(); j++) {
        if (adj_list[i][j] != adj_list[i][pos]) adj_list[i][++pos] = adj_list[i][j];
      }

      if (adj_list[i].size() > 0)
        colidx.insert(colidx.end(), adj_list[i].data(), adj_list[i].data() + pos + 1);  // adj_list is sorted

      adj_list[i].clear();
      rowptr.push_back(colidx.size());
    }
    nedges = colidx.size() / 2;

    std::cout << "Graph read complete. Number of vertex: " << nnodes << std::endl;
  }

}  // namespace euler
