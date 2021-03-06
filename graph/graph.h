#pragma once

#include <functional>
#include <vector>

namespace euler::graph {

class Graph {
 protected:
  int nnodes, nedges;

 public:
  virtual void read_graph(const char *filename) = 0;
  virtual bool is_neighbor(
      int x, int y, const std::function<void(int, int)> &f = nullptr) const = 0;
  virtual int get_vertex_label(int vid) const { return -1; };
  virtual int get_edge_label(int eid) const { return -1; };
  virtual void get_neighbors(int vid, std::vector<int> &neibs) const = 0;
  virtual  std::pair<const int*, size_t> get_neighbors(int vid) const = 0;

  virtual size_t max_degree() const = 0;
  virtual  std::vector<size_t> degree() const = 0;

  size_t num_nodes() const { return nnodes; }
  size_t num_edges() const { return nedges; }
};

class Graph_CSR_CPU : public Graph {
  std::vector<int> rowptr;
  std::vector<int> colidx;
  std::vector<int> vertex_labels;
  std::vector<int> edge_labels;

 public:
  void read_graph(const char *filename) override;
  bool is_neighbor(
      int x, int y,
      const std::function<void(int, int)> &f = nullptr) const override;
  int get_vertex_label(int vid) const override;
  void get_neighbors(int vid, std::vector<int> &neibs) const override;
  std::pair<const int*, size_t> get_neighbors(int vid) const override;
  size_t max_degree() const override;
  std::vector<size_t> degree() const override;
};
}  // namespace euler::graph