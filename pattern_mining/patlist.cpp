  #include "patlist.h"
  
  namespace euler::pattern_mining {

  std::vector<std::vector<std::pair<int, int>>> PatListing::spanning_tree_listing(int n) {
    graph_t g(2);
    add_edge(0, 1, g);
    std::list<graph_t> tmp_pat_list, tmp_pat_list2;
    tmp_pat_list.push_back(g);

    for (int i = 2; i < n; i++) {
      tmp_pat_list2.clear();
      for (auto p1 = tmp_pat_list.begin(); p1 != tmp_pat_list.end(); ++p1) {
        for (int j = 0; j < i; j++) {
          graph_t gt(*p1);
          add_edge(j, i, gt);
          tmp_pat_list2.push_back(gt);
        }
      }
     // std::cout << tmp_pat_list2.size() << std::endl;
      for (auto p1 = tmp_pat_list2.begin(); p1 != tmp_pat_list2.end(); ++p1) {
        auto p2 = p1;
        if (p2 == tmp_pat_list2.end()) break;
        ++p2;
        while (p2 != tmp_pat_list2.end()) {
          std::vector<graph_t::vertex_descriptor> f(n);
          if (boost::isomorphism(*p1, *p2, boost::isomorphism_map(f.data()))) {
            tmp_pat_list2.erase(p2++);
          } else {
            p2++;
          }
        }
      }
      tmp_pat_list = tmp_pat_list2;
    }

    std::vector<std::vector<std::pair<int, int>>> res;
    for (auto &g : tmp_pat_list) {
      res.push_back(std::vector<std::pair<int, int>>());
      typename boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        res.back().push_back(
            std::make_pair(boost::source(*ei, g), boost::target(*ei, g)));
      }
    }
    return res;
  }
  }