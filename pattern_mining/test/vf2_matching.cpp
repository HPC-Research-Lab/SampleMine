#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
using namespace boost;
using namespace std;


typedef adjacency_list< vecS, vecS, undirectedS > graph_type;

void print_vec(vector<int> &a) {
	for (int i: a) cout << i << " ";
	cout << endl;
}

size_t array_sort_and_hash(vector<int> &a) {
    sort(a.begin(), a.end());
	size_t result = 0;
	size_t const h = a[0];
    for (int p: a) {
        result ^= h + 0x9e3779b9 + (result << 6) + (result >> 2);
    }
    return result;
}

struct my_callback {
        static size_t count;

    my_callback(const graph_type& graph1, const graph_type& graph2)
    : graph1_(graph1), graph2_(graph2) {}

    template <typename CorrespondenceMap1To2,
            typename CorrespondenceMap2To1>
    bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) const {
        graph_traits<graph_type>::vertex_iterator vi, vi_end;
       // vector<int> m;
       // for (tie(vi, vi_end) = vertices(graph1_); vi != vi_end; ++vi) {
         //   m.push_back(get(f, *vi));
       // }
       // vector<int> d(m.begin(), m.end());
       // size_t t = array_sort_and_hash(d);
      //  if (hist.find(t) == hist.end()) hist[t] = map<vector<int>, vector<int>>();
      //  if (hist[t].find(d) == hist[t].end()) {
     //       hist[t][d] = m;
     //       count++;
            //print_vec(m);
     //   }
      //  if (count >= 1024l * 1024l * 1024l) hist.clear();
      //  #ifdef _COUNT
      //  if (count >= _COUNT) return false;
      //  #endif
      count++;
        return true;
    }
  private:
    const graph_type& graph1_;
    const graph_type& graph2_;
    static map<size_t, map<vector<int>, vector<int>>> hist;
 };

size_t my_callback::count = 0;
map<size_t, map<vector<int>, vector<int>>> my_callback::hist;

graph_type read_graph(char *filename) 
{
    ifstream fin(filename);	
	string line;
	while(getline(fin, line) && (line[0] == '%'));
	istringstream sin(line);
    int nnodes, nedges;
	sin >> nnodes >> nnodes >> nedges;
    graph_type g(nnodes);
	int v1, v2;
	string label;
	int t = 0;
	while(getline(fin, line)) {
		istringstream sin1(line);
		sin1 >> v1 >> v2 >> label;
		v1--; v2--;
        add_edge(v1, v2, g);
		t++;
	}
    assert(t == nedges);
    return g;
}

int main(int argc, char *argv[])
{
    graph_type graph1 = read_graph(argv[1]);

    // Build graph2
    graph_type graph2 = read_graph(argv[2]);


    // Create callback to print mappings
    my_callback callback(graph1, graph2);

    // Print out all subgraph isomorphism mappings between graph1 and graph2.
    // Vertices and edges are assumed to be always equivalent.
    vf2_subgraph_iso(graph1, graph2, callback);

    cout << my_callback::count << endl;


    return 0;
}

