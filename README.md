# SampleMine

## Overview

SampleMine is a general-purpose system for subgraph pattern mining based on subgraph enumeration and sampling.

It is built around two novel ideas: two-vertex exploration and multi-stage subgraph sampling. SampleMine achieves significant speedups against the state-of-the-art graph mining systems and supports larger pattern mining tasks that are difficult for the existing systems.

Our system supports conventional SPM tasks as well as arbitrary user-defined tasks.
It Supports the following tasks:

* **Subgraph Counting**: Counting the embeddings of different subgraph patterns and find the patterns with the largest counts.Subgraph counting considers labeled patterns
* **Frequent Subgraph Mining**: Obtaining all frequent subgraph patterns from a labeled input graph.
* **User-Defined Queries**: Processing any user-defined SPM tasks.



## Building SampleMine

### Installing Dependencies

Download Boost Library version 1.76
```Shell
wget https://boostorg.jfrog.io/native/main/release/1.76.0/source/boost_1_76_0.tar.gz
tar -xzvf boost_1_76_0.tar.gz
```
Install GCC version 7
```Shell
#For CentOS Users
#Link: https://www.softwarecollections.org/en/scls/rhscl/devtoolset-7/
sudo yum install centos-release-scl
sudo yum-config-manager --enable rhel-server-rhscl-7-rpms
sudo yum install devtoolset-7
source /opt/rh/devtoolset-7/enable

#For Ubuntu Users
sudo apt-get update
sudo apt install g++-7
```

### Compiling SampleMine

Download Source Code
```Shell
git clone SampleMine
```
Add Enviroment Variable BOOST_ROOT for compiling
```Shell
export BOOST_ROOT=PATH_TO_BOOST/boost_1_76_0/
```
Enter into the directory and Compiling SampleMine
```Shell
cd GPM/pattern_mining/
make -j 16
```



## Testing SampleMine

We provide some program examples for testing. You can use our programs to test directly, or you can use the relevant API to rewrite the program in the test folder. 

Before running the program, you need to export two environment variables
```Shell
#The path of database files that stores the intermediate data
export DB_PATH=./db
#The number of parallel OPENMP threads 
export OMP_NUM_THREADS=16
```

### Testing using program examples

Please modifying the EXE variable in Makefile to compile the corresponding program

**Subgraph Counting:**  
We provide a test example of 5-size Subgraph Counting
```Shell
#example: ./sc5_tv.exe data/citeseer 2
./sc5_tv.exe GraphPath sampling_ratio
```

**Frequent Subgraph Mining:**  
We provide a test example of 5-size Frequent Subgraph Mining:
```Shell
.#example: ./fsm5_two_vertex.exe data/citeseer 0.001 4
./fsm5_two_vertex.exe GraphPath support_threshold sampling_ratio
```

**User-Defined Queries:**
```Shell

```


### Writing your own program to perform the tasks

You can refer to the program examples provided by us to make modifications
We provide the following commonly used APIs
```cpp
//Define a graph and read in
graph::Graph_CSR_CPU g;
g.read_graph(onst char *filename);

//List all patterns of size n
std::vector<std::vector<std::pair<int, int>>> pattern_listing(int n)

//Encapsulating PatListing vector into pattern format
std::vector<std::shared_ptr<Pattern>> 
pattern_mining::PatListing::make_pattern(
	const std::vector<std::vector<std::pair<int, int>>> &pat)

//Using match function to get size=3 or size=2 embedding in the graph
//If the target task is fsm, set edge_induced to true, if it is subgraph counting, set edge_induced to false
SGList match(
    const graph::Graph &g, 
    const PatList &patterns,
    bool store_data = true, 
    bool edge_induced=false, 
    bool output_labeled=false,  
    double mni=-1, 
    bool testing = false, 
    bool pattern_labeled = false, 
    double sampling_threshold=0);

//Defining sampling ratio for each join element
class ProportionalSampler  {
    ProportionalSampler(const std::vector<double>& sample_ratio) 
};

//Specify the size of the vector according to the size of the target size. 
//For example: 
//If you need to get the embedding of size 5, you need to set it to {d3, d3} because 3+3-1=5. 
//If you need to get the embedding of size 6, you need to set it to {d3, d3, d2} because 3+3+2-2=6.
//For each target Size, the corresponding array can be derived 
vector<SGList> sgls = { d3, d3 };

/*
The essential function of the SPM task, joining the small embeddings obtained by the match function to obtain the larger embeddings.

If the target task is fsm, set edge_induced to true, if it is subgraph counting, set edge_induced to false

According to the target pattern size, you need to change the value of ncols1, ncols2...ncols
If you are joining two 3 embedding lists, you need to write as join<true, true, false, false, 2, 4, 4>
If you are joining two size-3 embedding lists and one size-2 embedding lists, you need to write as join<true, true, false, false, 3, 4, 4, 3>
*/
template <bool pat_agg, 
bool has_labels, bool edge_induced, bool mni, int K, size_t ncols1, size_t ncols2, size_t... ncols>
  std::tuple<SGList, std::map<std::string, double>> join(const graph::Graph& g, const std::vector<std::vector<std::shared_ptr<db::MyKV<int>>>>& H, const std::vector<SGList>& sgls, bool store, Sampler& sampler = default_sampler, double mni_threshold = -1, bool need_actual_pattern = false, bool est = false, bool adaptive_sampling = false, const std::pair<std::unordered_set<unsigned long>, std::unordered_set<unsigned long>>& sgl3 = join_dummy1, Query& query = default_query) {
  
```


**Size-5 Subgraph Counting example:**
```cpp
int main(int argc, char* argv[]) {
  graph::Graph_CSR_CPU g;
  g.read_graph(argv[1]);

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));

  double st2 = atof(argv[2]);

  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true);

  cout << "start join for pat3: " << endl;
  vector<SGList> sgls2 = { d2, d2 };

  util::Timer match_time;
  match_time.start();
  auto H2 = build_tables(sgls2);
  auto [d3, ess3] = join<true, true, false, false, 2, 3, 3>(g, H2, sgls2, true, default_sampler, -1, true);
  match_time.stop();

  cout << "join for pat3 time: " << match_time.get() << " sec" << endl;
  size_t npat3 = d3.sgl->size();
  cout << "num of size-3 patterns: " << npat3 << endl;

  vector<SGList> sgls = { d3, d3 };

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  cout << "build table done" << endl;

  Sampler *sm2;
  if (st2 > 0)
    sm2 = new ProportionalSampler({ st2, st2 });
  else sm2 = &default_sampler;

  util::Timer t;
  t.start();
  auto [d_res, ess] = join<true, true, false, false, 2, 4, 4>(g, H, sgls, false, *sm2, -1, false, st2 > 0);
  t.stop();

  cout << "Time: " << t.get() << " sec, ";
  if (d_res.sgl) {
    if (st2 == 0) {
      cout << "Num patterns: " << d_res.sgl->keys.size() << endl;

      vector<double> counts(d_res.sgl->count.begin(), d_res.sgl->count.end());
      sort(counts.begin(), counts.end(), std::greater<double>());

      for (int i = 0; i < (50 > counts.size() ? counts.size() : 50); i++) {
        cout << counts[i] << endl;
      }
    }
    else {
      cout << "Num patterns: " << ess.size() << endl;

      vector<double> counts;

      for (auto& [k, v] : ess) counts.push_back(v);

      sort(counts.begin(), counts.end(), std::greater<double>());

      for (int i = 0; i < (50 > counts.size() ? counts.size() : 50); i++) {
        cout << counts[i] << endl;
      }
    }
  }
  else {
    cout << "Num patterns: 0" << endl;
  }

  return 0;
}

```
**Size-5 Frequent Subgraph Mining example:**
```cpp
int main(int argc, char* argv[]) {
  graph::Graph_CSR_CPU g;
  g.read_graph(argv[1]);

  auto pat2 = pattern_mining::PatListing::make_pattern(
    pattern_mining::PatListing().pattern_listing(2));

  double thh = atof(argv[2]);

  double st2 = atof(argv[3]);

  double sup = (size_t)round(thh * g.num_nodes());

  cout << "support threshold: " << sup << endl;


  cout << "start matchings pat2: " << endl;
  auto d2 = match(g, pat2, true, false, true, sup);

  filter(d2, sup);
  cout << "num of size-2 frequent patterns: " << d2.sgl->size() << endl;

  double scaled_st2 = scale_sampling_param(d2, st2);

  cout << "start join for pat3: " << endl;
  vector<SGList> sgls2 = { d2, d2 };

  util::Timer match_time;
  match_time.start();
  auto H2 = build_tables(sgls2);

  auto [d3, ess3] = join<true, true, true, true, 2, 3, 3>(g, H2, sgls2, true, default_sampler, sup, true);

  match_time.stop();

  cout << "join for pat3 time: " << match_time.get() << " sec" << endl;

  filter(d3, sup);

  size_t npat3 = d3.sgl->size();

  cout << "num of size-3 frequent patterns: " << npat3 << endl;

  double st2_scaled = scale_sampling_param(d2, st2);

  cout << "scaled sampling param: " << st2_scaled << endl;

  vector<SGList> sgls = { d3, d3 };

  cout << "building tables..." << endl;
  auto H = build_tables(sgls);
  auto  subgraph_hist = get_subgraph_hist(sgls);
  cout << "build table done" << endl;

  Sampler *sm2;
  if (st2 > 0)
    sm2 = new BudgetSampler(subgraph_hist, {st2_scaled, st2_scaled});
  else sm2 = &default_sampler;

  util::Timer t;
  t.start();
  auto [d_res, ess] = join<true, true, true, true, 2, 4, 4>(g, H, sgls, false, *sm2, sup, false, false);
  t.stop();

  filter(d_res, sup);

  cout << "Time: " << t.get() << " sec, ";
  if (d_res.sgl) {
    cout << "Num patterns: " << d_res.sgl->keys.size() << endl;
  }
  else {
    cout << "Num patterns: 0" << endl;
  }

  return 0;
}
```

**User-Defined Queries:**
```cpp


```

