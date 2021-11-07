#include <boost/multiprecision/cpp_dec_float.hpp>

#include "pattern_mining/gmine.h"
#include "util.h"
typedef boost::multiprecision::cpp_dec_float_100 value_type;

using namespace std;
using namespace euler;
using namespace euler::pattern_mining;



int main(int argc, char* argv[])
{
    
    graph::Graph_CSR_CPU g;

    if(argc==2)
    {
        g.read_graph(argv[1]);
    }
    else
    {
        printf("input error.\n");
        printf("./test.exe matrix.lg\n");
    }
    
    int p1 = calc_join2_sampleratio(g);
    
    printf("The sampling ratio of 2-size subgraphs is %d.\n",p1);
    
    
    
    auto pat2 = pattern_mining::PatListing::make_pattern(pattern_mining::PatListing().pattern_listing(2));
    // match pat2
    auto d2 = match(g, pat2, true, false, true);
    // join for pat3
    vector<SGList> sgls2 = { d2, d2 };
    auto H2 = build_tables(sgls2);
    Sampler *sm_tmp = new ProportionalSampler({ (double)p1, (double)p1 });
    auto [d3, ess3] = join<true, true, false, false, 2, 3, 3>(g, H2, sgls2, true, *sm_tmp, -1, false);
    vector<SGList> sgls3 = { d3 };
    auto H3 = build_tables(sgls3);
    int p2 = calc_join3_sampleratio(g,H3);

    
    printf("The sampling ratio of 3-size subgraphs is %d.\n",p2);
    

  return 0;
}
