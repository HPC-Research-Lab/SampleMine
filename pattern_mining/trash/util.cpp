#include "util.h"
#include <algorithm>
using namespace std;

Ordered_set neighbor(int v, vector<int> &rowptr, vector<int> &colidx) {
	Ordered_set res;
	res.items = vector<int>(&colidx[rowptr[v]], &colidx[rowptr[v+1]]);
	//sort(res.items.begin(), res.items.end());
	return res;
}

bool connected(int v1, int v2, vector<int> &rowptr, vector<int> &colidx) {
	for (int i=rowptr[v1]; i<rowptr[v1+1]; i++) {
		if (v2 == colidx[i]) return true;
	}
	return false;
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

void print_vec(const vector<int> &a) {
	for (int i: a) cout << i << " ";
	cout << endl;
}