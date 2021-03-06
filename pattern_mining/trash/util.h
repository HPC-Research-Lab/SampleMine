#pragma once

#include <vector>
#include <iostream>
#include <algorithm>


// ordered set of indices
struct Ordered_set {
	std::vector<int> items;

	Ordered_set(const std::vector<int> &t) {
		for (int i:t) items.push_back(i);
		std::sort(items.begin(), items.end());
	}

	Ordered_set() {}

	// intersection of two ordered sets
	Ordered_set operator * (Ordered_set const &b) { 
    Ordered_set res; 
		size_t i = 0, j = 0;
		while (i != items.size() && j != b.items.size()) {
			if (items[i] < b.items[j]) i++;
			else if (items[i] > b.items[j]) j++;
			else {
				res.items.push_back(items[i]);
				i++; j++;
			}
		}
		return res;
    } 

	// difference of two ordered sets
	Ordered_set operator - (Ordered_set const &b) { 
         Ordered_set res; 
		std::vector<int> common_idx; 
		size_t i = 0, j = 0;
		while (i != items.size() && j != b.items.size()) {
			if (items[i] < b.items[j]) i++;
			else if (items[i] > b.items[j]) j++;
			else {
				common_idx.push_back(i);
				i++; j++;
			}
		}
		size_t pos = 0;
		for (i=0; i<items.size(); i++) {
			if (pos >= common_idx.size()) res.items.push_back(items[i]);
			else if (i == common_idx[pos]) {pos++; continue;}
			else {res.items.push_back(items[i]);}
		}
		
		return res;
    } 

	int& operator [] (int idx) { 
		return items[idx];
    }

	size_t size() {
		return items.size();
	}

	void print() {
		for (int t: items) std::cout << t << " ";
		std::cout << std::endl;
	}
};

Ordered_set neighbor(int v, std::vector<int> &rowptr, std::vector<int> &colidx);

bool connected(int v1, int v2, std::vector<int> &rowptr, std::vector<int> &colidx);

size_t array_sort_and_hash(std::vector<int> &a); 

void print_vec(const std::vector<int> &a);

void permutation(std::vector<std::vector<int>> &all, std::vector<int> &a,
                   int l, int r);


