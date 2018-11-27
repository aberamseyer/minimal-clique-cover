#include <fstream>
#include "Graph.h"
#include <iostream>
#include <math.h>
#include <omp.h>
#include <set>
#include <stdlib.h>
#include <string>
#include <vector>

std::vector<std::vector<int>> maximal_cliques;
Graph g;

std::vector<int> intersection(std::vector<int> &s1, std::vector<int> &s2) {
    std::vector<int> to_return;
    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(to_return));
    return to_return;
}
void bron_kerbosch (const std::vector<int> &r, std::vector<int> p, std::vector<int> x, int depth) {
	if (p.empty()) {
		if (x.empty()) maximal_cliques.push_back(r);
		return;
	}
	while (!p.empty()) {
		int v = p[0];
		std::vector<int> this_r = r;
		std::vector<int> this_p = p;
		std::vector<int> this_x = x;
		std::vector<int> neighborhood;
		auto neighbors = g.get_adjacent_vertices(v);
        neighborhood.reserve(neighborhood.size());
        for (auto vd : neighbors)
		    neighborhood.push_back(vd);

		this_r.push_back(v);
		// bron_kerbosch(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
		bron_kerbosch(this_r, intersection(this_p, neighborhood), intersection(this_x, neighborhood), depth + 1);
		// Move v from P to X
		p.erase(p.begin());
		x.push_back(v);
	}
}

std::vector<bool> nth_permutation (unsigned long size, unsigned num_trues, int index) {

	std::vector<bool> mask;
	mask.reserve(size);
	for (unsigned i = 0; i < size; i++) {
		bool val = i < num_trues;
		mask.push_back(val);
	}

	std::vector<bool> result;
	result.reserve(size);
    for (unsigned i = 0; i < size; i++) {
        unsigned long item = index % size;
        index /= size;
        result.push_back(mask[item]);
        mask.erase(mask.begin() + item);
    }

    return result;
}

unsigned n_choose_k (unsigned n, unsigned k) {
	if (k > n) return 0;
	if (k * 2 > n) k = n-k;
	if (k == 0) return 1;

	unsigned result = n;
	for (unsigned i = 2; i <= k; ++i )
		result *= (n-i+1) / i;
	return result;
}

int main(int argc, char* argv[]) {

	omp_set_num_threads(24);

	std::string current_exec_name = argv[0];
	std::string input_filename;
	
	if (argc != 2) {
		std::cout << "Usage: " << current_exec_name << " <input file>" << std::endl;
		return 0;
	}
	input_filename = argv[1];
	std::ifstream input_file(input_filename);
	if (!input_file.good()) {
		std::cout << "File " << input_filename << " couldn't be opened for reading" << std::endl;
		return 0;
	}
	int name1, name2;
	// format of file always has 2 nodes on every line
	std::vector<std::pair<int, int>> data;
	while(input_file >> name1 >> name2) {
	    std::pair<int, int> to_add (name1, name2);
	    data.push_back(to_add);
    }
    g.set_size(static_cast<unsigned long>(name1));
	for(auto i : data) {
		g.add_edge(i.first, i.second);
	}
    std::cout << "built graph" << std::endl;

	std::vector<int> all_vertices;
    all_vertices.reserve(g.get_num_vertices());
    for(unsigned i = 0; i < g.get_num_vertices(); ++i)
	    all_vertices.push_back(i);
	std::cout << "vertex count: " << g.get_num_vertices() << ", finding maximal cliques.." << std::endl;
	bron_kerbosch(std::vector<int>(), all_vertices, std::vector<int>(), 0);
	std::cout << "We found " << maximal_cliques.size() << " maximal cliques." << std::endl;
//	for (const std::vector<int> c: maximal_cliques) {
//		std::cout << "Clique: ";
//		for (const int v: c) std::cout << v << " ";
//		std::cout << std::endl;
//	}

	unsigned long size = maximal_cliques.size(); // how many maximal cliques we found
	std::vector<std::vector<int>*> result;
	bool done = false;

	// Loop through each size of clique cover, beginning at 1. i is the number of cliques that are to be in the cover
	for(unsigned num_included = 1; num_included <= size && !done; ++num_included) {
		if (num_included == 1) {
			std::cout << "checking for a complete Graph.." << std::endl;
		} else {
			std::cout << "checking combinations with " << num_included << " cliques.." << std::endl;
		}

		// Loop through all clique covers for a given size
		int nCk = n_choose_k(static_cast<unsigned int>(size), num_included);
#pragma omp parallel for
		for (int j = 0; j < nCk; j++) {
			// the set of cliques that is possibly a (minimal) clique cover
			std::set<int> ints_in_clique_set;
			// what might be our minimal clique cover
			std::vector<std::vector<int>*> candidate;
			std::vector<bool> mask = nth_permutation(size, num_included, j);
#pragma omp parallel for
			for (unsigned k = 0; k < size; ++k) {
				if (mask[k]) {
					// all the vertices that the candidate has within it
					for(const int a : maximal_cliques[k]) {
						ints_in_clique_set.insert(a);
					}
					candidate.push_back(&maximal_cliques[k]);
				}
			}
	        // determine if this is a clique cover
	        if (ints_in_clique_set.size() != all_vertices.size()) {
	            continue;
	        }
	        // && the size of this candidates element < fewest
	        if (result.empty() || candidate.size() < result.size()) {
	            result = candidate;
	            done = true;
//#pragma omp single
//	        	break;
	        }
		}
	}

	std::cout << std::endl;
	if (!result.empty()) {
		std::cout << "yay found a cover one with " << result.size() << " cliques in it" << std::endl;
		std::cout << "Minimum Clique Cover: " << std::endl;
		for (const std::vector<int>* clique : result) {
			for (const int &v: *clique)
				std::cout << v << " ";
			std::cout << std::endl;
		}
	}
	else {
		std::cout << "something went wrong, couldn't find a clique cover" << std::endl;
	}

    std::cout << std::endl << "Your clique cover is huge!" << std::endl;
    return 0;
}
