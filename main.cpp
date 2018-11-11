#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/range/irange.hpp>

using namespace boost;

typedef adjacency_list<setS, vecS, undirectedS> Graph;

std::vector<std::vector<int>> maximal_cliques;
Graph g;

std::vector<int> intersection (std::vector<int> v1, std::vector<int> v2) {
	std::vector<int> intersection;
	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(intersection));
	return intersection;
}

void bron_kerbosch (std::vector<int> r, std::vector<int> p, std::vector<int> x, int depth) {
	if (p.size() == 0) {
		if (x.size() == 0) maximal_cliques.push_back(r);
		return;
	}
	while (p.size() > 0) {
		int v = p[0];
		std::vector<int> this_r = r;
		std::vector<int> this_p = p;
		std::vector<int> this_x = x;
		std::vector<int> neighborhood;
		auto neighbors = adjacent_vertices(v, g);
		for (auto vd : make_iterator_range(neighbors)) neighborhood.push_back(vd);

		this_r.push_back(v);
		// bron_kerbosch(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
		bron_kerbosch(this_r, intersection(this_p, neighborhood), intersection(this_x, neighborhood), depth + 1);
		// Move v from P to X
		p.erase(p.begin());
		x.push_back(v);
	}
}

int main(int argc, char* argv[]) {
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
	while(input_file >> name1 >> name2) {
		add_edge(name1, name2, g);
	}

	std::vector<int> all_vertices;
	push_back(all_vertices, irange(0, (int)num_vertices(g)));
	std::cout << "vertex count: " << (int)num_vertices(g) << ", finding maximal cliques.." << std::endl;
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
	for(int i = 1; i <= size && !done; ++i) {
		if (i == 1) {
			std::cout << "checking for a complete graph.." << std::endl;
		} else {
			std::cout << "checking combinations with " << i << " cliques.." << std::endl;
		}
		std::vector<bool> mask(size);
		std::fill_n(mask.begin(), i, true);
		do {
			// the set of cliques that is possibly a (minimal) clique cover
			std::set<int> ints_in_clique_set;
			// what might be our minimal clique cover
			std::vector<std::vector<int>*> candidate;
			for (int j = 0; j < size; ++j) {
				if (mask[j]) {
					// all the vertices that the candidate has within it
					for(const int a : maximal_cliques[j]) {
						ints_in_clique_set.insert(a);
					}
					candidate.push_back(&maximal_cliques[j]);
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
	        	break;
	        }			
		} while (std::prev_permutation(mask.begin(), mask.end()));
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
