#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/range/irange.hpp>

using namespace boost;

typedef adjacency_list<setS, vecS, undirectedS> Graph;
//typedef graph_traits<Graph>::edge_iterator edge_iterator;

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

std::vector<std::vector<std::vector<int>>> combinations(std::vector<std::vector<int>> els) {
	unsigned long size = els.size();
	if (size == 0) {
		std::vector<std::vector<std::vector<int>>> ret;
		ret.emplace_back();
		return ret;
	}
	std::vector<int> p = els[size - 1];
	els.pop_back();
	std::vector<std::vector<std::vector<int>>> subCombs = combinations(els);
	els.push_back(p);
	std::vector<std::vector<std::vector<int>>> copy = std::vector<std::vector<std::vector<int>>>(subCombs);
	for (auto &subComb : subCombs) {
        subComb.push_back(p);
	}
	std::vector<std::vector<std::vector<int>>> un;
	std::set_union(copy.begin(), copy.end(), subCombs.begin(), subCombs.end(), std::inserter(un, un.end()));
	return un;
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
	int name1, name2;
	// format of file always has 2 nodes on every line
	while(input_file >> name1 >> name2) {
		add_edge(name1, name2, g);
	}

	std::vector<int> all_vertices;
	push_back(all_vertices, irange(0, (int)num_vertices(g)));
	std::cout << "vertex count: " << (int)num_vertices(g) << ", finding maximal cliques.." << std::endl;
	bron_kerbosch(std::vector<int>(), all_vertices, std::vector<int>(), 0);
	std::cout << "We found " << maximal_cliques.size() << " maximal cliques.\n";
//	for (const std::vector<int> c: maximal_cliques) {
//		std::cout << "Clique: ";
//		for (const int v: c) std::cout << v << " ";
//		std::cout << std::endl;
//	}

	std::vector<std::vector<std::vector<int>>> combos = combinations(maximal_cliques);
	combos.erase(combos.begin());
//	for(auto a : combos) {
//		std::cout << "set: " << std::endl;
//		for (auto b : a) {
//			for (auto c : b) {
//				std::cout << c;
//			}
//			std::cout << std::endl;
//		}
//		std::cout << std::endl;
//	}
	std::vector<std::vector<int>> *result = nullptr;
	std::cout << "Number of combinations to check: " << combos.size() << std::endl;
	for (std::vector<std::vector<int>> &f : combos) {
		// create a set of all the vertices covered by this set of maximal clique cover
		std::set<int> combination_verticies;
		for(std::vector<int> g : f) { // each graph in a single combination
			for (int h : g) { // each vertex in a single graph
				combination_verticies.insert(h);
			}
		}
		// check clique cover by computing the difference of verticies between those contained in this set of graphs and
		// all verticies in the entire graph (should be 0 if all vertices are in both)
		std::vector<int> diff;
		std::set_difference(all_vertices.begin(), all_vertices.end(), combination_verticies.begin(), combination_verticies.end(), std::inserter(diff, diff.begin()));
		// determine if this is a clique cover
		if (!diff.empty())
			continue;
		// && the size of this candidates element < fewest
		if (result == nullptr || f.size() < result->size()) {
			result = &f;
		}
	}
	if (result != nullptr) {
		std::cout << "yay found one of size " << result->size() << std::endl;
		std::cout << "Minimum Clique Cover: " << std::endl;
		for (const std::vector<int>& clique : *result) {
			for (const int &v: clique)
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
