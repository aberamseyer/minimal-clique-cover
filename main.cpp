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

typedef adjacency_list<vecS, vecS, undirectedS> Graph;
typedef std::pair<int, int> Edge;
typedef graph_traits<Graph>::edge_iterator edge_iterator;
typedef graph_traits<Graph>::adjacency_iterator adjacency_iterator;

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
		ret.push_back(std::vector<std::vector<int>>());
		return ret;
	}
	std::vector<int> p = els[size - 1];
	els.pop_back();
	std::vector<std::vector<std::vector<int>>> subCombs = combinations(els);
	els.push_back(p);
	std::vector<std::vector<std::vector<int>>> copy = std::vector<std::vector<std::vector<int>>>(subCombs);
	for (int i = 0; i < subCombs.size(); ++i) {
		subCombs[i].push_back(p);
	}
	std::vector<std::vector<std::vector<int>>> un;
	std::set_union(copy.begin(), copy.end(), subCombs.begin(), subCombs.end(), std::inserter(un, un.end()));
//	std::cout << "els.size() = " << size << "\n";
//	std::cout << "copy.size() = " << copy.size() << "; subCombs.size() = " << subCombs.size() << "; size = " << size << "\n";
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
	while(input_file >> name1 >> name2) add_edge(name1, name2, g);

	std::vector<int> all_vertices;
	push_back(all_vertices, irange(0, (int)num_vertices(g)));
	bron_kerbosch(std::vector<int>(), all_vertices, std::vector<int>(), 0);
	std::cout << "We found " << maximal_cliques.size() << " maximal cliques.\n";
	for (const std::vector<int> c: maximal_cliques) {
		std::cout << "Clique: ";
		for (const int v: c) std::cout << v << " ";
		std::cout << "\n";
	}

	std::vector<std::vector<std::vector<int>>> combos = combinations(maximal_cliques);
//	for (auto a : testcombos) {
//		for (auto b : a) std::cout << b << " ";
//		std::cout << "\n";
//	}
//		std::cout << std::endl;
//	return 0;

	//long fewest = LONG_MAX;
	std::vector<std::vector<int>> *result = nullptr;
	std::cout << "combos.size() = " << combos.size() << "\n";
	for (std::vector<std::vector<int>> &f : combos) {
		// create a set of all the vertices covered by this set of maximal clique cover
		std::set<int> combination_verticies;
		for(std::vector<int> g : f) { // each graph in a single combination
			for (int h : g) { // each vetex in a single graph
				combination_verticies.insert(h);
			}
		}
		// check clique cover by computing the difference of verticies between those contained in this set of graphs and
		// all verticies in the entire graph (should be 0 if all vertices are in both)
		std::vector<int> diff;
		std::set_difference(all_vertices.begin(), all_vertices.end(), combination_verticies.begin(), combination_verticies.end(), std::inserter(diff, diff.begin()));
		// determine if this is a clique cover
		if (diff.size() != 0)
			continue;
		// && the size of this candidates element < fewest
		if (result == nullptr || f.size() < result->size()) {
			result = &f;
//			std::cout << "Found a new smallest: " << result->size() << "\n";
//			for (const std::vector<int> clique : f) {
//				std::cout << "Clique (in cover): ";
//				for (int v: clique)
//					std::cout << v << " ";
//				std::cout << "\n";
//			}
			std::cout << "Inside if, result->size() = " << result->size() << "\n";
			//fewest = result->size();
		}
		std::cout << "End of iteration, result->size() = " << result->size() << "\n";
	}
	std::cout << "After for, result->size() = " << result->size() << "\n";
	if (result != nullptr) {
		std::cout << "yay found one of size " << result->size() << std::endl;
//		for (const std::vector<int> clique : *result) {
//			std::cout << "Clique: ";
//			for (int v: clique)
//				std::cout << v << " ";
//			std::cout << "\n";
//		}
	}
	else
		std::cout << "bummer" << std::endl;

	//	print_edges();
    // write_graphviz(std::cout, g);
    std::cout << "Your clique cover is huge!" << std::endl;

    return 0;
}

void print_edges() {
    // print out edges using iterators
    std::pair<edge_iterator, edge_iterator> ei = edges(g);
    for(edge_iterator edge_iter = ei.first; edge_iter != ei.second; ++edge_iter) {
        std::cout << "(" << source(*edge_iter, g) << ", " << target(*edge_iter, g) << ")" << std::endl;
    }
}