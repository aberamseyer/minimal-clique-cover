#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/range/irange.hpp>
#include <fstream>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <set>
#include <stdlib.h>
#include <string>
#include <vector>

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
	while(input_file >> name1 >> name2)
		add_edge(name1, name2, g);

	unsigned int num_maximal_cliques;

	MPI_Init(&argc, &argv);
	int my_rank, num_processes;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    std::vector<int> all_vertices;
    if (my_rank == 0) {
		push_back(all_vertices, irange(0, (int) num_vertices(g)));
		std::cout << "vertex count: " << (int) num_vertices(g) << ", finding maximal cliques.." << std::endl;
		bron_kerbosch(std::vector<int>(), all_vertices, std::vector<int>(), 0);
		std::cout << "We found " << maximal_cliques.size() << " maximal cliques." << std::endl;
		//	for (const std::vector<int> c: maximal_cliques) {
		//		std::cout << "Clique: ";
		//		for (const int v: c) std::cout << v << " ";
		//		std::cout << std::endl;
		//	}

		num_maximal_cliques = maximal_cliques.size(); // how many maximal cliques we found
	}
    MPI_Bcast(&num_maximal_cliques, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	std::vector<std::vector<int>*> result;

	bool done = false;
	// We need to go from 0...maximal_cliques - 1
	// That's num_maximal_cliques iterations
	for(unsigned int i = 1; i <= num_maximal_cliques && !done; ++i) {
		if (i == 1) {
			std::cout << "checking for a complete graph.." << std::endl;
		} else {
			std::cout << "checking combinations with " << i << " cliques.." << std::endl;
		}
		std::vector<bool> mask(num_maximal_cliques);
		std::fill_n(mask.begin(), i, true);
		do {
			// the set of cliques that is possibly a (minimal) clique cover
			std::set<int> ints_in_clique_set;
			// what might be our minimal clique cover
			std::vector<std::vector<int>*> candidate;
			for (unsigned int j = 0; j < num_maximal_cliques; ++j) {
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
		} while (std::prev_permutation(mask.begin(), mask.end()) && !done);
	}

	std::cout << std::endl;
	if (!result.empty()) {
		std::cout << "yay found a cover one with " << result.size() << " cliques in it" << std::endl;
		std::cout << "Mininum Clique Cover: " << std::endl;
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
