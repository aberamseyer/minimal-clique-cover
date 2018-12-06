#include <cstring>
#include <fstream>
#include "graph/Graph.h"
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <set>
#include <stdlib.h>
#include <string>
#include "ttmath/ttmath.h"
#include <vector>

#define TAG 0
#define MAX_NUM_MAXIMAL_CLIQUES 200
#define DEBUG true

std::vector<std::vector<int>> maximal_cliques;
Graph g;
typedef ttmath::UInt<TTMATH_BITS(MAX_NUM_MAXIMAL_CLIQUES)> bignum_t;
int num_processes;

std::vector<int> intersection(std::vector<int> &s1, std::vector<int> &s2) {
    std::vector<int> to_return;
    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(to_return));
    return to_return;
}

void sync_maximal_cliques (int my_rank){

	//convert to strings
	int num_cliques = maximal_cliques.size();
	std::string string_max_cliques[num_cliques];

	int i = 0;
	for(std::vector<int> clique : maximal_cliques){
		std::string string_clique = "";
		for(int v : clique){
			string_clique += std::to_string(v);
			if(v != clique.back())
				string_clique += ",";
		}
		string_max_cliques[i] = string_clique;
		i++;
	}
	std::string clique_library = "";
	for(std::string clique : string_max_cliques){
		clique_library += clique;
		if(clique != string_max_cliques[i - 1]) //last filled entry in array
				clique_library += "/";
	}

	if(my_rank != (num_processes-1) && maximal_cliques.size() > 0){
		clique_library += "/";
	}

	//std::cout << "My Rank: " <<my_rank << " my cliques: " << clique_library << "\n";


	//prep merge together
	int counts[num_processes];
	int displace[num_processes];
	int mySize = clique_library.size();
	MPI_Allgather(&mySize,1,MPI_INT,counts,1,MPI_INT,MPI_COMM_WORLD);

  	int displacement = 0;
	int vectors = 0;
	for(int i =0; i<num_processes ;i++){
		displace[i] = displacement;
		displacement += counts[i];
	}

	//merge
	char shared_clique_library[displacement]; //size of all characters

	MPI_Allgatherv(clique_library.c_str(),mySize,MPI_CHAR,&shared_clique_library,counts,displace,MPI_CHAR,MPI_COMM_WORLD);

	if(shared_clique_library[displacement-1] == '/'){
		shared_clique_library[displacement-1] = '\0';
	}
	/*
	if(my_rank == 0){
		std::cout <<"this right here";
		for(char c : shared_clique_library ){
			std::cout << c;
		}
	}*/

	//convert back
	maximal_cliques.clear();
	std::vector<char*> tempCliques;
	char *clique, *tofree, *temp_clique_library;
	tofree = temp_clique_library = strdup(shared_clique_library);

	while((clique = strsep(&temp_clique_library,"/")) != NULL){
		char *vertex;
		std::vector<int> tempClique;
		while((vertex = strsep(&clique,",")) != NULL){
			tempClique.push_back(std::atoi(vertex));
		}
		maximal_cliques.push_back(tempClique);
	}
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
		if(depth == 0) // only search the initialized node
			break;
	}
}

void bron_kerbosch (int my_rank, std::vector<int> all_vertices){

	int totalVertexNum = all_vertices.size();
	std::vector<int> my_vertices(all_vertices);
	std::vector<int> finished_vertices;

	for(int i=my_rank; i<totalVertexNum; i+=num_processes){

		int iterations;
		if(i == my_rank)
			iterations = my_rank;
		else
			iterations = num_processes;

		/*if (my_rank == 8) std::cout << "i=" << i << " j=" << j << " iterations=" << iterations << "\n";*/
		for(int j = 0; j<iterations;j++){
			finished_vertices.push_back(my_vertices[0]);
			my_vertices.erase(my_vertices.begin());
		}

		bron_kerbosch(std::vector<int>(), my_vertices, finished_vertices, 0);
	}
	/*
	for (const std::vector<int> c: maximal_cliques) {
		std::cout << "Clique: (rank " << my_rank << ") ";
		for (const int v: c) std::cout << v << " ";
		std::cout << std::endl;
	}*/

	if(num_processes > 1)
		sync_maximal_cliques(my_rank);
}

std::pair<long, long> my_range (unsigned long items, int threads, int my_rank) {
	double div = items / (double) threads;
	std::pair<long, long> to_return;
	//std::cout << "my_range(" << items << ", " << threads << ", " << my_rank << "); div = " << div << "\n";
	to_return.first = round(div * (double)my_rank);
	to_return.second = static_cast<long>(round(div * ++my_rank) - 1);
	return to_return;
}

// https://stackoverflow.com/a/9331125/863470
unsigned long long n_choose_k (unsigned n, unsigned k) {
	if (k > n) return 0;
	if (k * 2 > n) k = n - k;
	if (k == 0) return 1;
	long result = n;
	for (long i = 2; i <= k; ++i ) {
		result *= (n - i + 1);
		result /= i;
	}
	return result;
}

// https://math.stackexchange.com/a/1368570
std::vector<bool> nth_combination (long long n, long long r, int m) {
	std::vector<bool> out;
	out.reserve(n);
	for (unsigned i = 0; i < n; i++)
		out.push_back(false);

	int len = n;

	// In this loop, out.size() >= n >= 1
	for (unsigned long y; n > 0; n--) {
		auto n_minus_1 = n - 1;
		if (n > r && r >= 0) {
			y = n_choose_k(n_minus_1, r);
		} else {
			y = 0;
		}

		bool cond = m >= y;
		if (cond) {
			m -= y;
			out[len - n_minus_1 - 1] = true;
			r--;
		}
	}
	return out;
}

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
	if (my_rank == 0)
		std::cout << std::endl;

	double
		start_bron = 0,
		end_bron = 0,
		start_comb = 0,
		end_comb = 0;

	std::string current_exec_name = argv[0];
	std::string input_filename;

	int num_threads;
	if (argc != 3) {
		if (my_rank == 0)
			std::cout << "Usage: " << current_exec_name << " <input file> <num omp threads>" << std::endl;
		MPI_Finalize();
		return 0;
	}
	omp_set_num_threads(num_threads);

	input_filename = argv[1];
	std::ifstream input_file(input_filename);
	if (!input_file.good()) {
		if (my_rank == 0)
			std::cout << "File " << input_filename << " couldn't be opened for reading" << std::endl;
		MPI_Finalize();
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
	if (my_rank == 0 && !DEBUG) {
		std::cout << "built graph" << std::endl;
	}

	std::vector<int> all_vertices;
	all_vertices.reserve(g.get_num_vertices());
	for(unsigned i = 0; i < g.get_num_vertices(); ++i)
		all_vertices.push_back(i);

	if (my_rank == 0 && !DEBUG) {
		std::cout << "vertex count: " << (int) g.get_num_vertices() << ", finding maximal cliques.." << std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	start_bron = MPI_Wtime();
	bron_kerbosch(my_rank, all_vertices);
	MPI_Barrier(MPI_COMM_WORLD);
	end_bron = MPI_Wtime();
	
	if (my_rank == 0 && !DEBUG) {
		std::cout << "We found " << maximal_cliques.size() << " maximal cliques." << std::endl;
		for (const std::vector<int> c: maximal_cliques) {
			std::cout << "Clique: ";
			for (const int v: c) std::cout << v << " ";
			std::cout << std::endl;
		}
	}

	const unsigned long num_cliques = maximal_cliques.size(); // how many maximal cliques we found
	if (num_cliques > MAX_NUM_MAXIMAL_CLIQUES) {
		if (my_rank == 0)
			std::cout << MAX_NUM_MAXIMAL_CLIQUES << " is more than this program can handle, recompile it with a larger MAX_NUM_MAXIMAL_CLIQUES constant!" << std::endl;
		MPI_Finalize();
		return 0;
	}
	bignum_t global_combs = 0, local_combs = 0;
	std::vector<std::vector<int>*> result;
	int g_success_rank = -1, my_success_rank = -1;
	int g_done = 0, l_done = 0;

	MPI_Barrier(MPI_COMM_WORLD);
	start_comb = MPI_Wtime();

	// Loop through each num_cliques of clique cover, beginning at 1. num_included is the number of cliques that are to be in the cover
	for (unsigned cur_num_cliques = 1; cur_num_cliques <= num_cliques; ++cur_num_cliques) {

		unsigned long nCk = n_choose_k(static_cast<unsigned int>(num_cliques), cur_num_cliques);

		if (my_rank == 0 && !DEBUG)
			std::cout << (MPI_Wtime() - start_comb) << "s elapsed, checking " << nCk << " combinations with " << cur_num_cliques << " cliques.." << std::endl;

		// Loop through all clique covers for a given num_cliques
		// grab a range specific to my process rank
		std::pair<long, long> range = my_range(nCk, num_processes, my_rank);
		// starting mask will be different for each process
		auto orig_mask = nth_combination(num_cliques, cur_num_cliques, range.second);
		int num_iters = ceil((double)nCk / num_processes);

		volatile bool stop = false;
		#pragma omp parallel for shared(stop)
		for (int i = 0; i < num_iters; i++) {
			if (stop)
				i = num_iters; // break loop

			long permutation_num = range.second - i;

			// Does this thread have any more combinations to check?
			if (permutation_num >= range.first) {
				std::vector<bool> mask;
				#pragma omp critical
				{
					std::next_permutation(orig_mask.begin(), orig_mask.end());
					mask = std::vector<bool>(orig_mask);
				}
				local_combs++;
				// the set of cliques that is possibly a (minimal) clique cover
				std::set<int> ints_in_clique_set;
				// what might be our minimal clique cover
				std::vector<std::vector<int> *> candidate;
				for (unsigned k = 0; k < num_cliques; ++k)
					if (mask[k]) {
						// all the vertices that the candidate has within it
						for (const int a : maximal_cliques[k])
							ints_in_clique_set.insert(a);
						candidate.push_back(&maximal_cliques[k]);
					}

				// determine if this is a clique cover
				if (result.size() == 0 && ints_in_clique_set.size() == all_vertices.size()) {
					result = candidate;
					my_success_rank = my_rank;
					l_done = 1;
					stop = true; // break loop
				}
			}

	    } // end loop through all potential covers within a single num_cliques
		// has any thread found a cover?
		MPI_Allreduce(&l_done, &g_done, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
		if (g_done)
			break;
	} // end loop through all potential cover sizes

	MPI_Barrier(MPI_COMM_WORLD);
	end_comb = MPI_Wtime();
	// add up how much work we did
	MPI_Allreduce(&local_combs, &global_combs, TTMATH_BITS(MAX_NUM_MAXIMAL_CLIQUES), MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	// figure out which process found a solution
	MPI_Allreduce(&my_success_rank, &g_success_rank, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	if (my_rank == g_success_rank) {
		std::cout << std::endl;
		if (!result.empty()) {
			if (!DEBUG) {
				std::cout << "yay found a cover with " << result.size() << " cliques in it (inside thread " << my_rank << ")" << std::endl;
				std::cout << "Minimum Clique Cover: " << std::endl;
				for (const std::vector<int>* clique : result) {
					for (const int &v: *clique)
						std::cout << v << " ";
					std::cout << std::endl;
				}
			}
			bignum_t total_combs = 2; total_combs.Pow(num_cliques);
			std::cout << "Bron Kerbosch run time: " << (end_bron - start_bron) << "s" << std::endl;
			std::cout << "Combination run time: " << (end_comb - start_comb) << "s" << std::endl;
			std::cout << "Total run time: " << (end_bron + end_comb - start_bron - start_comb) << "s" << std::endl;
			std::cout << "Total number of combinations to check: " << total_combs << std::endl;
			std::cout << "We checked " << global_combs << " before finding the solution, skipping " << total_combs-global_combs << std::endl;
		}
		else {
			std::cout << "something went wrong, couldn't find a clique cover" << std::endl;
		}
		if (!DEBUG)
			std::cout << std::endl << "Your clique cover is huge!" << std::endl;
	}
	
	if (my_rank == 0)
		std::cout << std::endl;
	MPI_Finalize();
	return 0;
}
