#ifndef MINIMAL_CLIQUE_COVER_GRAPH_H
#define MINIMAL_CLIQUE_COVER_GRAPH_H

#include <algorithm>
#include <set>
#include <vector>

class Graph {
private:
    unsigned long num_vertices;
    std::vector<std::vector<bool>> matrix;
    Graph* setMatrix (int v1, int v2, bool connected);
public:
//    Graph(unsigned long size);
//    ~Graph();
    Graph* add_edge(int v1, int v2);
//    Graph* add_vertex(int vertex);
    bool are_connected(int v1, int v2);
    std::set<int> get_adjacent_vertices(int vertex);
    std::set<std::pair<int, int>> get_edges();
    unsigned long get_num_vertices();
    std::set<int> get_vertices();
    Graph* remove_vertex(int vertex);
    Graph* set_size(unsigned long size);
};


#endif //MINIMAL_CLIQUE_COVER_GRAPH_H
