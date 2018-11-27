#include "Graph.h"

// Public Methods

bool Graph::are_connected(int v1, int v2) {
    return matrix[v1][v2];
}

Graph* Graph::add_edge(int v1, int v2) {
    return setMatrix(v1, v2, true);
}

//Graph* Graph::add_vertex(int vertex) {
//    num_vertices++;
//    for (std::vector<bool> row : matrix)
//        row.insert(row.begin() + vertex, false);
//    std::vector<bool> to_add;
//    std::fill_n(to_add.begin(), num_vertices, false);
//    matrix.insert(matrix.begin() + vertex, to_add);
//    return this;
//}

std::set<int> Graph::get_adjacent_vertices(int vertex) {
    std::set<int> adjacencies;
    std::set<int> vertices = get_vertices();
    for(int v : vertices)
        if (are_connected(vertex, v))
            adjacencies.insert(v);
    return adjacencies;
}

//std::set<std::pair<int, int>> Graph::get_edges() {
//    std::set<std::pair<int, int>> edges;
//    for (int v1 = 0; v1 < num_vertices; v1++)
//        for (int v2 = v1 + 1; v2 < num_vertices; v2++)
//            if (are_connected(v1, v2))
//                edges.insert({v1, v2});
//    return edges;
//}

unsigned long Graph::get_num_vertices() {
    return num_vertices;
}

std::set<int> Graph::get_vertices() {
    std::set<int> vertices;
    for(unsigned int i = 0; i < num_vertices; ++i)
        vertices.insert(i);
    return vertices;
}

Graph* Graph::set_size(unsigned long size) {
    ++size;
    matrix.reserve(size);
    for(unsigned i = 0; i < size+1; ++i) {
        std::vector<bool> blank_row(size);
        std::fill(blank_row.begin(), blank_row.end(), 0);
        matrix.push_back(blank_row);
    }
    num_vertices = size;
    return this;
}

Graph* Graph::remove_vertex(int vertex) {
    matrix.erase(matrix.begin(), matrix.begin() + vertex);
    for (auto row : matrix)
        row.erase(row.begin(), row.begin() + vertex);
    --num_vertices;
    return this;
}

// Private Methods

Graph* Graph::setMatrix(int v1, int v2, bool connected) {
    matrix[v1][v2] = connected;
    matrix[v2][v1] = connected;
    return this;
}

