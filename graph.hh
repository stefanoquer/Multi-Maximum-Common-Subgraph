#include <limits.h>
#include <stdbool.h>
#include <functional>

#include <vector>
#include <string>

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <fstream>

struct Graph {
    std::string file_name;
    int n;
    std::vector<std::vector<unsigned int>> adjmat;
    std::vector<unsigned int> label;
    Graph(unsigned int n);
};

struct SolutionGraph {
    SolutionGraph* parent = nullptr;
    const Graph* g0;
    const Graph* g1;
    std::vector<int> map_g0;
    std::vector<int> map_g1;
    Graph* g = nullptr;
    SolutionGraph(unsigned int n, const Graph *g_v, const Graph *g_w) {
        g = new Graph(n);
        map_g0.resize(n);
        map_g1.resize(n);
        g0 = g_v;
        g1 = g_w;
    }
    ~SolutionGraph() {
        if (g != nullptr) {
            delete g;
            g = nullptr;
        }
    }
    SolutionGraph() {
        std::cout << "Non deve succedere, ma mi serve per allocare un vettore..." << std::endl;
    }
};

SolutionGraph* copy_solution(SolutionGraph* sg);

Graph induced_subgraph(struct Graph& g, std::vector<int> vv);

Graph readGraph(const char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);

