#include <limits.h>
#include <stdbool.h>

#include <vector>

typedef unsigned char uchar;
typedef unsigned int uint;

struct Graph {
    int n;
    std::vector<std::vector<unsigned int>> adjmat;
    std::vector<unsigned int> label;
    Graph(unsigned int n);
    Graph() {};
};

Graph induced_subgraph(struct Graph& g, std::vector<int> vv);

Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);

struct GraphData {
	Graph g;
	Graph g_sorted;
	GraphData *g0;
	std::vector<int> map_g0;
	GraphData *g1;
	std::vector<int> map_g1;
	int ordine;

	GraphData() {};
	GraphData(Graph graph) {
		g = graph;
		g0 = nullptr;
		g1 = nullptr;
		ordine = -1;
	}
	GraphData(Graph graph, GraphData* parent_g0, GraphData* parent_g1) {
		g = graph;
		g0 = parent_g0;
		map_g0.resize(g0->g.n);
		g1 = parent_g1;
		map_g1.resize(g1->g.n);
		ordine = -1;
	}
};
