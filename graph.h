#ifndef GRAPH_H
#define GRAPH_H

#define DEBUG 0
//#define GRAPH_MAX_DEF 40
//#define MAX_N_DEF 100
#define N_THREADS std::thread::hardware_concurrency()

#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <vector>



typedef unsigned long long ULL;

struct Graph {
    int n;
    std::vector<std::vector<unsigned char>> adjmat;
    std::vector<unsigned int> label;

    Graph(int size) {
        adjmat.resize(size, std::vector<unsigned char>(size, 0));

        label.resize(size, 0);
        n = size;
    }
};

// Precondition: *g is already zeroed out
Graph *readGraph(char* filename, char format);

// Precondition: *g is already zeroed out
Graph *readBinaryGraph(char* filename);

// Precondition: *g is already zeroed out
Graph *readLadGraph(char* filename);

#endif