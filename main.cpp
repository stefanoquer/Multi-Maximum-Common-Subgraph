#include "App.h"

#include <iostream>
#include <vector>
#include "CLI11.hpp"

using namespace std;

void set_default_arguments(App &app) {
    app.arguments.quiet = false;
    app.arguments.verbose = false;
    app.arguments.dimacs = false;
    app.arguments.lad = false;
    app.arguments.bin_enrico = false;
    app.arguments.ioi = false;
    app.arguments.connected = false;
    app.arguments.directed = false;
    app.arguments.edge_labelled = false;
    app.arguments.vertex_labelled = false;
    app.arguments.big_first = false;
    app.arguments.timeout = 0;
    app.arguments.threads = std::thread::hardware_concurrency();
    app.arguments.arg_num = 0;
    app.arguments.heuristic = min_max;
}

int main(const int argc, const char* argv[]) {
    //try {
        App app("Multiple graphs common subgraph", "mcsp");
        set_default_arguments(app);
        CLI11_PARSE(app, argc, argv);
        app.run();
    //}
    /*catch (std::exception &e)
    {
        cerr << e.what() << endl;
        return -1;
    }*/
	return 0;
}