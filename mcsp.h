#include <vector>
#include <string>
#include <thread>

enum Heuristic { min_max, min_product };

struct args {
    bool quiet;
    bool verbose;
    bool dimacs;
    bool lad;
    bool bin;
    bool bin_enrico;
    bool ioi;
    bool connected;
    bool directed;
    bool edge_labelled;
    bool vertex_labelled;
    bool big_first;
    Heuristic heuristic;
    std::vector<std::string> filenames;
    int timeout;
    int threads;
    int arg_num;
};

class mcsp {
public:
    void start(args &arguments);
};