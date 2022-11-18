#include "App.h"

#include <iostream>
#include <chrono>

using namespace std;

App::App(const std::string& desc, const std::string& filename) : CLI::App(desc, filename)
{
    add_flag("-q, --quiet", arguments.quiet, "Quiet output")->required(false);
    add_flag("-v, --verbose", arguments.verbose, "Verbose output")->required(false);
    auto format = this->add_option_group("Format");
    format->add_flag("-d, --dimacs", arguments.dimacs, "Read DIMACS format");
    format->add_flag("-l, --ladder", arguments.lad, "Read LAD format");
    format->add_flag("-b, --bin", arguments.bin, "Read binary format");
    format->add_flag("-e, --bin_e", arguments.bin_enrico, "Read binary alternate format");
    format->add_flag("-o, --ioi", arguments.ioi, "Read IOI format");
    format->required(true);
    add_flag("-c, --connected", arguments.connected, "Solve max common CONNECTED subgraph problem");
    add_flag("-i, --directed", arguments.directed, "Use directed graphs");
    add_flag("-a, --labelled", arguments.edge_labelled, "Use edge and vertex labels");
    add_flag("-x, --vertex-labelled-only", arguments.vertex_labelled, "Use vertex labels, but not edge labels");
    add_flag("-b, --big-first", arguments.big_first, "First try to find an induced subgraph isomorphism, then decrement the target size");
    add_option("-t, --timeout", arguments.timeout, "Specify a timeout (seconds)");
    add_option("-T, --thread", arguments.threads, "Specify how many threads to use");
    auto heuristic = this->add_option_group("Heuristic");
    heuristic->add_flag("-m, --min_max", _set_minmax, "Use min_max heuristic");
    heuristic->add_flag("-p, --min_product", _set_minproduct, "Use min_product heuristic");
    heuristic->required(true);
    add_option("files", arguments.filenames, "At least 2 input files")->check(CLI::ExistingFile);

    allow_config_extras(CLI::config_extras_mode::ignore);
}

App::~App() = default;

void App::run() {
    cout << "Checking files ... " << endl;
    arguments.arg_num = arguments.filenames.size();
    if (arguments.arg_num < 2) {
        cerr << "Not enough input files! " << arguments.arg_num << " args out of a minimum of 2!" << endl;
    }
    if (arguments.edge_labelled) {
        arguments.vertex_labelled = true;
    }
    if (_set_minmax) {
        arguments.heuristic = Heuristic::min_max;
    }
    else if (_set_minproduct) {
        arguments.heuristic = Heuristic::min_product;
    }
    cout << "Starting ... " << endl;
    
    mcsp m;
    m.start(arguments);
}