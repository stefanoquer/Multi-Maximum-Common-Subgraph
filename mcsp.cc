#include "graph.hh"

#include <algorithm>
#include <numeric>
#include <chrono>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <atomic>
#include <map>
#include <list>
#include <cassert>
#include <functional>
#include <array>

#include <argp.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

using std::vector;
using std::cout;
using std::endl;

using std::chrono::steady_clock;
using std::chrono::milliseconds;
using std::chrono::duration_cast;

using namespace std;

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

enum Heuristic { min_max, min_product };

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Find a maximum clique in a graph in DIMACS format\vHEURISTIC can be min_max or min_product";
static char args_doc[] = "HEURISTIC FILENAME1 FILENAME2";
static struct argp_option options[] = {
    {"quiet", 'q', 0, 0, "Quiet output"},
    {"verbose", 'v', 0, 0, "Verbose output"},
    {"dimacs", 'd', 0, 0, "Read DIMACS format"},
    {"lad", 'l', 0, 0, "Read LAD format"},
    {"connected", 'c', 0, 0, "Solve max common CONNECTED subgraph problem"},
    {"directed", 'i', 0, 0, "Use directed graphs"},
    {"labelled", 'a', 0, 0, "Use edge and vertex labels"},
    {"vertex-labelled-only", 'x', 0, 0, "Use vertex labels, but not edge labels"},
    {"big-first", 'b', 0, 0, "First try to find an induced subgraph isomorphism, then decrement the target size"},
    {"timeout", 't', "timeout", 0, "Specify a timeout (seconds)"},
    {"threads", 'T', "threads", 0, "Specify how many threads to use"},
    { 0 }
};

static struct {
    bool quiet;
    bool verbose;
    bool dimacs;
    bool lad;
    bool connected;
    bool directed;
    bool edge_labelled;
    bool vertex_labelled;
    bool big_first;
    Heuristic heuristic;
    vector<char*> filenames;
    int timeout;
    int threads;
    int arg_num;
} arguments;

static std::atomic<bool> abort_due_to_timeout;
atomic<int> global_position{ 0 };

void set_default_arguments() {
    arguments.quiet = false;
    arguments.verbose = false;
    arguments.dimacs = false;
    arguments.lad = false;
    arguments.connected = false;
    arguments.directed = false;
    arguments.edge_labelled = false;
    arguments.vertex_labelled = false;
    arguments.big_first = false;
    arguments.timeout = 0;
    arguments.threads = std::thread::hardware_concurrency();
    arguments.arg_num = 0;
    arguments.heuristic = min_max;
}

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    switch (key) {
        case 'd':
            if (arguments.lad)
                fail("The -d and -l options cannot be used together.\n");
            arguments.dimacs = true;
            break;
        case 'l':
            if (arguments.dimacs)
                fail("The -d and -l options cannot be used together.\n");
            arguments.lad = true;
            break;
        case 'q':
            arguments.quiet = true;
            break;
        case 'v':
            arguments.verbose = true;
            break;
        case 'c':
            if (arguments.directed)
                fail("The connected and directed options can't be used together.");
            arguments.connected = true;
            break;
        case 'i':
            if (arguments.connected)
                fail("The connected and directed options can't be used together.");
            arguments.directed = true;
            break;
        case 'a':
            if (arguments.vertex_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.edge_labelled = true;
            arguments.vertex_labelled = true;
            break;
        case 'x':
            if (arguments.edge_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.vertex_labelled = true;
            break;
        case 'b':
            arguments.big_first = true;
            break;
        case 't':
            arguments.timeout = std::stoi(arg);
            break;
        case 'T':
            arguments.threads = std::stoi(arg);
            break;
        case ARGP_KEY_ARG:
            if (arguments.arg_num == 0) {
                if (std::string(arg) == "min_max")
                    arguments.heuristic = min_max;
                else if (std::string(arg) == "min_product")
                    arguments.heuristic = min_product;
                else
                    fail("Unknown heuristic (try min_max or min_product)");
            } else {
                arguments.filenames.push_back(arg);
            }
            arguments.arg_num++;
            break;
        case ARGP_KEY_END:
            if (arguments.arg_num == 0)
                argp_usage(state);
            break;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

/*******************************************************************************
                                 MCS functions
*******************************************************************************/

struct VtxSet {
    vector<int> vv;
    VtxSet(vector<int> v) : vv(v) {};
};

struct Multidomain {
    vector<int> sets;
    vector<int> len;
    bool is_adjacent;

    Multidomain(vector<int> sets, vector<int> len, bool adj) :
        sets(sets),             // starting indices
        len(len),               // lengths
        is_adjacent(adj) { };   // is adjacent to last node
};

struct AtomicIncumbent
{
    std::atomic<unsigned> value;

    AtomicIncumbent()
    {
        value.store(0, std::memory_order_seq_cst);
    }

    bool update(unsigned v)
    {
        while (true) {
            unsigned cur_v = value.load(std::memory_order_seq_cst);
            if (v > cur_v) {
                if (value.compare_exchange_strong(cur_v, v, std::memory_order_seq_cst)) {
                    return true;
                }
            }
            else
                return false;
        }
    }
};

using PerThreadIncumbents = std::map<std::thread::id, vector<VtxSet> >;

const constexpr int split_levels = 4;

struct Position
{
    std::array<unsigned, split_levels + 1> values;
    unsigned depth;

    Position()
    {
        std::fill(values.begin(), values.end(), 0);
        depth = 0;
    }

    bool operator< (const Position & other) const
    {
        if (depth < other.depth)
            return true;
        else if (depth > other.depth)
            return false;

        for (unsigned p = 0 ; p < split_levels + 1 ; ++p)
            if (values.at(p) < other.values.at(p))
                return true;
            else if (values.at(p) > other.values.at(p))
                return false;

        return false;
    }

    void add(unsigned d, unsigned v)
    {
        depth = d;
        if (d <= split_levels)
            values[d] = v;
    }
};

struct HelpMe
{
    struct Task
    {
        const std::function<void (unsigned long long &)> * func;
        int pending;
    };

    std::mutex general_mutex;
    std::condition_variable cv;
    std::map<Position, Task> tasks;
    std::atomic<bool> finish;

    vector<std::thread> threads;

    std::list<milliseconds> times;
    std::list<unsigned long long> nodes;

    HelpMe(int n_threads) :
        finish(false)
    {
        for (int t = 0 ; t < n_threads ; ++t)
            threads.emplace_back([this, n_threads, t] {
                    milliseconds total_work_time = milliseconds::zero();
                    unsigned long long this_thread_nodes = 0;
                    while (! finish.load()) {
                        std::unique_lock<std::mutex> guard(general_mutex);
                        bool did_something = false;
                        for (auto task = tasks.begin() ; task != tasks.end() ; ++task) {
                            if (task->second.func) {
                                auto f = task->second.func;
                                ++task->second.pending;
                                guard.unlock();

                                auto start_work_time = steady_clock::now(); // local start time

                                (*f)(this_thread_nodes);

                                auto work_time = duration_cast<milliseconds>(steady_clock::now() - start_work_time);
                                total_work_time += work_time;

                                guard.lock();
                                task->second.func = nullptr;
                                if (0 == --task->second.pending)
                                    cv.notify_all();

                                did_something = true;
                                break;
                            }
                        }

                        if ((! did_something) && (! finish.load()))
                            cv.wait(guard);
                    }

                    std::unique_lock<std::mutex> guard(general_mutex);
                    times.push_back(total_work_time);
                    nodes.push_back(this_thread_nodes);
                    });
    }

    auto kill_workers() -> void
    {
        {
            std::unique_lock<std::mutex> guard(general_mutex);
            finish.store(true);
            cv.notify_all();
        }

        for (auto & t : threads)
            t.join();

        threads.clear();

        if (! times.empty()) {
            cout << "Thread work times";
            for (auto & t : times)
                cout << " " << t.count();
            cout << endl;
            times.clear();
        }
    }

    ~HelpMe()
    {
        kill_workers();
    }

    HelpMe(const HelpMe &) = delete;

    void get_help_with(
            const Position & position,
            const std::function<void (unsigned long long &)> & main_func,
            const std::function<void (unsigned long long &)> & thread_func,
            unsigned long long & main_nodes)
    {
        std::map<Position, HelpMe::Task>::iterator task;
        {
            std::unique_lock<std::mutex> guard(general_mutex);
            auto r = tasks.emplace(position, HelpMe::Task{ &thread_func, 0 });
            assert(r.second);
            task = r.first;
            cv.notify_all();
        }

        main_func(main_nodes);

        {
            std::unique_lock<std::mutex> guard(general_mutex);
            while (0 != task->second.pending)
                cv.wait(guard);
            tasks.erase(task);
        }
    }
};

bool check_sol(const vector<Graph> & g, const vector<VtxSet> & solution) {

    for (int ng = 1; ng < arguments.arg_num; ng++) {
        vector<bool> used_left(g[0].n, false);
        vector<bool> used_right(g[ng].n, false);

        for (int i = 0; i < solution[0].vv.size(); i++) {
            struct VtxSet p0 = solution[i];

            if (used_left[p0.vv[0]] || used_right[p0.vv[ng]])
                return false;

            used_left[p0.vv[0]] = true;
            used_right[p0.vv[ng]] = true;

            if (g[0].label[p0.vv[0]] != g[ng].label[p0.vv[ng]])
                return false;

            for (int j = i + 1; j < solution[0].vv.size(); j++) {
                struct VtxSet p1 = solution[j];
                if (g[0].adjmat[p0.vv[0]][p1.vv[0]] != g[ng].adjmat[p0.vv[ng]][p1.vv[ng]])
                    return false;
            }
        }
    }

    return true;
}

int calc_bound(const vector<Multidomain>& domains) {
    int bound = 0;
    for (const Multidomain &bd : domains) {
        bound += *min_element(bd.len.begin(), bd.len.end());
    }
    return bound;
}

int find_min_value(const vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i=0; i<len; i++)
        if (arr[start_idx + i] < min_v)
            min_v = arr[start_idx + i];
    return min_v;
}

int select_bidomain(const vector<Multidomain>& domains, const vector<int> & left,
        int current_matching_size)
{
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i=0; i<domains.size(); i++) {
        const Multidomain &bd = domains[i];
        if (arguments.connected && current_matching_size>0 && !bd.is_adjacent) continue;
        int len = arguments.heuristic == min_max ?
                *max_element(bd.len.begin(), bd.len.end()) :
                accumulate(bd.len.begin(), bd.len.end(), 1, std::multiplies<int>{});
        if (len < min_size) {
            min_size = len;
            min_tie_breaker = find_min_value(left, bd.sets[0], bd.len[0]);
            best = i;
        } else if (len == min_size) {
            int tie_breaker = find_min_value(left, bd.sets[0], bd.len[0]);
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;
}

// Returns length of left half of array
int partition(vector<int>& all_vv, int start, int len, const vector<unsigned int> & adjrow) {
    int i=0;
    for (int j=0; j<len; j++) {
        if (adjrow[all_vv[start+j]]) {
            std::swap(all_vv[start+i], all_vv[start+j]);
            i++;
        }
    }
    return i;
}

bool check_greater(vector<int>& lower, vector<int>& greater) {
    bool ret_val = true;
    for (int i = 0; i < lower.size(); i++) {
        ret_val &= greater[i] > lower[i];
    }
    return ret_val;
}

int max_elem(vector<unsigned int>& vet) {
    unsigned int max = 0;
    int counter = 0;
    for (int i = 0; i < arguments.arg_num; i++) {
        if (vet[i] > max) {
            max = vet[i];
            counter = 1;
        }
        else if (vet[i] == max) {
            counter++;
        }
    }
    if (counter == arguments.arg_num) {
        return -1;
    }
    return max;
}

// multiway is for directed and/or labelled graphs
vector<Multidomain> filter_domains(const vector<Multidomain>& d, vector<vector<int>> &vv,
    const vector<Graph>& g, vector<int>& vertex,
    bool multiway)
{
    vector<Multidomain> new_d;
    new_d.reserve(d.size());
    for (const Multidomain& old_bd : d) {
        vector<int> sets = old_bd.sets;
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        vector<int> len_edge(arguments.arg_num);
        vector<int> len_noedge(arguments.arg_num);
        for (int i = 0; i < vertex.size(); i++) {
            len_edge[i] = partition(vv[i], sets[i], old_bd.len[i], g[i].adjmat[vertex[i]]);
            len_noedge[i] = old_bd.len[i] - len_edge[i];
        }

        if (accumulate(len_noedge.begin(), len_noedge.end(), 1, multiplies<int>{})) {
            vector<int> new_d_sets (arguments.arg_num);
            transform(len_edge.begin(), len_edge.end(), sets.begin(), new_d_sets.begin(), std::plus<int>());
            new_d.push_back({ new_d_sets, len_noedge, old_bd.is_adjacent });
        }
        if (multiway && accumulate(len_edge.begin(), len_edge.end(), 1, multiplies<int>{})) {
            vector<const vector<unsigned int> *> adjrows(arguments.arg_num);
            vector<int> top(arguments.arg_num);
            for (int i = 0; i < arguments.arg_num; i++) {
                adjrows[i] = &g[i].adjmat[vertex[i]];
                vector<int>::iterator _begin = std::begin(vv[i]) + sets[i];
                std::sort(_begin, _begin + len_edge[i], [&](int a, int b)
                    { return adjrows.at(i)->at(a) < adjrows.at(i)->at(b); } );
                top[i] = sets[i] + len_edge[i];
            }
            while (check_greater(sets, top)) {
                vector<unsigned int> labels(arguments.arg_num);
                for (int i = 0; i < vertex.size(); i++) {
                    labels[i] = adjrows.at(i)->at(vv[i][sets[i]]);
                }
                int maximum = max_elem(labels);
                if (maximum != -1) {
                    for (int i = 0; i < arguments.arg_num; i++) {
                        if (labels[i] != maximum) {
                            sets[i]++;
                        }
                    }
                }
                else {
                    vector<int> min_sets = sets;
                    for (int i = 0; i < vertex.size(); i++) {
                        do { sets[i]++; } while (sets[i] < top[i] && adjrows.at(i)->at(vv[i][sets[i]]) == labels[0]);
                    }
                    vector<int> dif_sets(arguments.arg_num);
                    transform(sets.begin(), sets.end(), min_sets.begin(), dif_sets.begin(), std::minus<int>());
                    new_d.push_back({ min_sets, dif_sets, true });
                }
            }
        }
        else if (accumulate(len_edge.begin(), len_edge.end(), 1, multiplies<int>{})) {
            new_d.push_back({ sets, len_edge, true });
        }
    }
    return new_d;
}

// returns the index of the smallest value in arr that is >w.
// Assumption: such a value exists
// Assumption: arr contains no duplicates
// Assumption: arr has no values==INT_MAX
int index_of_next_smallest(const vector<int>& arr, int start_idx, int len, int w) {
    int idx = -1;
    int smallest = INT_MAX;
    for (int i=0; i<len; i++) {
        if (arr[start_idx + i]>w && arr[start_idx + i]<smallest) {
            smallest = arr[start_idx + i];
            idx = i;
        }
    }
    return idx;
}

void remove_vtx_from_domain(vector<int>& left, Multidomain& bd, int v, int idx)
{
    int i = 0;
    while(left[bd.sets[idx] + i] != v) i++;
    std::swap(left[bd.sets[idx]+i], left[bd.sets[idx]+bd.len[idx]-1]);
    bd.len[idx]--;
}

void remove_bidomain(vector<Multidomain>& domains, int idx) {
    domains[idx] = domains[domains.size()-1];
    domains.pop_back();
}

void show(vector<VtxSet>& current, vector<Multidomain>& domains, vector<vector<int>>& vv) {
    printf("Length of current assignment: %d\n", current.size());
    printf("Current assignment:");

    for (int i = 0; i < current.size(); i++) {
        for (int j = 0; j < arguments.arg_num; j++) {
            if (j == 0)
                printf("  %d", current[i].vv[j]);
            else
                printf("->%d", current[i].vv[j]);
        }
    }
    printf("\n");
    for (int i = 0; i < domains.size(); i++) {
        struct Multidomain bd = domains[i];
        for (int ng = 0; ng < arguments.arg_num; ng++) {
            printf("Graph %d  ", ng);
            for (int j = 0; j < bd.len[ng]; j++)
                printf("%d ", vv[ng][j+bd.sets[ng]]);
            printf("\n");
        }
    }
    printf("\n\n");
}

void solve_nopar(const unsigned depth, vector<Graph>& g,
    AtomicIncumbent& global_incumbent,
    vector<VtxSet>& my_incumbent,
    vector<VtxSet>& current, vector<Multidomain>& domains,
    vector<vector<int>>& vv, const unsigned int matching_size_goal,
    unsigned long long& my_thread_nodes);
void solve(const unsigned depth, vector<Graph>& g,
    AtomicIncumbent& global_incumbent,
    PerThreadIncumbents& per_thread_incumbents,
    vector<VtxSet>& current, vector<Multidomain>& domains,
    vector<vector<int>>& vv, const unsigned int matching_size_goal,
    const Position& position, HelpMe& help_me, unsigned long long& my_thread_nodes);

void solve_nopar_recursive(Multidomain& bd, std::vector<std::vector<int>>& vv, std::vector<Multidomain>& domains,
    std::vector<Graph>& g, std::vector<VtxSet>& current, const unsigned int& depth,
    AtomicIncumbent& global_incumbent, std::vector<VtxSet>& my_incumbent,
    const unsigned int& matching_size_goal, unsigned long long& my_thread_nodes,
    int bd_idx, vector<int> &nodi_inseriti, int n_nodi_inseriti)
{
    int w = -1;
    const int i_end = bd.len[n_nodi_inseriti] + 2; /* including the null */

    for (int i = 0; i < i_end /* not != */; i++) {
        if (i != i_end - 1) {
            int idx = index_of_next_smallest(vv[n_nodi_inseriti], bd.sets[n_nodi_inseriti], bd.len[n_nodi_inseriti] + 1, w);
            w = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx];

            // swap w to the end of its colour class
            vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx] = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]];
            vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]] = w;

            nodi_inseriti.push_back(w);

            if (n_nodi_inseriti == arguments.arg_num - 1) {
                auto new_domains = filter_domains(domains, vv, g, nodi_inseriti, arguments.directed || arguments.edge_labelled);
                
                current.push_back(VtxSet(nodi_inseriti));
                solve_nopar(depth + 1, g, global_incumbent, my_incumbent, current, new_domains, vv, matching_size_goal, my_thread_nodes);
                current.pop_back();
            }
            else {
                auto new_domains = domains;
                auto &new_bd = new_domains[bd_idx];
                solve_nopar_recursive(new_bd, vv, new_domains, g, current, depth, global_incumbent, my_incumbent, matching_size_goal, my_thread_nodes, bd_idx, nodi_inseriti, n_nodi_inseriti + 1);
            }
            nodi_inseriti.pop_back();
        }
        else {
            // Last assign is null. Keep it in the loop to simplify parallelism.
            if (n_nodi_inseriti == 1) {
                transform(bd.len.begin() + 1, bd.len.end(), bd.len.begin() + 1, [](int x) { return x + 1; });
                if (bd.len[0] == 0)
                    remove_bidomain(domains, bd_idx);

                solve_nopar(depth + 1, g, global_incumbent, my_incumbent, current, domains, vv, matching_size_goal, my_thread_nodes);
            }
        }
    }
}

void solve_nopar(const unsigned depth, vector<Graph> & g,
        AtomicIncumbent & global_incumbent,
        vector<VtxSet> & my_incumbent,
        vector<VtxSet> & current, vector<Multidomain> & domains,
        vector<vector<int>> &vv,
        const unsigned int matching_size_goal,
        unsigned long long & my_thread_nodes)
{
    if (arguments.verbose) show(current, domains, vv);

    if (abort_due_to_timeout)
        return;

    my_thread_nodes++;

    if (my_incumbent.size() < current.size()) {
        my_incumbent = current;
        global_incumbent.update(current.size());
    }

    unsigned int bound = current.size() + calc_bound(domains);
    if (bound <= global_incumbent.value || bound < matching_size_goal)
        return;

    if (arguments.big_first && global_incumbent.value == matching_size_goal)
        return;

    int bd_idx = select_bidomain(domains, vv[0], current.size());
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Multidomain &bd = domains[bd_idx];

    transform(bd.len.begin() + 1, bd.len.end(), bd.len.begin() + 1, [](int x) { return x - 1; });
    std::atomic<int> shared_i{ 0 };

    int v = find_min_value(vv[0], bd.sets[0], bd.len[0]);
    remove_vtx_from_domain(vv[0], domains[bd_idx], v, 0);

    vector<int> nodi_inseriti;
    nodi_inseriti.reserve(arguments.arg_num);
    nodi_inseriti.push_back(v);

    solve_nopar_recursive(bd, vv, domains, g, current, depth, global_incumbent, my_incumbent, matching_size_goal, my_thread_nodes, bd_idx, nodi_inseriti, 1);
}

void solve_recursive(const int& i_end, std::vector<std::vector<int>>& vv, Multidomain& bd, int which_i_should_i_run_next,
                     std::atomic_int& shared_i, std::vector<Multidomain>& domains, std::vector<Graph>& g,
                     std::vector<VtxSet>& current, const unsigned int& depth, AtomicIncumbent& global_incumbent,
                     PerThreadIncumbents& per_thread_incumbents, const unsigned int& matching_size_goal,
                     unsigned long long& main_thread_nodes, const Position& position, HelpMe& help_me,
                     int bd_idx, vector<int>& nodi_inseriti, int n_nodi_inseriti)
{
    int w = -1;

    for (int i = 0; i < i_end /* not != */; i++) {
        if (i != i_end - 1) {
            int idx = index_of_next_smallest(vv[n_nodi_inseriti], bd.sets[n_nodi_inseriti], bd.len[n_nodi_inseriti] + 1, w);
            w = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx];

            // swap w to the end of its colour class
            vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx] = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]];
            vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]] = w;


            nodi_inseriti.push_back(w);
            if (n_nodi_inseriti == arguments.arg_num - 1) {
                auto new_domains = filter_domains(domains, vv, g, nodi_inseriti, arguments.directed || arguments.edge_labelled);
                
                current.push_back(VtxSet(nodi_inseriti));
                if (depth > split_levels) {
                    solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, new_domains, vv, matching_size_goal, main_thread_nodes);
                }
                else {
                    auto new_position = position;
                    new_position.add(depth, ++global_position);
                    solve(depth + 1, g, global_incumbent, per_thread_incumbents, current, new_domains, vv, matching_size_goal, new_position, help_me, main_thread_nodes);
                }
                current.pop_back();
            }
            else {
                auto new_domains = domains;
                const int i_end = bd.len[n_nodi_inseriti + 1] + 2;
                auto& new_bd = new_domains[bd_idx];
                solve_recursive(i_end, vv, new_bd, which_i_should_i_run_next, shared_i, new_domains, g, current, depth,
                                global_incumbent, per_thread_incumbents, matching_size_goal, main_thread_nodes,
                                position, help_me, bd_idx, nodi_inseriti, n_nodi_inseriti + 1);
            }
            nodi_inseriti.pop_back();
        }
        else {
            // Last assign is null. Keep it in the loop to simplify parallelism.
            if (n_nodi_inseriti == 1) {
                Multidomain tmp = domains[bd_idx];
                transform(bd.len.begin() + 1, bd.len.end(), bd.len.begin() + 1, [](int x) { return x + 1; });
                if (bd.len[0] == 0)
                    remove_bidomain(domains, bd_idx);


                if (depth > split_levels) {
                    solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, domains, vv, matching_size_goal, main_thread_nodes);
                }
                else {
                    auto new_position = position;
                    new_position.add(depth, ++global_position);
                    solve(depth + 1, g, global_incumbent, per_thread_incumbents, current, domains, vv, matching_size_goal, new_position, help_me, main_thread_nodes);
                }
            }

        }
    }
}

void solve(const unsigned depth, vector<Graph> & g,
        AtomicIncumbent & global_incumbent,
        PerThreadIncumbents & per_thread_incumbents,
        vector<VtxSet> & current, vector<Multidomain> & domains,
        vector<vector<int>> & vv,
        const unsigned int matching_size_goal,
        const Position & position, HelpMe & help_me,
        unsigned long long & my_thread_nodes)
{

    if (arguments.verbose) show(current, domains, vv);
    if (abort_due_to_timeout)
        return;

    my_thread_nodes++;

    if (per_thread_incumbents.find(std::this_thread::get_id())->second.size() < current.size()) {
        per_thread_incumbents.find(std::this_thread::get_id())->second = current;
        global_incumbent.update(current.size());
    }

    unsigned int bound = current.size() + calc_bound(domains);
    if (bound <= global_incumbent.value || bound < matching_size_goal)
        return;

    if (arguments.big_first && global_incumbent.value == matching_size_goal)
        return;

    int bd_idx = select_bidomain(domains, vv[0], current.size());
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Multidomain &bd = domains[bd_idx];

    transform(bd.len.begin() + 1, bd.len.end(), bd.len.begin() + 1, [](int x) { return x - 1; });
    std::atomic<int> shared_i{ 0 };
    const int i_end = bd.len[1] + 2; /* including the null */

    // Version of the loop used by helpers
    std::function<void (unsigned long long &)> helper_function = [&shared_i, &g, &global_incumbent, &per_thread_incumbents, &position, &depth,
        i_end, matching_size_goal, current, domains, vv, &help_me] (unsigned long long & help_thread_nodes) {
        int which_i_should_i_run_next = shared_i++;

        if (which_i_should_i_run_next >= i_end)
            return; /* don't waste time recomputing */

        /* recalculate to this point */
        vector<VtxSet> help_current = current;
        vector<Multidomain> help_domains = domains;
        vector<vector<int>> help_vv = vv;

        /* rerun important stuff from before the loop */
        int help_bd_idx = select_bidomain(help_domains, help_vv[0], help_current.size());
        if (help_bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
            return;
        Multidomain &help_bd = help_domains[help_bd_idx];

        int help_v = find_min_value(help_vv[0], help_bd.sets[0], help_bd.len[0]);
        remove_vtx_from_domain(help_vv[0], help_domains[help_bd_idx], help_v, 0);

        vector<int> nodi_inseriti;
        nodi_inseriti.reserve(arguments.arg_num);
        nodi_inseriti.push_back(help_v);

        int n_nodi_inseriti = 1;

        int help_w = -1;

        for (int i = 0; i < i_end /* not != */; i++) {
            if (i != i_end - 1) {
                int idx = index_of_next_smallest(help_vv[n_nodi_inseriti], help_bd.sets[n_nodi_inseriti], help_bd.len[n_nodi_inseriti] + 1, help_w);
                help_w = help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + idx];

                // swap w to the end of its colour class
                help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + idx] = help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + help_bd.len[n_nodi_inseriti]];
                help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + help_bd.len[n_nodi_inseriti]] = help_w;

                if (i == which_i_should_i_run_next) {
                    which_i_should_i_run_next = shared_i++;
                    nodi_inseriti.push_back(help_w);
                    if (n_nodi_inseriti == arguments.arg_num - 1) {
                        auto new_domains = filter_domains(help_domains, help_vv, g, nodi_inseriti, arguments.directed || arguments.edge_labelled);
                        
                        help_current.push_back(VtxSet(nodi_inseriti));
                        if (depth > split_levels) {
                            solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, new_domains, help_vv, matching_size_goal, help_thread_nodes);
                        }
                        else {
                            auto new_position = position;
                            new_position.add(depth, ++global_position);
                            solve(depth + 1, g, global_incumbent, per_thread_incumbents, help_current, new_domains, help_vv, matching_size_goal, new_position, help_me, help_thread_nodes);
                        }
                        help_current.pop_back();
                    }
                    else {
                        auto new_domains = help_domains;
                        const int i_end = help_bd.len[n_nodi_inseriti + 1] + 2;
                        auto &new_help_bd = new_domains[help_bd_idx];
                        solve_recursive(i_end, help_vv, new_help_bd, which_i_should_i_run_next, shared_i,
                                             new_domains, g, help_current, depth, global_incumbent,
                                             per_thread_incumbents, matching_size_goal, help_thread_nodes, position,
                                             help_me, help_bd_idx, nodi_inseriti, n_nodi_inseriti + 1);
                    }
                    nodi_inseriti.pop_back();
                }
            }
            else {
                // Last assign is null. Keep it in the loop to simplify parallelism.
                if (n_nodi_inseriti == 1) {
                    transform(help_bd.len.begin() + 1, help_bd.len.end(), help_bd.len.begin() + 1, [](int x) { return x + 1; });
                    if (help_bd.len[0] == 0)
                        remove_bidomain(help_domains, help_bd_idx);

                    if (i == which_i_should_i_run_next) {
                        which_i_should_i_run_next = shared_i++;
                        if (depth > split_levels) {
                            solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, help_domains, help_vv, matching_size_goal, help_thread_nodes);
                        }
                        else {
                            auto new_position = position;
                            new_position.add(depth, ++global_position);
                            solve(depth + 1, g, global_incumbent, per_thread_incumbents, help_current, help_domains, help_vv, matching_size_goal, new_position, help_me, help_thread_nodes);
                        }
                    }
                }
            }
        }
    };

    // Grab this first, before advertising that we can get help
    int which_i_should_i_run_next = shared_i++;

    // Version of the loop used by the main thread
    std::function<void (unsigned long long &)> main_function = [&] (unsigned long long & main_thread_nodes) {
        int v = find_min_value(vv[0], bd.sets[0], bd.len[0]);
        remove_vtx_from_domain(vv[0], domains[bd_idx], v, 0);

        vector<int> nodi_inseriti;
        nodi_inseriti.reserve(arguments.arg_num);
        nodi_inseriti.push_back(v);

        int n_nodi_inseriti = 1;

        int w = -1;

        for (int i = 0; i < i_end /* not != */; i++) {
            if (i != i_end - 1) {

                int idx = index_of_next_smallest(vv[n_nodi_inseriti], bd.sets[n_nodi_inseriti], bd.len[n_nodi_inseriti] + 1, w);
                w = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx];

                // swap w to the end of its colour class
                vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx] = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]];
                vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]] = w;

                if (i == which_i_should_i_run_next) {

                    nodi_inseriti.push_back(w);
                    which_i_should_i_run_next = shared_i++;
                    if (n_nodi_inseriti == arguments.arg_num - 1) {
                        auto new_domains = filter_domains(domains, vv, g, nodi_inseriti, arguments.directed || arguments.edge_labelled);
                        current.push_back(VtxSet(nodi_inseriti));
                        if (depth > split_levels) {
                            solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, new_domains, vv, matching_size_goal, main_thread_nodes);
                        }
                        else {
                            auto new_position = position;
                            new_position.add(depth, ++global_position);
                            solve(depth + 1, g, global_incumbent, per_thread_incumbents, current, new_domains, vv, matching_size_goal, new_position, help_me, main_thread_nodes);
                        }
                        current.pop_back();
                    }
                    else {
                        auto new_domains = domains;
                        const int i_end = bd.len[n_nodi_inseriti + 1] + 2;
                        auto &new_bd = new_domains[bd_idx];
                        solve_recursive(i_end, vv, new_bd, which_i_should_i_run_next, shared_i, new_domains, g, current,
                                        depth, global_incumbent, per_thread_incumbents, matching_size_goal,
                                        main_thread_nodes, position, help_me, bd_idx, nodi_inseriti,
                                        n_nodi_inseriti + 1);
                    }
                    nodi_inseriti.pop_back();
                }
            }
            else {
                // Last assign is null. Keep it in the loop to simplify parallelism.
                if (n_nodi_inseriti == 1) {
                    Multidomain tmp = domains[bd_idx];
                    transform(bd.len.begin() + 1, bd.len.end(), bd.len.begin() + 1, [](int x) { return x + 1; });
                    if (bd.len[0] == 0)
                        remove_bidomain(domains, bd_idx);

                    if (i == which_i_should_i_run_next) {
                        which_i_should_i_run_next = shared_i++;
                        if (depth > split_levels) {
                            solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, domains, vv, matching_size_goal, main_thread_nodes);
                        }
                        else {
                            auto new_position = position;
                            new_position.add(depth, ++global_position);
                            solve(depth + 1, g, global_incumbent, per_thread_incumbents, current, domains, vv, matching_size_goal, new_position, help_me, main_thread_nodes);
                        }
                    }
                }

            }
        }
    };

    if (depth <= split_levels)
        help_me.get_help_with(position, main_function, helper_function, my_thread_nodes);
    else
        main_function(my_thread_nodes);
}

std::set<unsigned int> intersection(const vector<set<unsigned int>>& vecs) {

    auto last_intersection = vecs[0];
    set<unsigned int> curr_intersection;

    for (std::size_t i = 1; i < vecs.size(); ++i) {
        std::set_intersection(last_intersection.begin(), last_intersection.end(),
            vecs[i].begin(), vecs[i].end(),
            std::inserter(curr_intersection, std::begin(curr_intersection)));
        std::swap(last_intersection, curr_intersection);
        curr_intersection.clear();
    }
    return last_intersection;
}

std::pair<vector<VtxSet>, unsigned long long> mcs(vector<Graph> & gi) {
    
    // the buffer of vertex indices for the partitions
    vector<vector<int>> vtx_buf(arguments.arg_num);

    auto domains = vector<Multidomain>{};

    vector<set<unsigned int>> labels_vv (arguments.arg_num);
    for (int i = 0; i < arguments.arg_num; i++) {
        for (unsigned int label : gi[i].label) {
            labels_vv[i].insert(label);
        }
    }

    std::set<unsigned int> labels = intersection(labels_vv);
    
    // Create a bidomain for each label that appears in both graphs
    for (unsigned int label : labels) {
        vector<int> starts(arguments.arg_num);
        for (int i = 0; i < arguments.arg_num; i++) {
            starts[i] = vtx_buf[i].size();
            for (int j = 0; j < gi[i].n; j++) {
                if (gi[i].label[j] == label) {
                    vtx_buf[i].push_back(j);
                }
            }
        }

        vector<int> len(arguments.arg_num);
        for (int i = 0; i < arguments.arg_num; i++) {
            len[i] = vtx_buf[i].size() - starts[i];
        }
        domains.push_back({ starts, len, false });
    }

    AtomicIncumbent global_incumbent;
    vector<VtxSet> incumbent;
    unsigned long long global_nodes = 0;

    if (arguments.big_first) {
        for (int k=0; k<gi[0].n; k++) {
            unsigned int goal = gi[0].n - k;
            auto vtx_buf_copy = vtx_buf;
            auto domains_copy = domains;
            vector<VtxSet> current;
            PerThreadIncumbents per_thread_incumbents;
            per_thread_incumbents.emplace(std::this_thread::get_id(), vector<VtxSet>());
            Position position;
            HelpMe help_me(arguments.threads - 1);
            for (auto & t : help_me.threads)
                per_thread_incumbents.emplace(t.get_id(), vector<VtxSet>());
            solve(0, gi, global_incumbent, per_thread_incumbents, current, domains_copy, vtx_buf_copy, goal, position, help_me, global_nodes);
            help_me.kill_workers();
            for (auto & n : help_me.nodes) {
                global_nodes += n;
            }
            for (auto & i : per_thread_incumbents)
                if (i.second.size() > incumbent.size())
                    incumbent = i.second;
            if (global_incumbent.value == goal || abort_due_to_timeout) break;
            if (!arguments.quiet) cout << "Upper bound: " << goal-1 << std::endl;
        }

    } else {
        vector<VtxSet> current;
        PerThreadIncumbents per_thread_incumbents;
        per_thread_incumbents.emplace(std::this_thread::get_id(), vector<VtxSet>());
        Position position;
        HelpMe help_me(arguments.threads - 1);
        for (auto & t : help_me.threads)
            per_thread_incumbents.emplace(t.get_id(), vector<VtxSet>());
        solve(0, gi, global_incumbent, per_thread_incumbents, current, domains, vtx_buf, 1, position, help_me, global_nodes);
        help_me.kill_workers();
        for (auto & n : help_me.nodes)
            global_nodes += n;
        for (auto & i : per_thread_incumbents)
            if (i.second.size() > incumbent.size())
                incumbent = i.second;
    }

    return { incumbent, global_nodes };
}

vector<int> calculate_degrees(const Graph & g) {
    vector<int> degree(g.n, 0);
    for (int v=0; v<g.n; v++) {
        for (int w=0; w<g.n; w++) {
            unsigned int mask = 0xFFFFu;
            if (g.adjmat[v][w] & mask) degree[v]++;
            if (g.adjmat[v][w] & ~mask) degree[v]++;  // inward edge, in directed case
        }
    }
    return degree;
}

int sum(const vector<int> & vec) {
    return std::accumulate(std::begin(vec), std::end(vec), 0);
}

int main(int argc, char** argv) {
    set_default_arguments();
    argp_parse(&argp, argc, argv, 0, 0, 0);
    arguments.arg_num--;

    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';
    vector<Graph> gi;
    for (int i = 0; i < arguments.arg_num; i++) {
        gi.push_back(readGraph(arguments.filenames[i], format, arguments.directed,
            arguments.edge_labelled, arguments.vertex_labelled));
    }

    std::thread timeout_thread;
    std::mutex timeout_mutex;
    std::condition_variable timeout_cv;
    abort_due_to_timeout.store(false);
    bool aborted = false;

    if (0 != arguments.timeout) {
        timeout_thread = std::thread([&] {
                auto abort_time = steady_clock::now() + std::chrono::seconds(arguments.timeout);
                {
                    /* Sleep until either we've reached the time limit,
                     * or we've finished all the work. */
                    std::unique_lock<std::mutex> guard(timeout_mutex);
                    while (! abort_due_to_timeout.load()) {
                        if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                            /* We've woken up, and it's due to a timeout. */
                            aborted = true;
                            break;
                        }
                    }
                }
                abort_due_to_timeout.store(true);
                });
    }

    double begin = clock ();
    auto start = steady_clock::now();
    
    vector<vector<int>> gi_deg(arguments.arg_num);
    for (int i = 0; i < arguments.arg_num; i++) {
        gi_deg[i] = calculate_degrees(gi[i]);
    }

    // As implemented here, g1_dense and g0_dense are false for all instances
    // in the Experimental Evaluation section of the paper.  Thus,
    // we always sort the vertices in descending order of degree (or total degree,
    // in the case of directed graphs.  Improvements could be made here: it would
    // be nice if the program explored exactly the same search tree if both
    // input graphs were complemented.
    vector<vector<int>> vvi(arguments.arg_num);
    for (int i = 0; i < arguments.arg_num; i++) {
        vvi[i].resize(gi[i].n);
        iota(std::begin(vvi[i]), std::end(vvi[i]), 0);
        stable_sort(std::begin(vvi[i]), std::end(vvi[i]), [&](int a, int b) {
            return (gi_deg[i][a] > gi_deg[i][b]);
        });
    }

    vector<Graph> gi_sorted;
    for (int i = 0; i < arguments.arg_num; i++) {
        gi_sorted.push_back(induced_subgraph(gi[i], vvi[i]));
    }

    std::pair<vector<VtxSet>, unsigned long long> solution = mcs(gi_sorted);

    // Convert to indices from original, unsorted graphs
    for (auto& vtx_pair : solution.first) {
        for (int i = 0; i < arguments.arg_num; i++) {
            vtx_pair.vv[i] = vvi[i][vtx_pair.vv[i]];
        }
    }

    auto stop = steady_clock::now();
    auto time_elapsed = duration_cast<milliseconds>(stop - start).count();

    /* Clean up the timeout thread */
    if (timeout_thread.joinable()) {
        {
            std::unique_lock<std::mutex> guard(timeout_mutex);
            abort_due_to_timeout.store(true);
            timeout_cv.notify_all();
        }
        timeout_thread.join();
    }

    double end = clock ();

    cout << "Solution size " << solution.first.size() << std::endl;
    for (int i = 0; i < gi[0].n; i++) {
        for (unsigned int j = 0; j < solution.first.size(); j++) {
            if (solution.first[j].vv[0] == i) {
                cout << "(" << solution.first[j].vv[0];
                for (int k = 1; k < arguments.arg_num; k++) {
                    cout << " -> " << solution.first[j].vv[k];
                }
                cout << ") ";
            }
        }
    }
    cout << std::endl;

    cout << "Nodes:                      " << solution.second << endl;
    cout << "CPU time (ms):              " << time_elapsed << endl;

    fprintf (stdout, "Wall-Clock Time = %f sec\n", (double)(end-begin)/CLOCKS_PER_SEC);
    
    if (aborted)
        cout << "TIMEOUT" << endl;

    if (!check_sol(gi, solution.first))
        fail("\n\n*** Error: Invalid solution\n");
        
    printf(">>> %d - %f\n", solution.first.size(), (float)time_elapsed/1000);
}

