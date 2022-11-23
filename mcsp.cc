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
#include <condition_variable>
#include <atomic>
#include <map>
#include <list>
#include <cassert>
#include <functional>
#include <array>

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mcsp.h"

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

#define MAX_ARGS 10
#define UNSORTED 0

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static std::atomic<bool> abort_due_to_timeout;
atomic<int> global_position{ 0 };

struct args arguments;

/*******************************************************************************
                               Data structures
*******************************************************************************/

struct VtxPair {
    int vv[MAX_ARGS];
    VtxPair(int *v) {
        for(int i=0; i<arguments.arg_num; i++) {
            vv[i] = v[i];
        }
    }
    //vector<int> vv;
    //VtxPair(vector<int> v) : vv(v) {};

    /*int v;
    int w;
    VtxPair(int v, int w): v(v), w(w) {}*/
};

struct Bidomain {
    int sets[MAX_ARGS];
    int len[MAX_ARGS];
    bool is_adjacent;

    Bidomain(int *sets, int *len, bool adj)  {
        for (int i=0; i< arguments.arg_num; i++){
            this->sets[i] = sets[i];
            this->len[i] = len[i];
        }
        is_adjacent = adj;
    }

    /*vector<int> sets;
    vector<int> len;
    bool is_adjacent;

    Bidomain(vector<int> sets, vector<int> len, bool adj) :
        sets(sets),
        len(len),
        is_adjacent(adj) { };*/

    /*int l,        r;        // start indices of left and right sets
    int left_len, right_len;
    bool is_adjacent;
    Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent):
            l(l),
            r(r),
            left_len (left_len),
            right_len (right_len),
            is_adjacent (is_adjacent) { };*/
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

using PerThreadIncumbents = std::map<std::thread::id, vector<VtxPair> >;

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

/*******************************************************************************
                                    Utils
*******************************************************************************/

void show(vector<VtxPair>& current, vector<Bidomain>& domains, array<vector<int>, MAX_ARGS>& vv) {
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
        struct Bidomain bd = domains[i];
        for (int ng = 0; ng < arguments.arg_num; ng++) {
            printf("Graph %d  ", ng);
            for (int j = 0; j < bd.len[ng]; j++)
                printf("%d ", vv[ng][j+bd.sets[ng]]);
            printf("\n");
        }
    }
    printf("\n\n");
}

void string_show(vector<VtxPair>& current, vector<Bidomain>& domains, array<vector<int>, MAX_ARGS>& vv, int depth) {
    string s = "";
    //printf("Length of current assignment: %d\n", current.size());
    //printf("Current assignment:");

    for (int i = 0; i < current.size(); i++) {
        for (int j = 0; j < arguments.arg_num; j++) {
            if (j == 0) {
                s = s + to_string(current[i].vv[j]);
                //printf("  %d", current[i].vv[j])
            }
            else {
                s = s + "->" + to_string(current[i].vv[j]);
                //printf("->%d", current[i].vv[j]);
            }
        }
        s = s + " ";
    }
    s = s + ": " + to_string(depth);
    s = s + "\n";
    cout << s;
    /*printf("\n");
    for (int i = 0; i < domains.size(); i++) {
        struct Bidomain bd = domains[i];
        for (int ng = 0; ng < arguments.arg_num; ng++) {
            printf("Graph %d  ", ng);
            for (int j = 0; j < bd.len[ng]; j++)
                printf("%d ", vv[ng][j + bd.sets[ng]]);
            printf("\n");
        }
    }
    printf("\n\n");*/
}

/*******************************************************************************
                                Parallel queue
*******************************************************************************/

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

/*******************************************************************************
                                 MCS functions
*******************************************************************************/

bool check_sol(const vector<Graph> & g, const vector<VtxPair> & solution) {

    for (int ng = 1; ng < arguments.arg_num; ng++) {
        vector<bool> used_left(g[0].n, false);
        vector<bool> used_right(g[ng].n, false);

        for (int i = 0; i < arguments.arg_num; i++) {
            struct VtxPair p0 = solution[i];

            if (used_left[p0.vv[0]] || used_right[p0.vv[ng]])
                return false;

            used_left[p0.vv[0]] = true;
            used_right[p0.vv[ng]] = true;

            if (g[0].label[p0.vv[0]] != g[ng].label[p0.vv[ng]])
                return false;

            for (int j = i + 1; j < arguments.arg_num; j++) {
                struct VtxPair p1 = solution[j];
                if (g[0].adjmat[p0.vv[0]][p1.vv[0]] != g[ng].adjmat[p0.vv[ng]][p1.vv[ng]])
                    return false;
            }
        }
    }

    return true;
}

int calc_bound(const vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain &bd : domains) {
        bound += *min_element(bd.len, bd.len+arguments.arg_num);
    }
    return bound;
}

int find_min_value(int *arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i=0; i<len; i++)
        if (arr[start_idx + i] < min_v)
            min_v = arr[start_idx + i];
    return min_v;
}

int select_bidomain(const vector<Bidomain>& domains, int* left,
        int current_matching_size) // da modificare
{
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i=0; i<domains.size(); i++) {
        const Bidomain &bd = domains[i];
        if (arguments.connected && current_matching_size>0 && !bd.is_adjacent) continue;
        int len;
        switch (arguments.heuristic) {
            case min_max:
                len = *max_element(bd.len, bd.len+arguments.arg_num);
            break;
            case min_min:
                len = *min_element(bd.len, bd.len+arguments.arg_num);
            break;
            case min_sum:
                len = accumulate(bd.len, bd.len+arguments.arg_num, 0);
            break;
            case min_product:
                len = accumulate(bd.len, bd.len+arguments.arg_num, 1, std::multiplies<int>{});
            break;
            default:
                cout << "Error, not implemented heuristic!" << endl;
                exit (-1);
            break;
        }
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
int partition(int *all_vv, int start, int len, const vector<unsigned int> & adjrow) {
    int i=0;
    for (int j=0; j<len; j++) {
        if (adjrow[all_vv[start+j]]) {
            std::swap(all_vv[start+i], all_vv[start+j]);
            i++;
        }
    }
    return i;
}

bool controlla_maggiore(int *minori, int *maggiori) {
    bool ret_val = true;
    for (int i = 0; i < arguments.arg_num; i++) {
        ret_val &= maggiori[i] > minori[i];
    }
    return ret_val;
}

int min_elem(vector<unsigned int>& vet) {
    unsigned int min = UINT_MAX;
    int counter = 0;
    for (int i = 0; i < arguments.arg_num; i++) {
        if (vet[i] < min) {
            min = vet[i];
        }
    }
    return min;
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
vector<Bidomain> filter_domains(const vector<Bidomain>& d, array<vector<int>, MAX_ARGS> &vv /*left, right*/,
    const vector<Graph>& g, int *vertex /*v, w*/,
    bool multiway)
{
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain& old_bd : d) {
        int sets[MAX_ARGS];
        for(int i=0; i<arguments.arg_num; i++) {
            sets[i] = old_bd.sets [i];
        }
        //vector<int> sets = old_bd.sets; /*l, r*/
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int len_edge[MAX_ARGS]; /*left_len, right_len*/
        int len_noedge[MAX_ARGS]; /*left_len_noedge, right_len_noedge*/
        for (int i = 0; i < arguments.arg_num; i++) {
            len_edge[i] = partition(vv[i].data(), sets[i], old_bd.len[i], g[i].adjmat[vertex[i]]);
            len_noedge[i] = old_bd.len[i] - len_edge[i];
        }

        if (accumulate(len_noedge, len_noedge+arguments.arg_num, 1, multiplies<int>{})) {
            vector<int> new_vector_sets(arguments.arg_num);
            int *new_d_sets = new_vector_sets.data();
            //vector<int> new_d_sets (arguments.arg_num);
            transform(len_edge, len_edge+arguments.arg_num, sets, new_d_sets, std::plus<int>());
            new_d.push_back({ new_d_sets, len_noedge, old_bd.is_adjacent });
        }
        if (multiway && accumulate(len_edge, len_edge+arguments.arg_num, 1, multiplies<int>{})) {
            vector<const vector<unsigned int> *> adjrows(arguments.arg_num); /*adjrow_v, adjrow_w*/
            //vector<vector<int>::iterator> _begin(arguments.arg_num); /*l_begin, r_begin*/
            int top[MAX_ARGS]; /*l_top, r_top*/
            for (int i = 0; i < arguments.arg_num; i++) {
                adjrows[i] = &g[i].adjmat[vertex[i]];
                //_begin[i] = std::begin(vv[i]);
                int* _begin = vv[i].data() + sets[i];
                std::sort(_begin, _begin + len_edge[i], [&](int a, int b)
                    { return adjrows.at(i)->at(a) < adjrows.at(i)->at(b); } );
                top[i] = sets[i] + len_edge[i];
            }
            while (controlla_maggiore(sets, top)) {
                vector<unsigned int> labels(arguments.arg_num);
                for (int i = 0; i < arguments.arg_num; i++) {
                    labels[i] = adjrows.at(i)->at(vv[i][sets[i]]);
                }
                int massimo = max_elem(labels);
                if (massimo != -1) {
                    //sets[indice]++;
                    for (int i = 0; i < arguments.arg_num; i++) {
                        if (labels[i] != massimo) {
                            sets[i]++;
                        }
                    }
                }
                else {
                    //vector<int> min_sets = sets;
                    int min_sets[MAX_ARGS];
                    for (int  i=0; i<arguments.arg_num; i++) {
                        min_sets[i] = sets[i];
                    }
                    for (int i = 0; i < arguments.arg_num; i++) {
                        do { sets[i]++; } while (sets[i] < top[i] && adjrows.at(i)->at(vv[i][sets[i]]) == labels[0]);
                    }
                    int dif_sets[MAX_ARGS];
                    //vector<int> dif_sets(arguments.arg_num);
                    transform(sets, sets+arguments.arg_num, min_sets, dif_sets, std::minus<int>());
                    new_d.push_back({ min_sets, dif_sets, true });
                }
            }
        }
        else if (accumulate(len_edge, len_edge+arguments.arg_num, 1, multiplies<int>{})) {
            new_d.push_back({ sets, len_edge, true });
        }
    }
    return new_d;
}

// returns the index of the smallest value in arr that is >w.
// Assumption: such a value exists
// Assumption: arr contains no duplicates
// Assumption: arr has no values==INT_MAX
int index_of_next_smallest(int *arr, int start_idx, int len, int w) {
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

void remove_vtx_from_domain(int *left, Bidomain& bd, int v, int idx)
{
    int i = 0;
    while(left[bd.sets[idx] + i] != v) i++;
    std::swap(left[bd.sets[idx]+i], left[bd.sets[idx]+bd.len[idx]-1]);
    bd.len[idx]--;
}

void remove_bidomain(vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size()-1];
    domains.pop_back();
}

void solve_first_graph(const unsigned depth, vector<Graph>& g,
    AtomicIncumbent& global_incumbent,
    vector<VtxPair>& my_incumbent,
    vector<VtxPair>& current, vector<Bidomain>& domains,
    array<vector<int>, MAX_ARGS>& vv, const unsigned int matching_size_goal,
    unsigned long long& my_thread_nodes,
    array<int, MAX_ARGS>& nodi_inseriti,
    array<int, MAX_ARGS> order, Bidomain& bd)
{

    int pos = order[0];
    // riduciamo di 1 la dimensione di ciascun dominio non del primo grafo (assumendo di star prendendo un nodo)
    for (int i = 1; i < arguments.arg_num; i++) {
        bd.len[order[i]]--;
    }

    // scegliamo v come il minimo valore nel dominio
    int v = find_min_value(vv[pos].data(), bd.sets[pos], bd.len[pos]);
    // rimuoviamo v dal dominio, da ora in poi o si trova nella soluzione o no
    remove_vtx_from_domain(vv[pos].data(), bd, v, pos);

    // inseriamo v nella soluzione
    nodi_inseriti[pos] = v;
    
    // abbiamo un nodo, esploriamo
}

bool solve_other_graphs(const unsigned depth, vector<Graph>& g,
    AtomicIncumbent& global_incumbent,
    vector<VtxPair>& my_incumbent,
    vector<VtxPair>& current, vector<Bidomain>& domains,
    array<vector<int>, MAX_ARGS>& vv, const unsigned int matching_size_goal,
    unsigned long long& my_thread_nodes,
    int pos, Bidomain& bd, const int bd_idx,
    int& w)
{

    // prendiamo il più piccolo w del dominio
    int idx = index_of_next_smallest(vv[pos].data(), bd.sets[pos], bd.len[pos] + 1, w);
    if (idx == -1) {
        return false; // abbiamo finito i w del dominio
    }
    w = vv[pos][bd.sets[pos] + idx];

    // mettiamo w al fondo del dominio
    vv[pos][bd.sets[pos] + idx] = vv[pos][bd.sets[pos] + bd.len[pos]];
    vv[pos][bd.sets[pos] + bd.len[pos]] = w;

    return true; // abbiamo un nuovo w, esploriamo
}

void sorted_solve_nopar(const unsigned depth, vector<Graph> & g,
    AtomicIncumbent & global_incumbent,
    vector<VtxPair> & my_incumbent,
    vector<VtxPair> & current, vector<Bidomain> & domains,
    array<vector<int>, MAX_ARGS> &vv,
    const unsigned int matching_size_goal,
    unsigned long long & my_thread_nodes)
{
    if (arguments.verbose) {
        //show(current, domains, vv);
        string_show(current, domains, vv, depth);
    }

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

    
    // seleziona bidominio
    int bd_idx = select_bidomain(domains, vv[0].data(), current.size());
    if (bd_idx == -1)
        return; // abbiamo finito, più nulla da fare
    auto &bd = domains[bd_idx];

    array<int, MAX_ARGS> sorted_vv_idx;
    iota(sorted_vv_idx.begin(), sorted_vv_idx.begin() + MAX_ARGS, 0);
    //dobbiamo ordinare
    stable_sort(sorted_vv_idx.begin(), sorted_vv_idx.begin() + arguments.arg_num,
        [&](const int& a, const int& b) {
            return (bd.len[a] < bd.len[b]);
        }
    );

    array<int, MAX_ARGS> soluzione;
    for (int i = 0; i < MAX_ARGS; i++) {soluzione[i] = -1;}

    solve_first_graph(depth, g, global_incumbent, my_incumbent, current, domains, vv, matching_size_goal, my_thread_nodes, soluzione, sorted_vv_idx, bd);

    for (int i = 1; i > 0; ) {
        if (solve_other_graphs(depth, g, global_incumbent, my_incumbent, current, domains, vv, matching_size_goal, my_thread_nodes, sorted_vv_idx[i], bd, bd_idx, soluzione[sorted_vv_idx[i]]))
        {
            i ++;
            if (i == arguments.arg_num) {
                current.push_back(VtxPair(soluzione.data()));
                auto new_domains = filter_domains(domains, vv, g, soluzione.data(), arguments.directed || arguments.edge_labelled);
                sorted_solve_nopar(depth + 1, g, global_incumbent, my_incumbent, current, new_domains, vv, matching_size_goal, my_thread_nodes);
                i --;
                current.pop_back();
            }
        }
        else
        {
            soluzione[sorted_vv_idx[i]] = -1;
            i --;
        }
    }

    if (bd.len[sorted_vv_idx[0]] == 0)
    {
        remove_bidomain(domains, bd_idx);
    }
    else 
    {
        for (int i = 1; i < arguments.arg_num; i++) {
            bd.len[sorted_vv_idx[i]] ++;
        }
    }
    
    // adesso proviamo a proseguire senza prendere questo nodo (v)
    sorted_solve_nopar(depth + 1, g, global_incumbent, my_incumbent, current, domains, vv, matching_size_goal, my_thread_nodes);
}

void sorted_solve(const unsigned depth, vector<Graph>& g,
    AtomicIncumbent& global_incumbent,
    PerThreadIncumbents& per_thread_incumbents,
    vector<VtxPair>& current, vector<Bidomain>& domains,
    array<vector<int>, MAX_ARGS>& vv,
    const unsigned int matching_size_goal,
    const Position& position, HelpMe& help_me,
    unsigned long long& my_thread_nodes)
{


    if (arguments.verbose) {
        //show(current, domains, vv);
        string_show(current, domains, vv, depth);
    }
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

    int bd_idx = select_bidomain(domains, vv[0].data(), current.size());
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Bidomain& bd = domains[bd_idx];

    array<int, MAX_ARGS> sorted_vv_idx;
    iota(sorted_vv_idx.begin(), sorted_vv_idx.begin() + MAX_ARGS, 0);
    //dobbiamo ordinare
    stable_sort(sorted_vv_idx.begin(), sorted_vv_idx.begin() + arguments.arg_num,
        [&](const int& a, const int& b) {
            return (bd.len[a] < bd.len[b]);
        }
    );

    array<int, MAX_ARGS> soluzione = {};
    for (int i = 0; i < MAX_ARGS; i++) { soluzione[i] = -1; }

    std::atomic<int> shared_i{ 0 };
    const int i_end = bd.len[1] + 2; /* including the null */

    // Version of the loop used by helpers
    std::function<void(unsigned long long&)> helper_function =
        [&shared_i, &g, &global_incumbent, &per_thread_incumbents, &position,
        &depth, i_end, matching_size_goal, &help_me, current, domains, vv, bd_idx]
        (unsigned long long& help_thread_nodes) 
    {
        vector<VtxPair> help_current = current;
        vector<Bidomain> help_domains = domains;
        array<vector<int>, MAX_ARGS> help_vv = vv;
        array<int, MAX_ARGS> help_soluzione = {};
        for (int i = 0; i < MAX_ARGS; i++) { help_soluzione[i] = -1; }

        int w0_index = 0; //mi chiedo se la prima w è già stata esplorata da qualcun altro

        int which_i_should_i_run_next = shared_i++;

        if (which_i_should_i_run_next >= i_end)
            return; /* don't waste time recomputing */

        /* recalculate to this point */
        /* rerun important stuff from before the loop */
        int bd_idx = select_bidomain(help_domains, help_vv[0].data(), help_current.size());
        if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
            return;
        Bidomain& bd = help_domains[bd_idx];

        array<int, MAX_ARGS> sorted_vv_idx;
        iota(sorted_vv_idx.begin(), sorted_vv_idx.begin() + MAX_ARGS, 0);
        //dobbiamo ordinare
        stable_sort(sorted_vv_idx.begin(), sorted_vv_idx.begin() + arguments.arg_num,
            [&](const int& a, const int& b) {
            return (bd.len[a] < bd.len[b]);
        }
        );

        //gestiamo il primo grafo
        solve_first_graph(depth, g, global_incumbent, per_thread_incumbents.find(this_thread::get_id())->second, help_current, help_domains, help_vv, matching_size_goal, help_thread_nodes, help_soluzione, sorted_vv_idx, bd);

        //gestiamo il secondo grafo
        for (int i = 1; i > 0; ) {
            if (solve_other_graphs(depth, g, global_incumbent, per_thread_incumbents.find(this_thread::get_id())->second, help_current, help_domains, help_vv, matching_size_goal, help_thread_nodes, sorted_vv_idx[i], bd, bd_idx, help_soluzione[sorted_vv_idx[i]]))
            {
                bool should_i = (which_i_should_i_run_next == w0_index);
                w0_index += (i == 1);
                if (i != 1 || should_i) {
                    if (i == 1) {
                        which_i_should_i_run_next = shared_i++;
                    }
                    i++;
                    if (i == arguments.arg_num) {
                        help_current.push_back(VtxPair(help_soluzione.data()));
                        auto new_domains = filter_domains(help_domains, help_vv, g, help_soluzione.data(), arguments.directed || arguments.edge_labelled);
                        if (depth > split_levels) {
                            sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(this_thread::get_id())->second, help_current, new_domains, help_vv, matching_size_goal, help_thread_nodes);
                        }
                        else {
                            auto new_position = position;
                            new_position.add(depth, ++global_position);
                            sorted_solve(depth + 1, g, global_incumbent, per_thread_incumbents, help_current, new_domains, help_vv, matching_size_goal, new_position, help_me, help_thread_nodes);
                        }
                        i--;
                        help_current.pop_back();
                    }
                }
            }
            else
            {
                help_soluzione[sorted_vv_idx[i]] = -1;
                i--;
            }
        }

        if (bd.len[sorted_vv_idx[0]] == 0)
        {
            remove_bidomain(help_domains, bd_idx);
        }
        else
        {
            for (int i = 1; i < arguments.arg_num; i++) {
                bd.len[sorted_vv_idx[i]] ++;
            }
        }

        // adesso proviamo a proseguire senza prendere questo nodo (v)
        if (which_i_should_i_run_next == w0_index) {
            //which_i_should_i_run_next = shared_i++;
            if (depth > split_levels) {
                sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(this_thread::get_id())->second, help_current, help_domains, help_vv, matching_size_goal, help_thread_nodes);
            }
            else {
                auto new_position = position;
                new_position.add(depth, ++global_position);
                sorted_solve(depth + 1, g, global_incumbent, per_thread_incumbents, help_current, help_domains, help_vv, matching_size_goal, new_position, help_me, help_thread_nodes);
            }
        }

    };

    // Grab this first, before advertising that we can get help
    int which_i_should_i_run_next = shared_i++;

    // Version of the loop used by the main thread
    std::function<void(unsigned long long&)> main_function = [&](unsigned long long& main_thread_nodes) {

        solve_first_graph(depth, g, global_incumbent, per_thread_incumbents.find(this_thread::get_id())->second, current, domains, vv, matching_size_goal, my_thread_nodes, soluzione, sorted_vv_idx, bd);

        int w0_index = 0;

        for (int i = 1; i > 0; ) {
            if (solve_other_graphs(depth, g, global_incumbent, per_thread_incumbents.find(this_thread::get_id())->second, current, domains, vv, matching_size_goal, my_thread_nodes, sorted_vv_idx[i], bd, bd_idx, soluzione[sorted_vv_idx[i]]))
            {
                bool should_i = (which_i_should_i_run_next == w0_index);
                w0_index += (i == 1);
                if (i != 1 || should_i) {
                    if (i == 1) {
                        which_i_should_i_run_next = shared_i++;
                    }
                    i++;
                    if (i == arguments.arg_num) {
                        current.push_back(VtxPair(soluzione.data()));
                        auto new_domains = filter_domains(domains, vv, g, soluzione.data(), arguments.directed || arguments.edge_labelled);
                        if (depth > split_levels) {
                            sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(this_thread::get_id())->second, current, new_domains, vv, matching_size_goal, my_thread_nodes);
                        }
                        else {
                            auto new_position = position;
                            new_position.add(depth, ++global_position);
                            sorted_solve(depth + 1, g, global_incumbent, per_thread_incumbents, current, new_domains, vv, matching_size_goal, new_position, help_me, my_thread_nodes);
                        }
                        i--;
                        current.pop_back();
                    }
                }
            }
            else
            {
                soluzione[sorted_vv_idx[i]] = -1;
                i--;
            }
        }

        if (bd.len[sorted_vv_idx[0]] == 0)
        {
            remove_bidomain(domains, bd_idx);
        }
        else
        {
            for (int i = 1; i < arguments.arg_num; i++) {
                bd.len[sorted_vv_idx[i]] ++;
            }
        }

        // adesso proviamo a proseguire senza prendere questo nodo (v)
        if (which_i_should_i_run_next == w0_index) {
            //which_i_should_i_run_next = shared_i++;
            if (depth > split_levels) {
                sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(this_thread::get_id())->second, current, domains, vv, matching_size_goal, my_thread_nodes);
            }
            else {
                auto new_position = position;
                new_position.add(depth, ++global_position);
                sorted_solve(depth + 1, g, global_incumbent, per_thread_incumbents, current, domains, vv, matching_size_goal, new_position, help_me, my_thread_nodes);
            }
        }
    };

    if (depth <= split_levels) {
        help_me.get_help_with(position, main_function, helper_function, my_thread_nodes);
    }
    else {
        main_function(my_thread_nodes);
    }

}

// 1
void solve_nopar(const unsigned depth, vector<Graph>& g,
    AtomicIncumbent& global_incumbent,
    vector<VtxPair>& my_incumbent,
    vector<VtxPair>& current, vector<Bidomain>& domains,
    array<vector<int>, MAX_ARGS>& vv, const unsigned int matching_size_goal,
    unsigned long long& my_thread_nodes);
// 2
void solve(const unsigned depth, vector<Graph>& g,
    AtomicIncumbent& global_incumbent,
    PerThreadIncumbents& per_thread_incumbents,
    vector<VtxPair>& current, vector<Bidomain>& domains,
           array<vector<int>, MAX_ARGS>& vv, const unsigned int matching_size_goal,
    const Position& position, HelpMe& help_me, unsigned long long& my_thread_nodes);

// 3
void solve_nopar_recursive(Bidomain& bd, array<vector<int>, MAX_ARGS>& vv, std::vector<Bidomain>& domains,
    std::vector<Graph>& g, std::vector<VtxPair>& current, const unsigned int& depth,
    AtomicIncumbent& global_incumbent, std::vector<VtxPair>& my_incumbent,
    const unsigned int& matching_size_goal, unsigned long long& my_thread_nodes,
    int bd_idx, int *nodi_inseriti, int n_nodi_inseriti)
{
    if (abort_due_to_timeout)
        return;
    int w = -1;
    const int i_end = bd.len[n_nodi_inseriti] + 2; /* including the null */ //modificare

    for (int i = 0; i < i_end /* not != */; i++) {
        if (i != i_end - 1) {
            int idx = index_of_next_smallest(vv[n_nodi_inseriti].data(), bd.sets[n_nodi_inseriti], bd.len[n_nodi_inseriti] + 1, w);
            w = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx];
            //remove_vtx_from_domain(vv[n_nodi_inseriti], domains[bd_idx], w, n_nodi_inseriti);

            // swap w to the end of its colour class
            vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx] = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]];
            vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]] = w;

            nodi_inseriti[n_nodi_inseriti] = w;

            if (n_nodi_inseriti == arguments.arg_num - 1) {
                auto new_domains = filter_domains(domains, vv, g, nodi_inseriti, arguments.directed || arguments.edge_labelled);
                //auto new_domains = filter_domains(domains, vv[n_nodi_inseriti-1], vv[n_nodi_inseriti], g[n_nodi_inseriti-1], g[n_nodi_inseriti], nodi_inseriti.at(n_nodi_inseriti - 1), nodi_inseriti.at(n_nodi_inseriti), arguments.directed || arguments.edge_labelled, n_nodi_inseriti-1);
                current.push_back(VtxPair(nodi_inseriti));
                solve_nopar(depth + 1, g, global_incumbent, my_incumbent, current, new_domains, vv, matching_size_goal, my_thread_nodes);
                current.pop_back();
            }
            else {
                auto new_domains = domains;
                auto &new_bd = new_domains[bd_idx];
                solve_nopar_recursive(new_bd, vv, new_domains, g, current, depth, global_incumbent, my_incumbent, matching_size_goal, my_thread_nodes, bd_idx, nodi_inseriti, n_nodi_inseriti + 1);
            }
            //nodi_inseriti.pop_back();
        }
        else {
            // Last assign is null. Keep it in the loop to simplify parallelism.
            if (n_nodi_inseriti == 1) {
                //Bidomain tmp = bd;
                transform(bd.len + 1, bd.len+arguments.arg_num, bd.len + 1, [](int x) { return x + 1; });
                //bd.len[n_nodi_inseriti]++;
                if (bd.len[0] == 0)
                    remove_bidomain(domains, bd_idx);

                solve_nopar(depth + 1, g, global_incumbent, my_incumbent, current, domains, vv, matching_size_goal, my_thread_nodes); //attenzione potrebbe dover lanciare una funzione diversa
            }
        }
    }
}

void solve_nopar(const unsigned depth, vector<Graph> & g, /*g0, g1*/
        AtomicIncumbent & global_incumbent,
        vector<VtxPair> & my_incumbent,
        vector<VtxPair> & current, vector<Bidomain> & domains,
        array<vector<int>, MAX_ARGS> &vv,
        const unsigned int matching_size_goal,
        unsigned long long & my_thread_nodes)
{
    if (arguments.verbose) {
        //show(current, domains, vv);
        string_show(current, domains, vv, depth);
    }

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

    int bd_idx = select_bidomain(domains, vv[0].data(), current.size());
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Bidomain &bd = domains[bd_idx];

    transform(bd.len + 1, bd.len+arguments.arg_num, bd.len + 1, [](int x) { return x - 1; });
    std::atomic<int> shared_i{ 0 };

    int v = find_min_value(vv[0].data(), bd.sets[0], bd.len[0]);
    remove_vtx_from_domain(vv[0].data(), domains[bd_idx], v, 0); //0 � inteso

    int nodi_inseriti[MAX_ARGS];
    nodi_inseriti[0] = v;

    solve_nopar_recursive(bd, vv, domains, g, current, depth, global_incumbent, my_incumbent, matching_size_goal, my_thread_nodes, bd_idx, nodi_inseriti, 1);
}

// 4
void help_solve_recursive(const int& i_end, array<vector<int>, MAX_ARGS>& help_vv, Bidomain& help_bd, int& which_i_should_i_run_next,
                          std::atomic_int& shared_i, std::vector<Bidomain>& help_domains, std::vector<Graph>& g,
                          std::vector<VtxPair>& help_current, const unsigned int& depth, AtomicIncumbent& global_incumbent,
                          PerThreadIncumbents& per_thread_incumbents, const unsigned int& matching_size_goal,
                          unsigned long long& help_thread_nodes, const Position& position, HelpMe& help_me,
                          int help_bd_idx, int* nodi_inseriti, int n_nodi_inseriti)
{
    if (abort_due_to_timeout)
        return;
    int help_w = -1;

    for (int i = 0; i < i_end /* not != */; i++) {
        if (i != i_end - 1) {
            int idx = index_of_next_smallest(help_vv[n_nodi_inseriti].data(), help_bd.sets[n_nodi_inseriti], help_bd.len[n_nodi_inseriti] + 1, help_w);
            help_w = help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + idx];
            
            // swap w to the end of its colour class
            help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + idx] = help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + help_bd.len[n_nodi_inseriti]];
            help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + help_bd.len[n_nodi_inseriti]] = help_w;

            nodi_inseriti[n_nodi_inseriti] = help_w;
            if (n_nodi_inseriti == arguments.arg_num - 1) {
                auto new_domains = filter_domains(help_domains, help_vv, g, nodi_inseriti, arguments.directed || arguments.edge_labelled);

                help_current.push_back(VtxPair(nodi_inseriti));
                if (depth > split_levels) {
                    //solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, new_domains, help_vv, matching_size_goal, help_thread_nodes);
                    sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, new_domains, help_vv, matching_size_goal, help_thread_nodes);
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
                auto& new_help_bd = new_domains[help_bd_idx];
                help_solve_recursive(i_end, help_vv, new_help_bd, which_i_should_i_run_next, shared_i, new_domains, g,
                                     help_current, depth, global_incumbent,
                                     per_thread_incumbents, matching_size_goal, help_thread_nodes, position, help_me,
                                     help_bd_idx, nodi_inseriti, n_nodi_inseriti + 1);
            }
        }
        else {
            // Last assign is null. Keep it in the loop to simplify parallelism.
            if (n_nodi_inseriti == 1) {
                transform(help_bd.len + 1, help_bd.len+arguments.arg_num, help_bd.len + 1, [](int x) { return x + 1; });
                if (help_bd.len[0] == 0)
                    remove_bidomain(help_domains, help_bd_idx);

                if (depth > split_levels) {
                    //solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, help_domains, help_vv, matching_size_goal, help_thread_nodes);
                    sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, help_domains, help_vv, matching_size_goal, help_thread_nodes);
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

// 5
void solve_recursive(const int& i_end, array<vector<int>, MAX_ARGS>& vv, Bidomain& bd, int which_i_should_i_run_next,
                     std::atomic_int& shared_i, std::vector<Bidomain>& domains, std::vector<Graph>& g,
                     std::vector<VtxPair>& current, const unsigned int& depth, AtomicIncumbent& global_incumbent,
                     PerThreadIncumbents& per_thread_incumbents, const unsigned int& matching_size_goal,
                     unsigned long long& main_thread_nodes, const Position& position, HelpMe& help_me,
                     int bd_idx, int* nodi_inseriti, int n_nodi_inseriti)
{
    if (abort_due_to_timeout)
        return;
    int w = -1;

    for (int i = 0; i < i_end /* not != */; i++) {
        if (i != i_end - 1) {
            int idx = index_of_next_smallest(vv[n_nodi_inseriti].data(), bd.sets[n_nodi_inseriti], bd.len[n_nodi_inseriti] + 1, w);
            w = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx];
            
            // swap w to the end of its colour class
            vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx] = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]];
            vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]] = w;

            nodi_inseriti[n_nodi_inseriti] = w;
            if (n_nodi_inseriti == arguments.arg_num - 1) {
                auto new_domains = filter_domains(domains, vv, g, nodi_inseriti, arguments.directed || arguments.edge_labelled);
                current.push_back(VtxPair(nodi_inseriti));
                if (depth > split_levels) {
                    //solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, new_domains, vv, matching_size_goal, main_thread_nodes);
                    sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, new_domains, vv, matching_size_goal, main_thread_nodes);
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
        }
        else {
            // Last assign is null. Keep it in the loop to simplify parallelism.
            if (n_nodi_inseriti == 1) {
                Bidomain tmp = domains[bd_idx];
                transform(bd.len + 1, bd.len+arguments.arg_num, bd.len + 1, [](int x) { return x + 1; });
                if (bd.len[0] == 0)
                    remove_bidomain(domains, bd_idx);


                if (depth > split_levels) {
                    //solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, domains, vv, matching_size_goal, main_thread_nodes);
                    sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, domains, vv, matching_size_goal, main_thread_nodes);
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

// 1
void solve(const unsigned depth, vector<Graph> & g,
        AtomicIncumbent & global_incumbent,
        PerThreadIncumbents & per_thread_incumbents,
        vector<VtxPair> & current, vector<Bidomain> & domains,
        array<vector<int>, MAX_ARGS> & vv,
        const unsigned int matching_size_goal,
        const Position & position, HelpMe & help_me, unsigned long long & my_thread_nodes)
{

    if (arguments.verbose) {
        //show(current, domains, vv);
        string_show(current, domains, vv, depth);
    }
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

    int bd_idx = select_bidomain(domains, vv[0].data(), current.size());
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Bidomain &bd = domains[bd_idx];

    transform(bd.len + 1, bd.len+arguments.arg_num, bd.len + 1, [](int x) { return x - 1; });
    std::atomic<int> shared_i{ 0 };
    const int i_end = bd.len[1] + 2; /* including the null */

    array<vector<int>, MAX_ARGS> helper_vv;
    for (int i=0; i<arguments.arg_num; i++) {
        helper_vv[i].reserve(vv.size());
        for (int j=0; j<vv[i].size(); j++) {
            helper_vv[i].push_back(vv[i][j]);
        }
    }

    // Version of the loop used by helpers
    std::function<void (unsigned long long &)> helper_function = [&shared_i, &g, &global_incumbent, &per_thread_incumbents, &position, &depth,
        i_end, matching_size_goal, current, domains, helper_vv, &help_me] (unsigned long long & help_thread_nodes) {
        int which_i_should_i_run_next = shared_i++;

        if (which_i_should_i_run_next >= i_end)
            return; /* don't waste time recomputing */

        /* recalculate to this point */
        vector<VtxPair> help_current = current;
        vector<Bidomain> help_domains = domains;
        array<vector<int>, MAX_ARGS> help_vv;
        for (int i=0; i<arguments.arg_num; i++) {
            help_vv[i].reserve(helper_vv.size());
            for (int j=0; j<helper_vv[i].size(); j++) {
                help_vv[i].push_back(helper_vv[i][j]);
            }
        }

        /* rerun important stuff from before the loop */
        int help_bd_idx = select_bidomain(help_domains, help_vv[0].data(), help_current.size());
        if (help_bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
            return;
        Bidomain &help_bd = help_domains[help_bd_idx];

        int help_v = find_min_value(help_vv[0].data(), help_bd.sets[0], help_bd.len[0]);
        remove_vtx_from_domain(help_vv[0].data(), help_domains[help_bd_idx], help_v, 0); //0 � inteso

        int nodi_inseriti[MAX_ARGS];
        nodi_inseriti[0] = help_v;
        
        int n_nodi_inseriti = 1;
        int help_w = -1;

        for (int i = 0; i < i_end /* not != */; i++) {
            if (i != i_end - 1) {
                int idx = index_of_next_smallest(help_vv[n_nodi_inseriti].data(), help_bd.sets[n_nodi_inseriti], help_bd.len[n_nodi_inseriti] + 1, help_w);
                help_w = help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + idx];
                
                // swap w to the end of its colour class
                help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + idx] = help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + help_bd.len[n_nodi_inseriti]];
                help_vv[n_nodi_inseriti][help_bd.sets[n_nodi_inseriti] + help_bd.len[n_nodi_inseriti]] = help_w;

                if (i == which_i_should_i_run_next) {
                    which_i_should_i_run_next = shared_i++;
                    nodi_inseriti[n_nodi_inseriti] = help_w;
                    if (n_nodi_inseriti == arguments.arg_num - 1) {
                        auto new_domains = filter_domains(help_domains, help_vv, g, nodi_inseriti, arguments.directed || arguments.edge_labelled);
                        
                        help_current.push_back(VtxPair(nodi_inseriti));
                        if (depth > split_levels) {
                            //solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, new_domains, help_vv, matching_size_goal, help_thread_nodes);
                            sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, new_domains, help_vv, matching_size_goal, help_thread_nodes);
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
                        help_solve_recursive(i_end, help_vv, new_help_bd, which_i_should_i_run_next, shared_i,
                                             new_domains, g, help_current, depth, global_incumbent,
                                             per_thread_incumbents, matching_size_goal, help_thread_nodes, position,
                                             help_me, help_bd_idx, nodi_inseriti, n_nodi_inseriti + 1);
                    }
                }
            }
            else {
                // Last assign is null. Keep it in the loop to simplify parallelism.
                if (n_nodi_inseriti == 1) {
                    transform(help_bd.len + 1, help_bd.len+arguments.arg_num, help_bd.len + 1, [](int x) { return x + 1; });
                    if (help_bd.len[0] == 0)
                        remove_bidomain(help_domains, help_bd_idx);

                    if (i == which_i_should_i_run_next) {
                        which_i_should_i_run_next = shared_i++;
                        if (depth > split_levels) {
                            //solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, help_domains, help_vv, matching_size_goal, help_thread_nodes);
                            sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, help_domains, help_vv, matching_size_goal, help_thread_nodes);
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
        int v = find_min_value(vv[0].data(), bd.sets[0], bd.len[0]);
        remove_vtx_from_domain(vv[0].data(), domains[bd_idx], v, 0); //0 � inteso

        int nodi_inseriti[MAX_ARGS];
        nodi_inseriti[0] = v;
        int n_nodi_inseriti = 1;

        int w = -1;

        for (int i = 0; i < i_end /* not != */; i++) {
            if (i != i_end - 1) {
                int idx = index_of_next_smallest(vv[n_nodi_inseriti].data(), bd.sets[n_nodi_inseriti], bd.len[n_nodi_inseriti] + 1, w);
                w = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx];
                //remove_vtx_from_domain(vv[n_nodi_inseriti], domains[bd_idx], w, n_nodi_inseriti);

                // swap w to the end of its colour class
                vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + idx] = vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]];
                vv[n_nodi_inseriti][bd.sets[n_nodi_inseriti] + bd.len[n_nodi_inseriti]] = w;

                if (i == which_i_should_i_run_next) {

                    nodi_inseriti[n_nodi_inseriti] = w;
                    which_i_should_i_run_next = shared_i++;
                    if (n_nodi_inseriti == arguments.arg_num - 1) {
                        auto new_domains = filter_domains(domains, vv, g, nodi_inseriti, arguments.directed || arguments.edge_labelled);
                        
                        current.push_back(VtxPair(nodi_inseriti));
                        if (depth > split_levels) {
                            //solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, new_domains, vv, matching_size_goal, main_thread_nodes);
                            sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, new_domains, vv, matching_size_goal, main_thread_nodes);
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
                }
            }
            else {
                // Last assign is null. Keep it in the loop to simplify parallelism.
                if (n_nodi_inseriti == 1) {
                    Bidomain tmp = domains[bd_idx];
                    transform(bd.len + 1, bd.len+arguments.arg_num, bd.len + 1, [](int x) { return x + 1; });
                    if (bd.len[0] == 0) //1 inteso
                        remove_bidomain(domains, bd_idx);

                    if (i == which_i_should_i_run_next) {
                        which_i_should_i_run_next = shared_i++;
                        if (depth > split_levels) {
                            //solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, domains, vv, matching_size_goal, main_thread_nodes);
                            sorted_solve_nopar(depth + 1, g, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, domains, vv, matching_size_goal, main_thread_nodes);
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

std::pair<vector<VtxPair>, unsigned long long> mcs(vector<Graph> & gi) {
    
    // the buffer of vertex indices for the partitions
    array<vector<int>, MAX_ARGS> vtx_buf;

    auto domains = vector<Bidomain>{};

    vector<set<unsigned int>> labels_vv (arguments.arg_num);
    for (int i = 0; i < arguments.arg_num; i++) {
        for (unsigned int label : gi[i].label) {
            labels_vv[i].insert(label);
        }
    }

    std::set<unsigned int> labels = intersection(labels_vv);
    
    // Create a bidomain for each label that appears in both graphs
    for (unsigned int label : labels) {
        int starts[MAX_ARGS];
        int len[MAX_ARGS];
        for (int i = 0; i < arguments.arg_num; i++) {
            int jj = 0;
            starts[i] = vtx_buf[i].size();
            for (int j = 0; j < gi[i].n; j++) {
                if (gi[i].label[j] == label) {
                    vtx_buf[i].push_back(j);
                }
            }
            len[i] = vtx_buf[i].size() - starts[i];
        }

        domains.push_back({ starts, len, false });
    }

    AtomicIncumbent global_incumbent;
    vector<VtxPair> incumbent;
    unsigned long long global_nodes = 0;

    if (arguments.big_first) {
        for (int k=0; k<gi[0].n; k++) {
            unsigned int goal = gi[0].n - k;
            array<vector<int>, MAX_ARGS> vtx_buf_copy;
            //vtx_buf_copy.resize(vtx_buf.size());
            for (int i=0; i<arguments.arg_num; i++) {
                vtx_buf_copy[i].reserve(vtx_buf[i].size());
                for (int j=0; j<vtx_buf[i].size(); j++) {
                    vtx_buf_copy[i].push_back(vtx_buf[i][j]);
                }
            }
            auto domains_copy = domains;
            vector<VtxPair> current;
            PerThreadIncumbents per_thread_incumbents;
            per_thread_incumbents.emplace(std::this_thread::get_id(), vector<VtxPair>());
            Position position;
            HelpMe help_me(arguments.threads - 1);
            for (auto & t : help_me.threads)
                per_thread_incumbents.emplace(t.get_id(), vector<VtxPair>());
#if UNSORTED
            solve(0, gi, global_incumbent, per_thread_incumbents, current, domains_copy, vtx_buf_copy, goal, position, help_me, global_nodes);
#else
            sorted_solve(0, gi, global_incumbent, per_thread_incumbents, current, domains_copy, vtx_buf_copy, goal, position, help_me, global_nodes);
#endif
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
        vector<VtxPair> current;
        PerThreadIncumbents per_thread_incumbents;
        per_thread_incumbents.emplace(std::this_thread::get_id(), vector<VtxPair>());
        Position position;
        HelpMe help_me(arguments.threads - 1);
        for (auto & t : help_me.threads)
            per_thread_incumbents.emplace(t.get_id(), vector<VtxPair>());
#if UNSORTED
        solve(0, gi, global_incumbent, per_thread_incumbents, current, domains, vtx_buf, 1, position, help_me, global_nodes);
#else
        sorted_solve(0, gi, global_incumbent, per_thread_incumbents, current, domains, vtx_buf, 1, position, help_me, global_nodes);
#endif
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

void mcsp::start(args &arg) {
    arguments = arg;

    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : arguments.bin_enrico ? 'E' : arguments.ioi ? 'I' : 'B';
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

    std::pair<vector<VtxPair>, unsigned long long> solution = mcs(gi_sorted);

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
        
    cout << ">>> " << solution.first.size() << " - " << solution.second << " - " << (double) time_elapsed/1000 << endl;
}

