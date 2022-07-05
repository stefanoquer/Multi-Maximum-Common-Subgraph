#include <iostream>
#include <thread>
#include <vector>

#include "graph.h"
#include "SecureQueue.h"

#if !USING_WINDOWS
#include <argp.h>
#endif
#include <limits.h>
#include <locale.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <signal.h>
//#include <io.h>

using std::chrono::steady_clock;
using std::chrono::milliseconds;
using std::chrono::duration_cast;

using namespace std;

#if USING_WINDOWS
using namespace Concurrency;
//#include <windows.h>
#include <ppl.h>
#include <amp.h>
#endif
#ifndef _WIN32
#include <unistd.h>
#endif // !_WIN32


thread t;
int GRAPH_MAX;
int MAX_N = 0;
int biggest_graph = 0;


typedef unsigned long long ULL;

void array_copy(int [], int [], int);

int timeout = 0;
//cancellation_token_source cts;


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

/*
 * Signal handler (for timeout)
 */

void sig_alarm(int sig) {
    fprintf(stderr, "Time out reached.\n");
    timeout = 1;

    return;
}

void swap(int *a, int *b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

static void fail(char *msg) {
    printf("%s\n", msg);
    exit(1);
}

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Find a maximum clique in a graph in DIMACS format";
static char args_doc[] = "FILENAME 1..n";
#if !USING_WINDOWS
static struct argp_option options[] = {
        {"quiet",     'q', 0,         0, "Quiet output"},
        {"verbose",   'v', 0,         0, "Verbose output"},
        {"dimacs",    'd', 0,         0, "Read DIMACS format"},
        {"lad",       'l', 0,         0, "Read LAD format"},
        {"connected", 'c', 0,         0, "Solve max common CONNECTED subgraph problem"},
        {"timeout",   't', "timeout", 0, "Specify a timeout (seconds)"},
        {0}
};
#endif

static struct {
    bool quiet;
    bool verbose;
    bool connected;
    bool dimacs;
    bool lad;
    vector<char *> filename;
    int timeout;
    int arg_num;
} arguments;

void set_default_arguments() {
    int i;
    arguments.quiet = false;
    arguments.verbose = false;
    arguments.connected = false;
    arguments.dimacs = false;
    arguments.lad = false;
    /*for (i = 0; i < GRAPH_MAX; i++)
        arguments.filename[i] = NULL;*/
    arguments.timeout = 0;
    arguments.arg_num = 0;
}
#if !USING_WINDOWS
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
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
            arguments.connected = true;
            break;
        case 't':
            arguments.timeout = atoi(arg);
            break;
        case ARGP_KEY_ARG:
            arguments.filename.push_back(arg);
            //argp_usage(state);
            arguments.arg_num++;
            break;
        case ARGP_KEY_END:
            if (arguments.arg_num == 0)
                argp_usage(state);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};
#endif
/*******************************************************************************
                                     Stats
*******************************************************************************/

static long nodes = 0;

/*******************************************************************************
                                 MCS functions
*******************************************************************************/

struct VtxPair {
    vector<int> v;
    VtxPair() {
        v.resize(GRAPH_MAX);
    }
};

struct VtxPairList {
    vector<struct VtxPair> vals;
    int len;
    VtxPairList() {
        vals.resize(MAX_N);
    }
};

struct Bidomain {
    vector<vector<int>> vv;
    vector<int> len;
    bool is_adjacent;
    Bidomain() {
        vv.resize(GRAPH_MAX, vector<int>(MAX_N));
        len.resize(GRAPH_MAX);
    }
};

struct BidomainList {
    vector<struct Bidomain> vals;
    int len = 0;
    BidomainList() {
        vals.resize(MAX_N * 2);
    }
};

struct D {
    vector<struct Graph *> g;
    struct VtxPairList *incumbent;
    struct VtxPairList *current;
    struct BidomainList *domains;
    D() {
        g.resize(GRAPH_MAX);
    }
};

struct dati_Solve {
    struct D d;
    struct D new_d;
    struct Bidomain *bd;
    int level;
    vector<int> vertex = vector<int> (GRAPH_MAX+1);
    vector<struct BidomainList> preallocated_lists = vector<struct BidomainList> (biggest_graph+1);
    VtxPairList incombenti;
    VtxPairList correnti;
};

struct AtomicIncumbent {
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
                if (value.compare_exchange_strong(cur_v, v, std::memory_order_seq_cst))
                    return true;
            }
            else
                return false;
        }
    }
} largest_incumbent;

int calc_bound(struct BidomainList *domains) {

    int bound = 0;
    for (int i = 0; i < domains->len; i++) {
        struct Bidomain *bd = &domains->vals[i];
        int min_size = INT_MAX;
        for (int j = 0; j < arguments.arg_num; j++) {
            if (bd->len[j] < min_size)
                min_size = bd->len[j];
        }
        bound += min_size;
    }
    return bound;
}

struct Bidomain *select_bidomain(struct BidomainList *domains, int current_matching_size) {
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    struct Bidomain *best = NULL;
    for (int i = 0; i < domains->len; i++) {
        struct Bidomain *bd = &domains->vals[i];
        if (arguments.connected && current_matching_size > 0 && !bd->is_adjacent) continue;
        int len = bd->len[0];
        for (int j = 1; j < arguments.arg_num; j++) {
            if (len < bd->len[j])
                len = bd->len[j];
        }
        if (len < min_size) {
            min_size = len;
            best = bd;
        }
    }
    return best;
}

int select_bidomain_position(struct BidomainList* domains, int current_matching_size) {
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int best = -1;
    for (int i = 0; i < domains->len; i++) {
        struct Bidomain* bd = &domains->vals[i];
        if (arguments.connected && current_matching_size > 0 && !bd->is_adjacent) continue;
        int len = bd->len[0];
        for (int j = 1; j < arguments.arg_num; j++) {
            if (len < bd->len[j])
                len = bd->len[j];
        }
        if (len < min_size) {
            min_size = len;
            best = i;
        }
    }
    return best;
}

// Returns length of left half of array
// DA NON CAMBIARE
int partition(int *vv, int vv_len, unsigned char *adjrow) {
    int i = 0;
    for (int j = 0; j < vv_len; j++) {
        int adj = adjrow[vv[j]];
        if (adj) {
            swap(&vv[i], &vv[j]);
            i++;
        }
    }
    return i;
}

void add_bidomain(struct BidomainList* bd_list, vector<vector<int>> &vvMat, int* lenV, bool is_adjacent) {

    int len = bd_list->len++;
    for (int ng = 0; ng < arguments.arg_num; ng++) {
        bd_list->vals[len].vv[ng] = vvMat[ng];
        bd_list->vals[len].len[ng] = lenV[ng];
    }
    bd_list->vals[len].is_adjacent = is_adjacent;
}

void add_bidomain(
        struct BidomainList *bd_list, int **vvMat, int *lenV, bool is_adjacent
) {
    for (int ng = 0; ng < arguments.arg_num; ng++) {
        bd_list->vals[bd_list->len].vv[ng] = vector<int>(vvMat[ng], vvMat[ng] + lenV[ng]);
        bd_list->vals[bd_list->len].len[ng] = lenV[ng];
    }
    bd_list->vals[bd_list->len].is_adjacent = is_adjacent;
    bd_list->len++;
}


void filter_domains(
        struct BidomainList *d, struct BidomainList *new_d,
        struct Graph *g[], int vertex[]
) {
    //parallel_filter_domains(d, new_d, g, vertex);
    //return;

    int flag_len_edge, flag_len_noedge;
    vector<int> len_edge(GRAPH_MAX), len_noedge(GRAPH_MAX);

    new_d->len = 0;
    for (int j = 0; j < d->len; j++) {
        struct Bidomain *old_bd = &d->vals[j];

        flag_len_edge = flag_len_noedge = 0;
        for (int i = 0; i < arguments.arg_num; i++) {
            len_edge[i] = partition(&old_bd->vv[i][0], old_bd->len[i], &g[i]->adjmat[vertex[i]][0]);
            len_noedge[i] = old_bd->len[i] - len_edge[i];

            if (len_edge[i])
                flag_len_edge++;
            if (len_noedge[i])
                flag_len_noedge++;
        }

        // Se ci sono vertici NON adiancenti su g e h aggiungo il bidomain
        if (flag_len_noedge == arguments.arg_num) {
            vector<int*>vv(GRAPH_MAX);
            for (int i = 0; i < arguments.arg_num; i++) {
                //vv[i] = old_bd->vv[i] + len_edge[i];
                vv[i] = &old_bd->vv[i][0] + len_edge[i];
            }
            add_bidomain(new_d, &vv[0], &len_noedge[0], old_bd->is_adjacent);
        }

        // Se ci sono vertici adiacenti su g e h aggiungo il bidomain
        if (flag_len_edge == arguments.arg_num) {
            add_bidomain(new_d, old_bd->vv, &len_edge[0], true);
        }
    }
}

void remove_bidomain(struct BidomainList *list, struct Bidomain *b) {
    int i = b - &list->vals[0];
    list->vals[i] = list->vals[list->len - 1];
    list->len--;
}

void show(struct D &d) {
    printf("Nodes: %ld\n", nodes);
    printf("Length of current assignment: %d\n", d.current->len);
    printf("Current assignment:");

    for (int i = 0; i < d.current->len; i++) {
        for (int j = 0; j < arguments.arg_num; j++) {
            if (j == 0)
                printf("  %d", d.current->vals[i].v[j]);
            else
                printf("->%d", d.current->vals[i].v[j]);
        }
    }
    printf("\n");
    for (int i = 0; i < d.domains->len; i++) {
        struct Bidomain bd = d.domains->vals[i];
        for (int ng = 0; ng < arguments.arg_num; ng++) {
            printf("Graph %d  ", ng);
            for (int j = 0; j < bd.len[ng]; j++)
                printf("%d ", bd.vv[ng][j]);
            printf("\n");
        }
    }
    printf("\n\n");
}

void set_incumbent(struct VtxPairList *current, struct VtxPairList *incumbent) {
    incumbent->len = current->len;
    for (int i = 0; i < current->len; i++)
        incumbent->vals[i] = current->vals[i];
}

int find_and_remove_min_value(vector<int> &arr, int &len) {
    int min_v = INT_MAX;
    int idx = -1;
    for (int i = 0; i < len; i++) {
        if (arr.at(i) < min_v) {
            min_v = arr.at(i);
            idx = i;
        }
    }
    swap(arr.at(idx), arr.at(len - 1));
    (len)--;
    return min_v;
}

// returns the index of the smallest value in arr that is >w.
// Assumption: such a value exists
// Assumption: arr contains no duplicates
// Assumption: arr has no values==INT_MAX
int index_of_next_smallest(vector<int> &arr, int len, int w) {
    int idx = -1;
    int smallest = INT_MAX;
    for (int i = 0; i < len; i++) {
        if (arr.at(i) > w && arr.at(i) < smallest) {
            smallest = arr.at(i);
            idx = i;
        }
    }
    return idx;
}


void solve(struct D &d, struct D &new_d, struct Bidomain *bd, int vertex[], int level, struct BidomainList *preallocated_lists) {
    int ng, r, c;

    if (timeout == 1/* || cts.get_token().is_canceled()*/)
        return;

    r = (level - 1) % arguments.arg_num;
    c = (level - 1) / arguments.arg_num;

    if (r == 0) {
        if (arguments.verbose) show(d);
        nodes++;

        /*if (d.domains->len == 2 && d.domains->vals[0].len[0] == 2 && d.domains->vals[0].vv[0][0] == 8 && d.domains->vals[0].vv[0][1] == 7) {
            cout << "controllo" << endl;
        }*/

        //if (d.current->len > d.incumbent->len) {
        if (largest_incumbent.update(d.current->len)) {
            set_incumbent(d.current, d.incumbent);
            if (!arguments.quiet) {
                printf("Incumbent size: %d\n", d.incumbent->len);
            }
        }

        //if (d.current->len + calc_bound(d.domains) <= d.incumbent->len)
        if (d.current->len + calc_bound(d.domains) <= largest_incumbent.value)
            return;

        struct Bidomain *bd = select_bidomain(d.domains, d.current->len);
        if (bd == NULL)   // In the MCCS case, there may be nothing we can branch on
            return;
        struct Bidomain bd_copy;

        vertex[0] = find_and_remove_min_value(bd->vv[0], bd->len[0]);

        if (bd->len[0] == 0) {
            bd_copy = *bd;
            remove_bidomain(d.domains, bd);
            bd = &bd_copy;
        }

        struct D new_d = d;
        new_d.domains = &preallocated_lists[c + 1];

        solve(d, new_d, bd, vertex, level + 1, preallocated_lists);
        solve(d, new_d, bd, vertex, level + arguments.arg_num, preallocated_lists);
#if DEBUG
        printf ("0 >>> D AFTER FILTER DOMAINS\n"); show(d);
        printf ("0 >>> NEW_D AFTER FILTER DOMAINS\n"); show(new_d);
#endif
    }

    if (r > 0) {
        // Try assigning v to each vertex w in bd->right_vv, in turn
        bd->len[r]--;
        vertex[r] = -1;
        for (int i = 0; i <= bd->len[r]; i++) {
            /*if (r == arguments.arg_num - 1) {
                cout << "controllo" << endl;
            }*/
            int idx = index_of_next_smallest(bd->vv[r], bd->len[r] + 1, vertex[r]);
            vertex[r] = bd->vv[r][idx];

            // swap w with the value just past the end of the right_vv array
            bd->vv[r][idx] = bd->vv[r][bd->len[r]];
            bd->vv[r][bd->len[r]] = vertex[r];

            if (r == arguments.arg_num - 1) {
                filter_domains(d.domains, new_d.domains, &d.g[0], vertex);
                for (ng = 0; ng < arguments.arg_num; ng++) {
                    d.current->vals[d.current->len].v[ng] = vertex[ng];
                }
                d.current->len++;

                vector<int> new_vertex (GRAPH_MAX+1);
                array_copy(&new_vertex[0], vertex, arguments.arg_num);
                solve(new_d, new_d, bd, &new_vertex[0], level + 1, preallocated_lists);

                d.current->len--;
            } else {
                solve(d, new_d, bd, vertex, level + 1, preallocated_lists);
            }

#if DEBUG
            printf ("1 >>> D AFTER FILTER DOMAINS\n"); show(d);
            printf ("1 >>> NEW_D AFTER FILTER DOMAINS\n"); show(*dP);
#endif
        }
        bd->len[r]++;

        // Qui occorrerebbe ricorrere evitando di fare match su un nodo del grafo
        // corrente non solo saltare il nodo del primo grafo
        //solve(d, &d, &bd, vertex, level+1);
#if DEBUG
        printf ("2 >>> D AFTER FILTER DOMAINS\n"); show(d);
        printf ("2 >>> NEW_D AFTER FILTER DOMAINS\n"); show(*dP);
#endif
    }
}
vector<dati_Solve *> coda_dati;
SecureQueue<int> sq;
#if USING_WINDOWS
vector<task<void>> tasks;
#endif

void win_parallel_solve(int indice) {

    dati_Solve* data = coda_dati.at(indice);

    solve(data->d, data->new_d, data->bd, &data->vertex[0], data->level, &data->preallocated_lists[0]);

    return;

    /*struct D d = data->d;
    struct D new_d = data->new_d;
    struct Bidomain* bd = &data->bd;
    int* vertex = &data->vertex[0];
    int level = data->level;
    struct BidomainList* preallocated_lists = &data->preallocated_lists[0];

    int ng, r, c;

    if (timeout == 1)
        return;

    r = (level - 1) % arguments.arg_num;
    c = (level - 1) / arguments.arg_num;

    if (r == 0) {
        if (arguments.verbose) show(d);
        nodes++;

        //if (d.current->len > d.incumbent->len) {
        if (largest_incumbent.update(d.current->len)) {
            set_incumbent(d.current, d.incumbent);
            if (!arguments.quiet) printf("Incumbent size: %d\n", d.incumbent->len);
        }

        //if (d.current->len + calc_bound(d.domains) <= d.incumbent->len)
        if (d.current->len + calc_bound(d.domains) <= largest_incumbent.value)
            return;

        struct Bidomain* bd = select_bidomain(d.domains, d.current->len);
        if (bd == NULL)   // In the MCCS case, there may be nothing we can branch on
            return;
        struct Bidomain bd_copy;

        vertex[0] = find_and_remove_min_value(bd->vv[0], bd->len[0]);
        if (bd->len[0] == 0) {
            bd_copy = *bd;
            remove_bidomain(d.domains, bd);
            bd = &bd_copy;
        }

        struct D new_d = d;
        new_d.domains = &preallocated_lists[c + 1];

        solve(d, new_d, bd, vertex, level + 1, preallocated_lists);
        solve(d, new_d, bd, vertex, level + arguments.arg_num, preallocated_lists);
#if DEBUG
        printf("0 >>> D AFTER FILTER DOMAINS\n"); show(d);
        printf("0 >>> NEW_D AFTER FILTER DOMAINS\n"); show(new_d);
#endif
    }

    if (r > 0) {
        // Try assigning v to each vertex w in bd->right_vv, in turn
        bd->len[r]--;
        vertex[r] = -1;
        for (int i = 0; i <= bd->len[r]; i++) {
            int idx = index_of_next_smallest(bd->vv[r], bd->len[r] + 1, vertex[r]);
            vertex[r] = bd->vv[r][idx];

            // swap w with the value just past the end of the right_vv array
            bd->vv[r][idx] = bd->vv[r][bd->len[r]];
            bd->vv[r][bd->len[r]] = vertex[r];

            if (r == arguments.arg_num - 1) {
                filter_domains(d.domains, new_d.domains, &d.g[0], vertex);
                for (ng = 0; ng < arguments.arg_num; ng++) {
                    d.current->vals[d.current->len].v[ng] = vertex[ng];
                }
                d.current->len++;

                vector<int> new_vertex(GRAPH_MAX + 1);
                array_copy(&new_vertex[0], vertex, arguments.arg_num);
                solve(new_d, new_d, bd, &new_vertex[0], level + 1, preallocated_lists);

                d.current->len--;
            }
            else {
                solve(d, new_d, bd, vertex, level + 1, preallocated_lists);
            }

#if DEBUG
            printf("1 >>> D AFTER FILTER DOMAINS\n"); show(d);
            printf("1 >>> NEW_D AFTER FILTER DOMAINS\n"); show(*dP);
#endif
        }
        bd->len[r]++;

        // Qui occorrerebbe ricorrere evitando di fare match su un nodo del grafo
        // corrente non solo saltare il nodo del primo grafo
        //solve(d, &d, &bd, vertex, level+1);
#if DEBUG
        printf("2 >>> D AFTER FILTER DOMAINS\n"); show(d);
        printf("2 >>> NEW_D AFTER FILTER DOMAINS\n"); show(*dP);
#endif
    }*/


}

void parallel_solve() {
    int indice;
    while ((indice = sq.pop()) != -1) {
        win_parallel_solve(indice);
    }
}

void crea_dato(const struct D &d, const struct D &new_d, const struct Bidomain &bd_noRef, int vertex[], int level, int d_domain, int new_d_domain, vector<struct BidomainList> preallocated_lists, int b_pos) {
    coda_dati.push_back(new dati_Solve());
    dati_Solve *dati = coda_dati.back();

    dati->preallocated_lists.at(d_domain).len = d.domains->len;
    for (int j = 0; j < d.domains->vals.size(); j++) {
        for (int k = 0; k < min(GRAPH_MAX, (int)d.domains->vals[j].len.size()); k++) {
            dati->preallocated_lists[d_domain].vals[j].len[k] = d.domains->vals[j].len[k];
        }
        dati->preallocated_lists[d_domain].vals[j].is_adjacent = d.domains->vals[j].is_adjacent;
        for (int k = 0; k < d.domains->vals[j].vv.size(); k++) {
            dati->preallocated_lists[d_domain].vals[j].vv[k] = d.domains->vals[j].vv[k];
        }
    }
    dati->preallocated_lists.at(new_d_domain).len = new_d.domains->len;
    for (int j = 0; j < new_d.domains->vals.size(); j++) {
        for (int k = 0; k < min(GRAPH_MAX, (int)new_d.domains->vals[j].len.size()); k++) {
            dati->preallocated_lists[new_d_domain].vals[j].len[k] = new_d.domains->vals[j].len[k];
        }
        dati->preallocated_lists[new_d_domain].vals[j].is_adjacent = new_d.domains->vals[j].is_adjacent;
        for (int k = 0; k < new_d.domains->vals[j].vv.size(); k++) {
            dati->preallocated_lists[new_d_domain].vals[j].vv[k] = new_d.domains->vals[j].vv[k];
        }
    }
    dati->incombenti.len = d.incumbent->len;
    for(int i=0; i<MAX_N; i++) {
        for(int j=0; j<d.incumbent->vals[i].v.size(); j++) {
            dati->incombenti.vals[i].v[j] = d.incumbent->vals[i].v[j];
        }
    }
    dati->correnti.len = d.current->len;
    for(int i=0; i<MAX_N; i++) {
        for(int j=0; j<d.current->vals[i].v.size(); j++) {
            dati->correnti.vals[i].v[j] = d.current->vals[i].v[j];
        }
    }



    for(int i=0; i<d.g.size(); i++) {
        dati->new_d.g[i] = dati->d.g[i] = d.g[i];
    }
    dati->new_d.incumbent = &dati->incombenti;
    dati->d.incumbent = &dati->incombenti;
    dati->new_d.current = &dati->correnti;
    dati->d.current = &dati->correnti;
    dati->d.domains = &dati->preallocated_lists[d_domain];
    dati->new_d.domains = &dati->preallocated_lists[new_d_domain];

    dati->level = level;
    array_copy(&dati->vertex[0], vertex, arguments.arg_num);

    dati->bd = &dati->d.domains->vals[b_pos];
    /*
    for (int i=0; i < min(GRAPH_MAX, (int)bd_noRef.len.size()); i++) {
        dati->bd.len[i] = bd_noRef.len[i];
    }

    for(int i=0; i < bd_noRef.vv.size(); i++) {
        dati->bd.vv[i] = bd_noRef.vv[i];
    }

    dati->bd.is_adjacent = bd_noRef.is_adjacent;
    */

    int dim = coda_dati.size() - 1;
#if !USING_WINDOWS
    sq.push(dim);
#endif

#if USING_WINDOWS
    tasks.push_back(
        create_task(
            [dim]() -> void {
                return win_parallel_solve(dim);
            }
        )
    );
#endif
}

void nuova_starting_solve(struct D &d, struct D &new_d, struct Bidomain *bd, int vertex[], int level, vector<struct BidomainList> &preallocated_lists, int bidomain_position) {
    int ng, r, c;

    if (timeout == 1/* || cts.get_token().is_canceled()*/)
        return;

    r = (level - 1) % arguments.arg_num;
    c = (level - 1) / arguments.arg_num;

    if (r == 0) {
        if (arguments.verbose) show(d);
        nodes++;

        if (largest_incumbent.update(d.current->len)) {
            set_incumbent(d.current, d.incumbent);
            if (!arguments.quiet) printf("Incumbent size: %d\n", d.incumbent->len);
        }

        if (d.current->len + calc_bound(d.domains) <= largest_incumbent.value)
            return;

        int bidomain_position = select_bidomain_position(d.domains, d.current->len);
        if (bidomain_position == -1)   // In the MCCS case, there may be nothing we can branch on
            return;
        struct Bidomain* bd = &d.domains->vals[bidomain_position];
        struct Bidomain bd_copy;

        vertex[0] = find_and_remove_min_value(bd->vv[0], bd->len[0]);
        if (bd->len[0] == 0) {
            bd_copy = *bd;
            remove_bidomain(d.domains, bd);
            bd = &bd_copy;
        }


        struct D new_d = d;
        new_d.domains = &preallocated_lists[c + 1];

        nuova_starting_solve(d, new_d, bd, vertex, level+1, preallocated_lists, bidomain_position);
        //crea_dato(d, new_d, *bd, vertex, level + arguments.arg_num, 0, 1, preallocated_lists);
        nuova_starting_solve(d, new_d, bd, vertex, level + arguments.arg_num, preallocated_lists, bidomain_position);
#if DEBUG
        printf ("0 >>> D AFTER FILTER DOMAINS\n"); show(d);
        printf ("0 >>> NEW_D AFTER FILTER DOMAINS\n"); show(new_d);
#endif
    }

    if (r>0) {
        // Try assigning v to each vertex w in bd->right_vv, in turn
        bd->len[r]--;
        vertex[r] = -1;
        for (int i = 0; i <= bd->len[r]; i++) {
            int idx = index_of_next_smallest(bd->vv[r], bd->len[r] + 1, vertex[r]);
            vertex[r] = bd->vv[r][idx];

            // swap w with the value just past the end of the right_vv array
            bd->vv[r][idx] = bd->vv[r][bd->len[r]];
            bd->vv[r][bd->len[r]] = vertex[r];

            if (r == arguments.arg_num - 1) {

                filter_domains(d.domains, new_d.domains, &d.g[0], vertex);
                for (ng = 0; ng < arguments.arg_num; ng++) {
                    d.current->vals[d.current->len].v[ng] = vertex[ng];
                }
                d.current->len++;


                new_d.current = d.current;
                crea_dato(new_d, new_d, *bd, vertex, level+1, 1, 1, preallocated_lists, bidomain_position);

                d.current->len--;
            } else {
                //nuova_starting_solve(d, new_d, bd, vertex, level+1, preallocated_lists);
                crea_dato(d, new_d,  *bd, vertex, level+1, 0, 1, preallocated_lists, bidomain_position);
            }

#if DEBUG
            printf ("1 >>> D AFTER FILTER DOMAINS\n"); show(d);
            printf ("1 >>> NEW_D AFTER FILTER DOMAINS\n"); show(*dP);
#endif
        }
        bd->len[r]++;

        // Qui occorrerebbe ricorrere evitando di fare match su un nodo del grafo
        // corrente non solo saltare il nodo del primo grafo
        //solve(d, &d, &bd, vertex, level+1);
#if DEBUG
printf ("2 >>> D AFTER FILTER DOMAINS\n"); show(d);
printf ("2 >>> NEW_D AFTER FILTER DOMAINS\n"); show(*dP);
#endif
    }
}

std::vector <thread> threads;


void startParallel(int nThread, struct D &d, struct D &new_d, int vertex[], int level) {
#if !USING_WINDOWS
    for (int i = 0; i < nThread-1; i++) {
        threads.emplace_back(&parallel_solve);
    }
#endif

    vector<BidomainList> pre_alloc (biggest_graph+1);
    pre_alloc.at(0).vals = d.domains->vals;
    pre_alloc.at(0).len = d.domains->len;
    nuova_starting_solve(d, new_d, NULL, vertex, level, pre_alloc, -1);
#if USING_WINDOWS
    for (int i = 0; i < tasks.size(); i++) {
        tasks[i].wait();
    }
    return;
#else
    for (int i = 0; i < nThread; i++) {
        sq.push(-1);
    }
    parallel_solve();
    for (int i = 0; i < nThread-1; i++) {
        threads.at(i).join();
    }
#endif
}

void mcs(vector<struct Graph *> g) {
    int incumbent_size = 0;
    struct VtxPairList incumbent;
    incumbent.len=incumbent_size;
    struct BidomainList domains;

    vector<vector<int>> leftright;
    leftright.resize(GRAPH_MAX);
    vector<int> lr(GRAPH_MAX, 0);
    vector<int> start_lr(GRAPH_MAX);
    vector<int> len(GRAPH_MAX, 0);
    vector<int *> pointer(GRAPH_MAX);

    // Create a bidomain for vertices without loops (label 0),
    // and another for vertices with loops (label 1)
    for (int label = 0; label <= 1; label++) {
        for (int ng = 0; ng < arguments.arg_num; ng++) {
            start_lr[ng] = lr[ng];
            leftright[ng].resize(g.at(ng)->n + 1);
            for (int i = 0; i < g.at(ng)->n; i++)
                if (g.at(ng)->label[i] == label)
                    leftright[ng][lr[ng]++] = i;
        }

        int flag = 0;
        for (int ng = 0; ng < arguments.arg_num; ng++) {
            len[ng] = lr[ng] - start_lr[ng];
            pointer[ng] = &leftright[ng][start_lr[ng]];
            if (len[ng] != 0)
                flag++;
        }
        if (flag == arguments.arg_num) {
            add_bidomain(&domains, &pointer[0], &len[0], false);
        }
    }

    struct D d;
    for (int ng = 0; ng < arguments.arg_num; ng++) {
        d.g[ng] = g[ng];
    }
    struct VtxPairList vpl;
    vpl.len = 0;

    d.incumbent = &incumbent;
    d.current = &vpl;
    d.domains = &domains;

    vector<int> vertex(GRAPH_MAX+1);
    startParallel(N_THREADS, d, d, &vertex[0], 1);
}

bool check_sol(struct Graph *g[], struct VtxPairList *solution) {

    for (int ng = 1; ng < arguments.arg_num; ng++) {
        vector<bool> used_left(g[0]->n, false);
        vector<bool> used_right(g[ng]->n, false);

        for (int i = 0; i < solution->len; i++) {
            struct VtxPair p0 = solution->vals[i];

            if (used_left[p0.v[0]] || used_right[p0.v[ng]])
                return false;

            used_left[p0.v[0]] = true;
            used_right[p0.v[ng]] = true;

            if (g[0]->label[p0.v[0]] != g[ng]->label[p0.v[ng]])
                return false;

            for (int j = i + 1; j < solution->len; j++) {
                struct VtxPair p1 = solution->vals[j];
                if (g[0]->adjmat[p0.v[0]][p1.v[0]] != g[ng]->adjmat[p0.v[ng]][p1.v[ng]])
                    return false;
            }
        }
    }

    return true;
}

void parse(char key, char* val, int& counter) {
    switch (key) {
    case 'v':
        arguments.verbose = true;
        break;
    case 'l':
        arguments.lad = true;
        if (arguments.lad && arguments.dimacs) {
            exit(-200);
        }
        break;
    case 'd':
        arguments.dimacs = true;
        if (arguments.lad && arguments.dimacs) {
            exit(-200);
        }
        break;
    case 't':
        arguments.timeout = atoi(val);
        counter++;
        break;
    case 'a':
        arguments.filename.push_back(val);
        arguments.arg_num++;
        break;
    }
}

void arg_parse(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {

        if (argv[i][0] == '-') {
            parse(argv[i][1], argv[i + 1], i);
        }
        else
        {
            parse('a', argv[i], i);
        }

    }
}

int main(int argc, char **argv) {
    int i;

    fprintf(stdout, "### StQ Multi-Graph Version Running.\n");

#ifdef SIGALRM
    // Signal handler
    if (signal(SIGALRM, sig_alarm) == SIG_ERR) {
        fprintf(stderr, "Error signal alarm instantiation.\n");
        exit(0);
    }
#endif

    set_default_arguments();
#if USING_WINDOWS
    arg_parse(argc, argv);
#else
    argp_parse(&argp, argc, argv, 0, 0, 0);
#endif //USING_WINDOWS

    if (arguments.arg_num < 2) {
        return 1;
    }

    GRAPH_MAX = arguments.arg_num;

    vector <struct Graph *> g(GRAPH_MAX);

    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';

    for (i = 0; i < arguments.arg_num; i++) {
        //g[i] = (Graph *) calloc(1, sizeof(struct Graph));
        g[i] = readGraph(arguments.filename[i], format);
    }
    biggest_graph = MAX_N = g[0]->n;
    for (i = 1; i < arguments.arg_num; i++) {
        MAX_N = MAX_N > g[i]->n ? g[i]->n : MAX_N;
        biggest_graph = biggest_graph > g[i]->n ? biggest_graph : g[i]->n;
    }


    clock_t start = clock();

    if (arguments.timeout > 0) {
#if USING_WINDOWS
        create_task(
            []() -> void {
                Sleep(1000 * arguments.timeout);
                timeout = 1;
            }
        );
#else
#ifdef SIGALRM
        alarm(arguments.timeout);
#endif
#endif
    }

    auto beginning = steady_clock::now();

    mcs(g);

    auto stop = steady_clock::now();
    auto total_time_elapsed = duration_cast<milliseconds>(stop - beginning).count();

    int indice = -1, best = -1, size = coda_dati.size();
    for(int i=0; i<coda_dati.size(); i++) {
        if(coda_dati.at(i)->d.incumbent->len > best) {
            best = coda_dati.at(i)->d.incumbent->len;
            indice = i;
        }
    }
    if (timeout == 1) {
        cout << "Esecuzione interrotta in seguito a timeout!" << endl;
    }
    if (indice == -1) {
        cout << "Nessuna soluzione trovata!" << endl;
        return 1;
    }

    struct VtxPairList *solution = coda_dati.at(indice)->d.incumbent;

    clock_t time_elapsed = clock() - start;

    if (!check_sol(&g[0], solution))
        fprintf(stderr, "*** Error: Invalid solution\n");
    //fail("*** Error: Invalid solution\n");

    printf("Solution size %d\n", solution->len);
    for (int i = 0; i < g[0]->n; i++)
        for (int j = 0; j < solution->len; j++)
            if (solution->vals[j].v[0] == i)
                for (int k = 0; k < arguments.arg_num; k++) {
                    if (k == 0)
                        printf("(%d", solution->vals[j].v[k]);
                    else if (k == arguments.arg_num - 1)
                        printf(" -> %d) ", solution->vals[j].v[k]);
                    else
                        printf(" -> %d", solution->vals[j].v[k]);
                }
    printf("\n");

    setlocale(LC_NUMERIC, "");
    printf("Nodes:                      %15ld\n", nodes);
    printf("CPU time (ms):              %15ld\n", time_elapsed * 1000 / CLOCKS_PER_SEC);

    printf(">> ");
    for (i = 0; i < arguments.arg_num; i++) {
        printf("%s ", arguments.filename[i]);
    }
    if (timeout == 0) {
        printf("- NGraph=%d - OK - len=%d - time=%015.10f\n", arguments.arg_num, solution->len,
               (float) time_elapsed / CLOCKS_PER_SEC);
    } else {
        printf("- NGraph=%d - TimeOut - len=%d - time=%015.10f\n", arguments.arg_num, solution->len,
               (float) time_elapsed / CLOCKS_PER_SEC);
    }

    for (i = 0; i < arguments.arg_num; i++) {
        free(g[i]);
    }
    for (i=0; i<coda_dati.size(); i++) {
        free(coda_dati[i]);
    }

    printf(">>> %d - %f\n", solution->len, (float) total_time_elapsed/1000);

    return 0;
}

void array_copy(int dst[], int src[], int size) {
    int i;

    for (i = 0; i <= size; i++)
        dst[i] = src[i];
}
