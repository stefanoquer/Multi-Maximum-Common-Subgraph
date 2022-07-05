#define _GNU_SOURCE
#define _POSIX_SOURCE

#include "graph.h"

#include <argp.h>
#include <limits.h>
#include <locale.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <signal.h>
#include <unistd.h>

#include <locale.h>

typedef unsigned long long ULL;

int **malloc2d (int, int);
void free2d (int **, int);
void array_copy (int [], int [], int);

int timeout = 0;

#define DEBUG 0
#define GRAPH_MAX 20

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

/*
 * Signal handler (for timeout)
 */

void sig_alarm (int sig) {
  fprintf (stderr, "Time out reached.\n");
  timeout = 1;

  return;
}

void swap(int *a, int *b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}
static void fail(char* msg) {
    printf("%s\n", msg);
    exit(1);
}

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Find a maximum clique in a graph in DIMACS format";
static char args_doc[] = "FILENAME 1..n";
static struct argp_option options[] = {
    {"quiet", 'q', 0, 0, "Quiet output"},
    {"verbose", 'v', 0, 0, "Verbose output"},
    {"dimacs", 'd', 0, 0, "Read DIMACS format"},
    {"lad", 'l', 0, 0, "Read LAD format"},
    {"connected", 'c', 0, 0, "Solve max common CONNECTED subgraph problem"},
    {"timeout", 't', "timeout", 0, "Specify a timeout (seconds)"},
    { 0 }
};

static struct {
    bool quiet;
    bool verbose;
    bool connected;
    bool dimacs;
    bool lad;
    char *filename[GRAPH_MAX];
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
    for (i=0; i<GRAPH_MAX; i++)
      arguments.filename[i] = NULL;
    arguments.timeout = 0;
    arguments.arg_num = 0;
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
            arguments.connected = true;
            break;
        case 't':
            arguments.timeout = atoi(arg);
            break;
        case ARGP_KEY_ARG:
            arguments.filename[arguments.arg_num] = arg;
            //argp_usage(state);
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
                                     Stats
*******************************************************************************/

static long nodes = 0;

/*******************************************************************************
                                 MCS functions
*******************************************************************************/

struct VtxPair {
    int v[GRAPH_MAX];
};

struct VtxPairList {
    struct VtxPair vals[MAX_N];
    int len;
};

struct Bidomain {
    int *vv[GRAPH_MAX];
    int len[GRAPH_MAX];
    bool is_adjacent;
};

struct BidomainList {
    struct Bidomain vals[MAX_N*2];
    int len;
};

struct D {
    struct Graph *g[GRAPH_MAX];
    struct VtxPairList *incumbent;
    struct VtxPairList *current;
    struct BidomainList *domains;
};

struct BidomainList preallocated_lists[MAX_N];

int calc_bound(struct BidomainList *domains) {
    int bound = 0;
    for (int i=0; i<domains->len; i++) {
        struct Bidomain *bd = &domains->vals[i];
	int min_size = INT_MAX;
        for (int j=0; j<arguments.arg_num; j++) {
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
    for (int i=0; i<domains->len; i++) {
        struct Bidomain *bd = &domains->vals[i];
        if (arguments.connected && current_matching_size>0 && !bd->is_adjacent) continue;
        int len = bd->len[0];
        for (int j=1; j<arguments.arg_num; j++) {
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

// Returns length of left half of array
// DA NON CAMBIARE
int partition(int *vv, int vv_len, unsigned char *adjrow) {
    int i=0;
    for (int j=0; j<vv_len; j++) {
        int adj = adjrow[vv[j]];
        if (adj) {
            swap(&vv[i], &vv[j]);
            i++;
        }
    }
    return i;
}

void add_bidomain (
  struct BidomainList *bd_list, int *vvMat[GRAPH_MAX], int lenV[GRAPH_MAX], bool is_adjacent
) {
    for (int ng=0; ng<arguments.arg_num; ng++) {
      bd_list->vals[bd_list->len].vv[ng] = vvMat[ng];
      bd_list->vals[bd_list->len].len[ng] = lenV[ng];
    }
    bd_list->vals[bd_list->len].is_adjacent = is_adjacent;
    bd_list->len++;
}

/* dati domain d crea domanin new_d dato l'accoppiamento v-w */
void filter_domains(
  struct BidomainList *d, struct BidomainList *new_d,
  struct Graph *g[], int vertex[]
) {
    int flag_len_edge, flag_len_noedge;
    int len_edge[GRAPH_MAX], len_noedge[GRAPH_MAX];
    
    new_d->len=0;
    for (int j=0; j<d->len; j++) {
        struct Bidomain *old_bd = &d->vals[j];

	flag_len_edge = flag_len_noedge = 0;
	for (int i=0; i<arguments.arg_num; i++) {
          len_edge[i] = partition (old_bd->vv[i], old_bd->len[i], g[i]->adjmat[vertex[i]]);
          len_noedge[i] = old_bd->len[i] - len_edge[i];

	  if (len_edge[i])
	    flag_len_edge++;
 	  if (len_noedge[i])
	    flag_len_noedge++;
        }
	
	// Se ci sono vertici NON adiancenti su g e h aggiungo il bidomain
	if (flag_len_noedge==arguments.arg_num) {
	   int *vv[GRAPH_MAX];
           for (int i=0; i<arguments.arg_num; i++) {
              vv[i] = old_bd->vv[i] + len_edge[i];
	   }
           add_bidomain (new_d, vv, len_noedge, old_bd->is_adjacent);
        }

	// Se ci sono vertici adiacenti su g e h aggiungo il bidomain	
        if (flag_len_edge==arguments.arg_num) {
           add_bidomain(new_d, old_bd->vv, len_edge, true);
	}
    }
}

void remove_bidomain(struct BidomainList *list, struct Bidomain *b) {
    int i = b - &list->vals[0];
    list->vals[i] = list->vals[list->len-1];
    list->len--;
}

void show(struct D d) {
    printf("Nodes: %ld\n", nodes);
    printf("Length of current assignment: %d\n", d.current->len);
    printf("Current assignment:");
    for (int i=0; i<d.current->len; i++) {
        for (int j=0; j<arguments.arg_num; j++) {
          if (j==0) 
	    printf("  %d", d.current->vals[i].v[j]);
 	  else
	    printf("->%d", d.current->vals[i].v[j]);
        }
    }
    printf("\n");
    for (int i=0; i<d.domains->len; i++) {
        struct Bidomain bd = d.domains->vals[i];
        for (int ng=0; ng<arguments.arg_num; ng++) {
          printf("Graph %d  ", ng);
          for (int j=0; j<bd.len[ng]; j++)
            printf("%d ", bd.vv[ng][j]);
          printf("\n");
	}
    }
    printf("\n\n");
}

void set_incumbent(struct VtxPairList *current, struct VtxPairList *incumbent) {
    incumbent->len = current->len;
    for (int i=0; i<current->len; i++)
        incumbent->vals[i] = current->vals[i];
}

int find_and_remove_min_value(int *arr, int *len) {
    int min_v = INT_MAX;
    int idx = -1;
    for (int i=0; i<*len; i++) {
        if (arr[i] < min_v) {
            min_v = arr[i];
            idx = i;
        }
    }
    swap(&arr[idx], &arr[*len-1]);
    (*len)--;
    return min_v;
}

// returns the index of the smallest value in arr that is >w.
// Assumption: such a value exists
// Assumption: arr contains no duplicates
// Assumption: arr has no values==INT_MAX
int index_of_next_smallest(int *arr, int len, int w) {
    int idx = -1;
    int smallest = INT_MAX;
    for (int i=0; i<len; i++) {
        if (arr[i]>w && arr[i]<smallest) {
            smallest = arr[i];
            idx = i;
        }
    }
    return idx;
}

void solve(struct D d, struct D new_d, struct Bidomain *bd, int vertex[], int level) {
  int ng, r, c;

  if (timeout==1)
    return;

  r = (level-1) % arguments.arg_num;
  c = (level-1) / arguments.arg_num;

  if (r==0) {
    if (arguments.verbose) show(d);
    nodes++;

    if (d.current->len > d.incumbent->len) {
      set_incumbent(d.current, d.incumbent);
      if (!arguments.quiet) printf("Incumbent size: %d\n", d.incumbent->len);
    }

    if (d.current->len + calc_bound(d.domains) <= d.incumbent->len)
      return;

    struct Bidomain *bd = select_bidomain(d.domains, d.current->len);;
    if (bd == NULL)   // In the MCCS case, there may be nothing we can branch on
       return;
    struct Bidomain bd_copy;

    vertex[0] = find_and_remove_min_value(bd->vv[0], &bd->len[0]);
    if (bd->len[0] == 0) {
        bd_copy = *bd;
        remove_bidomain(d.domains, bd);
        bd = &bd_copy;
    }

    struct D new_d = d;
    new_d.domains = &preallocated_lists[c+1];

    solve(d, new_d, bd, vertex, level+1);
    solve(d, new_d, bd, vertex, level+arguments.arg_num);
#if DEBUG
    printf ("0 >>> D AFTER FILTER DOMAINS\n"); show(d);
    printf ("0 >>> NEW_D AFTER FILTER DOMAINS\n"); show(new_d);
#endif
  }
  
  if (r>0) {
    // Try assigning v to each vertex w in bd->right_vv, in turn
    bd->len[r]--;
    vertex[r] = -1;
    for (int i=0; i<=bd->len[r]; i++) {
      int idx = index_of_next_smallest(bd->vv[r], bd->len[r]+1, vertex[r]);
      vertex[r] = bd->vv[r][idx];

      // swap w with the value just past the end of the right_vv array
      bd->vv[r][idx] = bd->vv[r][bd->len[r]];
      bd->vv[r][bd->len[r]] = vertex[r];

      if (r==arguments.arg_num-1) {
        filter_domains(d.domains, new_d.domains, d.g, vertex);
        for (ng=0; ng<arguments.arg_num; ng++) {
          d.current->vals[d.current->len].v[ng] = vertex[ng];
	    }
        d.current->len++;

        int new_vertex[GRAPH_MAX] = {0};
        array_copy (new_vertex, vertex, arguments.arg_num);
        solve(new_d, new_d, bd, new_vertex, level+1);

        d.current->len--;
      } else {
        solve (d, new_d, bd, vertex, level+1);
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

struct VtxPairList mcs(struct Graph *g[]) {
    int incumbent_size = 0;
    struct VtxPairList incumbent = {.len=incumbent_size};
    struct BidomainList *domains = &preallocated_lists[0];

    int leftright[GRAPH_MAX][MAX_N];
    int lr[GRAPH_MAX] = {0};
    int start_lr[GRAPH_MAX];
    int len[GRAPH_MAX] = {0};
    int *pointer[GRAPH_MAX];

    // Create a bidomain for vertices without loops (label 0),
    // and another for vertices with loops (label 1)
    for (int label=0; label<=1; label++) {
        for (int ng=0; ng<arguments.arg_num; ng++) {
	  start_lr[ng] = lr[ng];
          for (int i=0; i<g[ng]->n; i++)
            if (g[ng]->label[i]==label)
                leftright[ng][lr[ng]++] = i;
	}

        int flag = 0;
        for (int ng=0; ng<arguments.arg_num; ng++) {
	  len[ng] = lr[ng] - start_lr[ng];
          pointer[ng] = &leftright[ng][start_lr[ng]];
	  if (len[ng]!=0)
	    flag++;
	}
        if (flag==arguments.arg_num) {
          add_bidomain (domains, pointer, len, false);
        }
    }

    struct D d;
    for (int ng=0; ng<arguments.arg_num; ng++) {
      d.g[ng] = g[ng];
    }
    d.incumbent=&incumbent,
    d.current = &(struct VtxPairList){.len = 0};
    d.domains = domains;

    int vertex[GRAPH_MAX] = {0};
    solve(d, d, NULL, vertex, 1);

    return incumbent;
}

bool check_sol(struct Graph *g[], struct VtxPairList *solution) {
  bool used_left[MAX_N];
  bool used_right[MAX_N];

    for (int ng=1; ng<arguments.arg_num; ng++) {
      for (int i=0; i<MAX_N; i++)
        used_left[i] = used_right[i] = false;
      
      for (int i=0; i<solution->len; i++) {
        struct VtxPair p0 = solution->vals[i];
	
        if (used_left[p0.v[0]] || used_right[p0.v[ng]])
            return false;

        used_left[p0.v[0]] = true;
        used_right[p0.v[ng]] = true;

        if (g[0]->label[p0.v[0]] != g[ng]->label[p0.v[ng]])
            return false;

        for (int j=i+1; j<solution->len; j++) {
            struct VtxPair p1 = solution->vals[j];
            if (g[0]->adjmat[p0.v[0]][p1.v[0]] != g[ng]->adjmat[p0.v[ng]][p1.v[ng]])
                return false;
        }
      }	
    }
    
    return true;
}

int main(int argc, char** argv) {
    setlocale(LC_ALL, "en_GB");
    int i;

    fprintf (stdout, "### StQ Multi-Graph Version Running.\n");

     // Signal handler
    if (signal (SIGALRM, sig_alarm) == SIG_ERR) {
      fprintf (stderr, "Error signal alarm instantiation.\n");
      exit (0);
    }

    set_default_arguments();
    argp_parse(&argp, argc, argv, 0, 0, 0);

    struct Graph *g[GRAPH_MAX];

    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';

    for (i=0; i<arguments.arg_num; i++) {
      g[i] = calloc (1, sizeof (struct Graph));
      readGraph (arguments.filename[i], g[i], format);
    }

    
    clock_t start = clock();
    if (arguments.timeout > 0) {
      alarm (arguments.timeout);
    }
    
    struct VtxPairList solution = mcs(g);

    clock_t time_elapsed = clock() - start;

    if (!check_sol(g, &solution))
      fprintf(stderr, "*** Error: Invalid solution\n");
      //fail("*** Error: Invalid solution\n");

    printf("Solution size %d\n", solution.len);
    for (int i=0; i<g[0]->n; i++)
        for (int j=0; j<solution.len; j++)
	  if (solution.vals[j].v[0] == i)
	    for (int k=0; k<arguments.arg_num; k++) {
	      if (k==0)
                printf("(%d", solution.vals[j].v[k]);
	      else if (k==arguments.arg_num-1)
                printf(" -> %d) ", solution.vals[j].v[k]);
	      else
                printf(" -> %d", solution.vals[j].v[k]);
	    }
    printf("\n");

    setlocale(LC_NUMERIC, "");
    printf("Nodes:                      %'15ld\n", nodes);
    printf("CPU time (ms):              %15ld\n", time_elapsed * 1000 / CLOCKS_PER_SEC);

    printf(">> ");
    for (i=0; i<arguments.arg_num; i++) {
      printf ("%s ", arguments.filename[i]);
    }
    if (timeout==0) {
      printf("- NGraph=%d - OK - len=%d - time=%015.10f\n", arguments.arg_num, solution.len, (float)time_elapsed/CLOCKS_PER_SEC);
    } else {
      printf("- NGraph=%d - TimeOut - len=%d - time=%015.10f\n", arguments.arg_num, solution.len, (float)time_elapsed/CLOCKS_PER_SEC);
    }

    for (i=0; i<arguments.arg_num; i++) {
      free(g[i]);
    }

    printf(">>> %d - %d.%d\n", solution.len, time_elapsed/CLOCKS_PER_SEC,time_elapsed%CLOCKS_PER_SEC);
    
    return 0;
}

int **malloc2d (int r, int c) {
  int i;
  int **mat;

  mat = (int **) malloc (r * sizeof (int *));
  if (mat == NULL) {
    fprintf (stderr, "Memory allocation error.\n");
    exit(EXIT_FAILURE);
  }
  for (i=0; i<r; i++) {
    mat[i] = (int *) malloc(c * sizeof (int));
    if (mat[i]==NULL) {
      fprintf (stderr, "Memory allocation error.\n");
      exit(EXIT_FAILURE);
    }
  }

  return mat;
}

void free2d (int **mat, int r) {
  int i;

  for (i=0; i<r; i++) {
    free (mat[i]);
  }
  free (mat);

  return;
}

void array_copy (int dst[], int src[], int size) {
  int i;

  for (i=0; i<=size; i++)
    dst[i] = src[i];
}
