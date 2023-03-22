#include "graph.h"

#include <algorithm>
#include <functional>
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
#include <cstdlib>
#include <argp.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>


using std::vector;
using std::cout;
using std::endl;

using std::chrono::steady_clock;
using std::chrono::milliseconds;
using std::chrono::duration_cast;

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
    { "debug", 'D', 0, 0, "Print more data for debug" },
    { "buffer", 'B', "buffer_size", 0, "Size of the buffer, more RAM for better performances (on the GPU side)" },
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
    bool debug;
    int buffer_size;
} arguments;

static std::atomic<bool> abort_due_to_timeout;

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
    arguments.debug = false;
    arguments.buffer_size = 1;
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
        break;case 'D':
	    arguments.debug = true;
	    break;
	case 'B':
	    arguments.buffer_size = std::stoi(arg);
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
                                 GPU definitions
*******************************************************************************/


#define L   0
#define R   1
#define LL  2
#define RL  3
#define ADJ 4
#define P   5
#define W   6
#define IRL 7

#define BDS 8

#define START 0
#define END 1

#define MIN(a, b) (a < b)? a : b
#define MAX(a, b) (a > b)? a : b

#define N_BLOCKS 368
//#define N_BLOCKS 184
#define BLOCK_SIZE 1024
//#define N_BLOCKS 8
//#define BLOCK_SIZE 32
#define MAX_GRAPH_SIZE 64
#define checkCudaErrors(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)

#define CLEAN 100
uchar time_to_clean = 0;

__constant__ uchar d_adjmat0[MAX_GRAPH_SIZE][MAX_GRAPH_SIZE];
__constant__ uchar d_adjmat1[MAX_GRAPH_SIZE][MAX_GRAPH_SIZE];
__constant__ uchar d_n0;
__constant__ uchar d_n1;

uchar adjmat0[MAX_GRAPH_SIZE][MAX_GRAPH_SIZE];
uchar adjmat1[MAX_GRAPH_SIZE][MAX_GRAPH_SIZE];
uchar n0;
uchar n1;

const uint __gpu_level = 5;
struct timespec start;
double ttcd = 0, ttsg = 0;
Graph *g0, *g1;

std::thread::id main_th_id;

/*******************************************************************************
                                 GPU functions
*******************************************************************************/



__host__ __device__
void uchar_swap(uchar *a, uchar *b){
	uchar tmp = *a;
	*a = *b;
	*b = tmp;
}

__host__ __device__
uchar select_next_v(uchar *left, uchar *bd){
	uchar min = UCHAR_MAX, idx = UCHAR_MAX;
	if(bd[RL] != bd[IRL])
		return left[bd[L] + bd[LL]];
	for (uchar i = 0; i < bd[LL]; i++)
		if (left[bd[L] + i] < min) {
			min = left[bd[L] + i];
			idx = i;
		}
	uchar_swap(&left[bd[L] + idx], &left[bd[L] + bd[LL] - 1]);
	bd[LL]--;
	bd[RL]--;
	return min;
}


__host__ __device__
uchar select_next_w(uchar *right, uchar *bd) {
	uchar min = UCHAR_MAX, idx = UCHAR_MAX;
	for (uchar i = 0; i < bd[RL]+1; i++)
		if ((right[bd[R] + i] > bd[W] || bd[W] == UCHAR_MAX)
				&& right[bd[R] + i] < min) {
			min = right[bd[R] + i];
			idx = i;
		}
	if(idx == UCHAR_MAX)
		bd[RL]++;
	return idx;
}


__host__  __device__ uchar index_of_next_smallest(const uchar *arr,
		uchar start_idx, uchar len, uchar w) {
	uchar idx = UCHAR_MAX;
	uchar smallest = UCHAR_MAX;
	for (uchar i = 0; i < len; i++) {
		if ((arr[start_idx + i] > w || w == UCHAR_MAX)
				&& arr[start_idx + i] < smallest) {
			smallest = arr[start_idx + i];
			idx = i;
		}
	}
	return idx;
}

__host__  __device__ uchar find_min_value(const uchar *arr, uchar start_idx,
		uchar len) {
	uchar min_v = UCHAR_MAX;
	for (int i = 0; i < len; i++) {
		if (arr[start_idx + i] < min_v)
			min_v = arr[start_idx + i];
	}
	return min_v;
}

__host__ __device__
void remove_from_domain(uchar *arr, const uchar *start_idx, uchar *len,
		uchar v) {
	int i = 0;
	for (i = 0; arr[*start_idx + i] != v; i++)
		;
	uchar_swap(&arr[*start_idx + i], &arr[*start_idx + *len - 1]);
	(*len)--;
}

__host__ __device__
void update_incumbent(uchar cur[][2], uchar inc[][2], uchar cur_pos,
		uchar *inc_pos) {
	if (cur_pos > *inc_pos) {
		*inc_pos = cur_pos;
		for (int i = 0; i < cur_pos; i++) {
			inc[i][L] = cur[i][L];
			inc[i][R] = cur[i][R];
		}
	}
}

// BIDOMAINS FUNCTIONS /////////////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__
void add_bidomain(uchar domains[][BDS], uint *bd_pos, uchar left_i,
		uchar right_i, uchar left_len, uchar right_len, uchar is_adjacent,
		uchar cur_pos) {
	domains[*bd_pos][L] 	= left_i;
	domains[*bd_pos][R] 	= right_i;
	domains[*bd_pos][LL] 	= left_len;
	domains[*bd_pos][RL] 	= right_len;
	domains[*bd_pos][ADJ] 	= is_adjacent;
	domains[*bd_pos][P] 	= cur_pos;
	domains[*bd_pos][W] 	= UCHAR_MAX;
	domains[*bd_pos][IRL] 	= right_len;

	(*bd_pos)++;
}

__host__  __device__ uint calc_bound(uchar domains[][BDS], uint bd_pos,
		uint cur_pos, uint *bd_n) {
	uint bound = 0;
	int i;
	for (i = bd_pos - 1; i >= 0 && domains[i][P] == cur_pos; i--)
		bound += MIN(domains[i][LL], domains[i][IRL]);
	*bd_n = bd_pos - 1 - i;
	return bound;
}

__host__  __device__ uchar partition(uchar *arr, uchar start, uchar len,
		const uchar *adjrow) {
	uchar i = 0;
	for (uchar j = 0; j < len; j++) {
		if (adjrow[arr[start + j]]) {
			uchar_swap(&arr[start + i], &arr[start + j]);
			i++;
		}
	}
	return i;
}

__host__  __device__
uchar find_min_value(uchar *arr, uchar start_idx, uchar len){
	uchar min_v = UCHAR_MAX;
	for(int i = 0; i < len; i++){
		if(arr[start_idx+i] < min_v)
			min_v = arr[start_idx + i];
	}
	return min_v;
}

__host__  __device__
void select_bidomain(uchar domains[][BDS], uint bd_pos,  uchar *left, int current_matching_size, bool connected){
	int i;
	uint min_size = UINT_MAX;
	uint min_tie_breaker = UINT_MAX;
	uint best = UINT_MAX;
	uchar *bd;
	for (i = bd_pos - 1, bd = &domains[i][L]; i >= 0 && bd[P] == current_matching_size; i--, bd = &domains[i][L]) {
		if (connected && current_matching_size>0 && !bd[ADJ]) continue;
		int len = bd[LL] > bd[RL] ? bd[LL] : bd[RL];
		if (len < min_size) {
			min_size = len;
			min_tie_breaker = find_min_value(left, bd[L], bd[LL]);
			best = i;
		} else if (len == min_size) {
			int tie_breaker = find_min_value(left, bd[L], bd[LL]);
			if (tie_breaker < min_tie_breaker) {
				min_tie_breaker = tie_breaker;
				best = i;
			}
		}
	}
	if(best != UINT_MAX && best != bd_pos-1){
		uchar tmp[BDS];
		for(i = 0; i < BDS; i++) tmp[i] = domains[best][i];
		for(i = 0; i < BDS; i++) domains[best][i] = domains[bd_pos-1][i];
		for(i = 0; i < BDS; i++) domains[bd_pos-1][i] = tmp[i];

	}
}



__device__
void d_generate_next_domains(uchar domains[][BDS], uint *bd_pos, uint cur_pos, uchar *left, uchar *right, uchar v, uchar w, uint inc_pos) {
	int i;
	uint bd_backup = *bd_pos;
	uint bound = 0;
	uchar *bd;
	for (i = *bd_pos - 1, bd = &domains[i][L]; i >= 0 && bd[P] == cur_pos - 1; i--, bd = &domains[i][L]) {

		uchar l_len = partition(left, bd[L], bd[LL], d_adjmat0[v]);
		uchar r_len = partition(right, bd[R], bd[RL], d_adjmat1[w]);

		if (bd[LL] - l_len && bd[RL] - r_len) {
			add_bidomain(domains, bd_pos, bd[L] + l_len, bd[R] + r_len, bd[LL] - l_len, bd[RL] - r_len, bd[ADJ], (uchar) (cur_pos));
			bound += MIN(bd[LL] - l_len, bd[RL] - r_len);
		}
		if (l_len && r_len) {
			add_bidomain(domains, bd_pos, bd[L], bd[R], l_len, r_len, true, (uchar) (cur_pos));
			bound += MIN(l_len, r_len);
		}
	}
	if (cur_pos + bound <= inc_pos)
		*bd_pos = bd_backup;
}

__global__
void d_mcs(uchar *args, int n_threads, uchar a_size, uint *args_i, uint actual_inc, uchar *device_solutions, uint max_sol_size, uint last_arg, bool verbose, bool connected, uint *glo_sh_inc) {
	uint my_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
	uchar cur[MAX_GRAPH_SIZE][2], incumbent[MAX_GRAPH_SIZE][2],
	domains[MAX_GRAPH_SIZE * MAX_GRAPH_SIZE / 2][BDS], left[MAX_GRAPH_SIZE],
	right[MAX_GRAPH_SIZE], v, w;
	uint bd_pos = 0, bd_n = 0;
	uchar inc_pos = 0;
	__shared__ uint sh_inc;
	sh_inc = actual_inc;
	uint loc_sh_inc = actual_inc;
	
	__syncthreads();
	if (my_idx < n_threads) {
		for (int i = args_i[my_idx]; ( my_idx < n_threads-1 &&  i < args_i[my_idx +1] ) || ( my_idx == n_threads-1 && i < last_arg );) {
			add_bidomain(domains, &bd_pos, args[i++], args[i++], args[i++], args[i++], args[i++], args[i++]);
			for (int p = 0; p < domains[bd_pos - 1][P]; p++) {
				cur[p][L] = args[i++];
			}
			for (int p = 0; p < domains[bd_pos - 1][P]; p++) {
				cur[p][R] = args[i++];
			}
			for (int l = 0; l < d_n0; l++) {
				left[l] = args[i++];
			}
			for (int r = 0; r < d_n1; r++) {
				right[r] = args[i++];
			}
		}
		while (bd_pos > 0) {
			uchar *bd = &domains[bd_pos - 1][L];
			if (calc_bound(domains, bd_pos, bd[P], &bd_n) + bd[P] + (bd[RL] != bd[IRL]) <= loc_sh_inc || (bd[LL] == 0 && bd[RL] == bd[IRL])) {
				bd_pos--;
			} else {
				select_bidomain(domains, bd_pos, left, domains[bd_pos - 1][P], connected);
				if (bd[RL] == bd[IRL]) {
					v = find_min_value(left, bd[L], bd[LL]);
					remove_from_domain(left, &bd[L], &bd[LL], v);
					bd[RL]--;
				} else {
					v = left[bd[L] + bd[LL]];
				}
				if ((bd[W] = index_of_next_smallest(right, bd[R], bd[RL] + (uchar) 1, bd[W])) == UCHAR_MAX) {
					bd[RL]++;
				} else {
					w = right[bd[R] + bd[W]];
					right[bd[R] + bd[W]] = right[bd[R] + bd[RL]];
					right[bd[R] + bd[RL]] = w;
					bd[W] = w;
					cur[bd[P]][L] = v;
					cur[bd[P]][R] = w;
					update_incumbent(cur, incumbent, bd[P] + 1, &inc_pos);
					loc_sh_inc = MAX(atomicMax(&sh_inc, inc_pos), loc_sh_inc);
					loc_sh_inc = MAX(atomicMax(glo_sh_inc, loc_sh_inc), loc_sh_inc);
					d_generate_next_domains(domains, &bd_pos, bd[P] + 1, left, right, v, w, inc_pos);
				}
			}
		}
	}
	device_solutions[blockIdx.x* max_sol_size] = 0;

	__syncthreads();
	if (atomicCAS(&sh_inc, inc_pos, 0) == inc_pos && inc_pos > 0) {
		if(verbose) printf("Th_%d found new solution of size %d\n", my_idx, inc_pos);
		bd_pos = 0;
		device_solutions[blockIdx.x* max_sol_size + bd_pos++] = inc_pos;
		for (int i = 0; i < inc_pos; i++)
			device_solutions[blockIdx.x* max_sol_size + bd_pos++] = incumbent[i][L];
		for (int i = 0; i < inc_pos; i++)
			device_solutions[blockIdx.x* max_sol_size + bd_pos++] = incumbent[i][R];
	}
}

double compute_elapsed_sec(struct timespec strt){
	struct timespec now;
	double time_elapsed;

	clock_gettime(CLOCK_MONOTONIC, &now);
	time_elapsed = (now.tv_sec - strt.tv_sec);
	time_elapsed += (double)(now.tv_nsec - strt.tv_nsec) / 1000000000.0;

	return time_elapsed;
}

static void CheckCudaErrorAux(const char *file, unsigned line,
		const char *statement, cudaError_t err) {
	if (err == cudaSuccess)
		return;
	fprintf(stderr, "%s returned %s(%d) at %s:%d\n", statement,
			cudaGetErrorString(err), err, file, line);
	exit(1);
}

void move_graphs_to_gpu(Graph *graph0, Graph *graph1) {
	n0 = graph0->n;
	n1 = graph1->n;
	g0 = graph0;
	g1 = graph1;
	checkCudaErrors(cudaMemcpyToSymbol(d_n0, &graph0->n, sizeof(uchar)));
	checkCudaErrors(cudaMemcpyToSymbol(d_n1, &graph1->n, sizeof(uchar)));
	checkCudaErrors(cudaMemcpyToSymbol(d_adjmat0, adjmat0, MAX_GRAPH_SIZE*MAX_GRAPH_SIZE));
	checkCudaErrors(cudaMemcpyToSymbol(d_adjmat1, adjmat1, MAX_GRAPH_SIZE*MAX_GRAPH_SIZE));
}


void h_generate_next_domains(uchar domains[][BDS], uint *bd_pos, uint cur_pos,
		uchar *left, uchar *right, uchar v, uchar w, uint inc_pos) {
	int i;
	uint bd_backup = *bd_pos;
	uint bound = 0;
	uchar *bd;
	for (i = *bd_pos - 1, bd = &domains[i][L]; i >= 0 && bd[P] == cur_pos - 1;
			i--, bd = &domains[i][L]) {

		uchar l_len = partition(left, bd[L], bd[LL], adjmat0[v]);
		uchar r_len = partition(right, bd[R], bd[RL], adjmat1[w]);

		if (bd[LL] - l_len && bd[RL] - r_len) {
			add_bidomain(domains, bd_pos, bd[L] + l_len, bd[R] + r_len,
					bd[LL] - l_len, bd[RL] - r_len, bd[ADJ], (uchar) (cur_pos));
			bound += MIN(bd[LL] - l_len, bd[RL] - r_len);
		}
		if (l_len && r_len) {
			add_bidomain(domains, bd_pos, bd[L], bd[R], l_len, r_len, true,
					(uchar) (cur_pos));
			bound += MIN(l_len, r_len);
		}
	}
	if (cur_pos + bound <= inc_pos)
		*bd_pos = bd_backup;
}


bool check_sol(Graph *g0, Graph *g1, uchar sol[][2], uint sol_len) {
	for (int i = 0; i < sol_len; i++) {
		if (g0->label[sol[i][L]] != g1->label[sol[i][R]]) {
			printf("g0:%d and g1:%d have different labels\n", sol[i][L],
					sol[i][R]);
			return false;
		}
		for (int j = i + 1; j < sol_len; j++) {
			if (g0->adjmat[sol[i][L]][sol[j][L]]
			                          != g1->adjmat[sol[i][R]][sol[j][R]]) {
				printf("g0(%d-%d) is different than g1(%d-%d)\n", sol[i][L],
						sol[j][L], sol[i][R], sol[j][R]);
				return false;
			}
		}
	}
	return true;
}

void launch_kernel(uchar *args, int n_threads, uchar a_size, uint sol_size, uint *args_i,
		uchar incumbent[][2], uchar *inc_pos, uint total_args_size, uint last_arg) {
		
	uchar *device_args;
	uchar *device_solutions;
	uchar *host_solutions;
	uint *device_args_i;
	uint max_sol_size = 1 + 2 * (MIN(n0, n1));
	struct timespec sleep;
	sleep.tv_sec = 0;
	sleep.tv_nsec = 2000;
	cudaEvent_t stop;
	
	uint *glo_sh_inc;

	host_solutions = (uchar*) malloc(N_BLOCKS * max_sol_size * sizeof *host_solutions);

	checkCudaErrors(cudaEventCreate(&stop));
	
	checkCudaErrors(cudaMalloc(&glo_sh_inc, sizeof *glo_sh_inc));
	checkCudaErrors(cudaMemset(glo_sh_inc, 0, sizeof *glo_sh_inc));

	checkCudaErrors(cudaMalloc(&device_args, total_args_size * sizeof *device_args));
	checkCudaErrors(cudaMalloc(&device_solutions, N_BLOCKS * max_sol_size * sizeof *device_solutions));


	checkCudaErrors(cudaMemcpy(device_args, args, total_args_size * sizeof *device_args, cudaMemcpyHostToDevice));

	checkCudaErrors(cudaMalloc(&device_args_i, N_BLOCKS * BLOCK_SIZE * sizeof *device_args_i));
	checkCudaErrors(cudaMemcpy(device_args_i, args_i, N_BLOCKS * BLOCK_SIZE * sizeof *device_args_i, cudaMemcpyHostToDevice));

	if(arguments.verbose) printf("Launching kernel... %d threads\n", n_threads);

	d_mcs<<<N_BLOCKS, BLOCK_SIZE>>>(device_args, n_threads, a_size, device_args_i, *inc_pos, device_solutions, max_sol_size, last_arg, arguments.verbose, arguments.connected, glo_sh_inc);
	checkCudaErrors(cudaEventRecord(stop));

	while(cudaEventQuery(stop) == cudaErrorNotReady){
		nanosleep(&sleep, NULL);
		if(arguments.timeout && compute_elapsed_sec(start) > arguments.timeout) {
			break;
        }
	}

	if(arguments.verbose) printf("Kernel executed...\n");

	checkCudaErrors(cudaMemcpy(host_solutions, device_solutions, N_BLOCKS * max_sol_size * sizeof *device_solutions, cudaMemcpyDeviceToHost));

	checkCudaErrors(cudaFree(device_args));
	checkCudaErrors(cudaFree(device_args_i));
	checkCudaErrors(cudaFree(glo_sh_inc));

	for(int b = 0; b < N_BLOCKS; b++){
		if (*inc_pos < host_solutions[b*max_sol_size]) {
			*inc_pos = host_solutions[b*max_sol_size];
			for (int i = 1; i < *inc_pos + 1; i++) {
				incumbent[i - 1][L] = host_solutions[b*max_sol_size + i];
				incumbent[i - 1][R] = host_solutions[b*max_sol_size + *inc_pos + i];
				if(arguments.verbose) printf("|%d %d| ", incumbent[i-1][L], incumbent[i-1][R]);
			}if(arguments.verbose) printf("\n");
		}
	}
	free(host_solutions);
	checkCudaErrors(cudaFree(device_solutions));
}

struct Data {
	uchar domains[BDS - 2];
	uchar cur[__gpu_level][2];
	uchar left[MAX_GRAPH_SIZE];
	uchar right[MAX_GRAPH_SIZE];
	
	Data() {}
	Data(uchar domains[][BDS], uchar cur[][2], uchar left[], uchar right[], uint bd_pos) {
		uint i;
		for (i = 0; i < BDS - 2; i++) {
			this->domains[i] = domains[bd_pos - 1][i];
		}
		for (i = 0; i < __gpu_level; i++) {
			this->cur[i][L] = cur[i][L];
			this->cur[i][R] = cur[i][R];
		}
		for (i = 0; i < n0; i++) {
			this->left[i] = left[i];
		}
		for (i = 0; i < n1; i++) {
			this->right[i] = right[i];
		}
	}
};
struct wrapperData  {
	uint i;
	uint number_of_bidomains;
	uint bound_size;
};

vector<Data> global_arguments;
vector<wrapperData> global_arguments_i;
vector<vector<uchar>> global_solution;

void *safe_realloc(void* old, uint new_size){
	void *tmp = realloc(old, new_size);
	if (tmp != NULL) return tmp;
	else exit(-1);
}

bool compareData(wrapperData d1, wrapperData d2)
{
	return (d1.bound_size <= d2.bound_size);
}

void clean_arguments_v2(vector<Data>& arguments_to_clean, vector<wrapperData>& arguments_to_clean_i, uchar inc_pos, uint& n_args, int& n_threads) {
	if(arguments.debug)
		cout << "Start cleaning ... " << n_threads << " thread already executed, " << +inc_pos << " size found" << endl;
		
	n_threads -= min(N_BLOCKS * BLOCK_SIZE, n_threads);
	if (n_threads == 0) {
		n_args = 0;
		arguments_to_clean.resize(0);
		arguments_to_clean_i.resize(n_threads);
		return;
	}
	uint lower_bound = 0;
	uint higher_bound = n_threads-1;
	uint i;
	uint insertion_pos_args_i = 0;
	
	vector<Data> args(0);
	args.reserve(arguments_to_clean.size());
	
	while (lower_bound != higher_bound) {
		i = (lower_bound + higher_bound) / 2;
		if (arguments_to_clean_i[i].bound_size <= inc_pos) {
			lower_bound = i + 1;
		}
		else {
			higher_bound = i;
		}
	}
	
	n_threads -= lower_bound;
	
	if (n_threads == 1 && arguments_to_clean_i[lower_bound].bound_size <= inc_pos) {
		n_threads = 0;
		n_args = 0;
		arguments_to_clean.resize(0);
		arguments_to_clean_i.resize(n_threads);
		return;
	}
	if ((lower_bound == 0) && ((time_to_clean++) != CLEAN)) {
		arguments_to_clean_i.resize(n_threads);
		return;
	}
	time_to_clean = 0;
	
	for (i = 0; i < n_threads; i++) {
		arguments_to_clean_i[insertion_pos_args_i] = arguments_to_clean_i[lower_bound + i];
		uint pos_args = args.size();
		for (int j = 0; j < arguments_to_clean_i[lower_bound + i].number_of_bidomains; j++) {
			args.push_back(arguments_to_clean[arguments_to_clean_i[lower_bound + i].i + j]);
		}
		arguments_to_clean_i[insertion_pos_args_i++].i = pos_args;
	}
	
	arguments_to_clean = args;
	n_args = args.size();
	if(arguments.debug)
		cout << "End cleaning" << endl;
	
	arguments_to_clean_i.resize(n_threads);
}


/*******************************************************************************
                                 MCS definitions
*******************************************************************************/

struct VtxPair {
    int v;
    int w;
    VtxPair(int v, int w): v(v), w(w) {}
};

struct Bidomain {
    int l,        r;        // start indices of left and right sets
    int left_len, right_len;
    bool is_adjacent;
    Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent):
            l(l),
            r(r),
            left_len (left_len),
            right_len (right_len),
            is_adjacent (is_adjacent) { };
};

struct AtomicIncumbent{
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
                        for (std::map<Position, Task>::iterator task = tasks.begin() ; task != tasks.end() ; ++task) {
                            // std::cout<< task->first.depth << " - ";
                             //for(int m = 0; m < task->first.values.size(); m++)
                                //std::cout << task->first.values[m] << " ";
                            //std::cout << std::endl;
                            if (task->second.func) { // whait for a function to be associated to this task by help_me_with()
                                auto f = task->second.func;
                                ++task->second.pending;
                                guard.unlock(); // so now other threads can stars executing this same function? why

                                auto start_work_time = steady_clock::now(); // local start time

                                (*f)(this_thread_nodes);

                                auto work_time = duration_cast<milliseconds>(steady_clock::now() - start_work_time);
                                total_work_time += work_time;

                                guard.lock();
                                task->second.func = nullptr; // why only now the function reference is removed?
                                if (0 == --task->second.pending) //decrement pending functions for this task
                                    cv.notify_all();            // notify all threads waiting on CV that this thread has finished his task

                                did_something = true;
                                break;
                            }
                        }
                        
                        if ((! did_something) && (! finish.load()))
                            cv.wait(guard); // if nothing has be done in the previous if scope, just wait for it
                    }

                    std::unique_lock<std::mutex> guard(general_mutex);
                    times.push_back(total_work_time);
                    nodes.push_back(this_thread_nodes);
                    });
    }
    auto kill_workers() -> void{
        {
            std::unique_lock<std::mutex> guard(general_mutex);
            finish.store(true);
            cv.notify_all();
        }

        for (std::thread & t : threads)
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
    ~HelpMe(){
        kill_workers();
    }
    HelpMe(const HelpMe &) = delete;
    void get_help_with(
            const Position & position,
            const std::function<void (unsigned long long &)> & main_func,
            const std::function<void (unsigned long long &)> & thread_func,
            unsigned long long & main_nodes){
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
std::atomic<int> indice_help_me(0);





/*******************************************************************************
                                 Merging functions
*******************************************************************************/

void prepara_argomenti_globali (AtomicIncumbent & global_incumbent, vector<VtxPair> & sol_corrente) { //////////////////////////////////////////////// NUOVE FUNZIONI MCSPLIT-MOSCA ////////////////////////////////////////////////////////////////////////

	uint sol_size = 1 + 2*(MIN(n0, n1));
	uint args_num = N_BLOCKS * BLOCK_SIZE * 2 * arguments.buffer_size;
	uint a_size = (BDS - 2 + 2 * __gpu_level + n0 + n1);
	int n_threads = global_arguments_i.size();
	uint n_args = global_arguments.size();
	uint arg_i = 0, i = 0;
	int min_size = MIN(g0->n, g1->n);
	
	uint args_i[N_BLOCKS * BLOCK_SIZE];
	
	uint args_size = args_num * a_size;
	uchar *args = (uchar*) malloc(args_size * sizeof *args);
	
	stable_sort(global_arguments_i.begin(), global_arguments_i.begin() + n_threads, compareData);
	
	uint dim_kernel = 0;
	for (uint b = 0; b < n_threads; b++)  {
		dim_kernel += global_arguments_i[n_threads-b-1].number_of_bidomains;
	}
	dim_kernel = dim_kernel * a_size;
    
	if (dim_kernel > args_size) {
		args = (uchar*) safe_realloc(args, dim_kernel * sizeof *args);
	}
	args_size = dim_kernel;
		
	for (uint b = 0; b < n_threads; b++) {
		args_i[b] = arg_i;
		for(uint c = 0; c < global_arguments_i[n_threads-b-1].number_of_bidomains; c++) {
			for (i = 0; i < BDS - 2; i++, arg_i++)
				args[arg_i] = global_arguments[global_arguments_i[n_threads-b-1].i + c].domains[i];
			for (i = 0; i < __gpu_level; i++, arg_i++)
				args[arg_i] = global_arguments[global_arguments_i[n_threads-b-1].i + c].cur[i][L];
			for (i = 0; i < __gpu_level; i++, arg_i++)
				args[arg_i] = global_arguments[global_arguments_i[n_threads-b-1].i + c].cur[i][R];
			for (i = 0; i < n0; i++, arg_i++)
				args[arg_i] = global_arguments[global_arguments_i[n_threads-b-1].i + c].left[i];
			for (i = 0; i < n1; i++, arg_i++)
				args[arg_i] = global_arguments[global_arguments_i[n_threads-b-1].i + c].right[i];
		}
	}
    
	if(arguments.debug)
		cout << "launching kernel:" << endl << 
		"\tthreads:              " << n_threads << endl <<
		"\targuments:            " << n_args << endl <<
		"\tthreads kernel:       " << N_BLOCKS * BLOCK_SIZE << endl <<
		"\targuments kernel:     " << arg_i/a_size << endl <<
		"\ttargs_size:           " << args_size << endl <<
		"\tbest bound:           " << global_arguments_i[n_threads-1].bound_size << endl <<
		"\tworst bound:          " << global_arguments_i[MAX(n_threads-1-N_BLOCKS*BLOCK_SIZE, 0)].bound_size<<endl<<
		"\tlast bound:           " << global_arguments_i[0].bound_size << endl;
	
	uchar inc_pos = global_incumbent.value;
    
	uchar local_solution[min_size][2];
	
	launch_kernel(args, n_threads, a_size, sol_size, args_i, local_solution, &inc_pos, args_size, arg_i);
	
	if(global_incumbent.update(inc_pos)) {
		sol_corrente.clear();
		for(int i=0; i<inc_pos; i++) {
			VtxPair pair ((int)local_solution[i][L], (int)local_solution[i][R]);
			sol_corrente.push_back(pair);
		}
	}
	
	if(arguments.debug)
		cout << "Soluzione migliore: " << +(uchar)inc_pos << endl;
	
	clean_arguments_v2(global_arguments, global_arguments_i, inc_pos, n_args, n_threads);
	
	if(arguments.debug)
		cout << "\tthreads:             " << n_threads << endl <<
		"\tdim args_i:          " << global_arguments_i.size() << endl <<
		"\targomenti:           " << n_args << endl <<
		"\tdim args:            " << global_arguments.size() << endl;
}

void compila_argomenti_globali (vector<Bidomain> & domains, uint bound_size, vector<VtxPair> & current_sol, vector<uchar> & left, vector<uchar> & right, AtomicIncumbent & global_incumbent, vector<VtxPair> & my_incumbent) { /////////////////////////////////////// NUOVA FUNZIONE ///////////////////////////////////////////
	Data new_arg;
	wrapperData pos_new_arg;
	pos_new_arg.i = global_arguments.size();
	pos_new_arg.number_of_bidomains = domains.size();
	pos_new_arg.bound_size = bound_size;
	global_arguments_i.push_back(pos_new_arg);
	for (int i=0; i<domains.size(); i++) {
		new_arg.domains[L] = domains[i].l;
		new_arg.domains[R] = domains[i].r;
		new_arg.domains[LL] = domains[i].left_len;
		new_arg.domains[RL] = domains[i].right_len;
		new_arg.domains[ADJ] = domains[i].is_adjacent;
		new_arg.domains[P] = split_levels+1;
		for (int j=0; j<__gpu_level; j++) {
			new_arg.cur[j][L] = current_sol[j].v;
			new_arg.cur[j][R] = current_sol[j].w;
		}
		for (int j=0; j<n0; j++) {
			new_arg.left[j] = left[j];
		}
		for (int j=0; j<n1; j++) {
			new_arg.right[j] = right[j];
		}
		global_arguments.push_back(new_arg);
	}
	if(global_arguments_i.size() == N_BLOCKS * BLOCK_SIZE * arguments.buffer_size){
        // we need to clean and discard all results that do not satisfy the bound
        // since the GPU has improved the solution
        // if after cleaning the size is not large enough we go back collecting new data
		prepara_argomenti_globali(global_incumbent, my_incumbent);
	}
	
}





/*******************************************************************************
                                 MCS functions
*******************************************************************************/



bool check_sol(const Graph & g0, const Graph & g1 , const vector<VtxPair> & solution) {
    vector<bool> used_left(g0.n, false);
    vector<bool> used_right(g1.n, false);
    unsigned int sol_size = solution.size();
    for (unsigned int i=0; i<sol_size; i++) {
        struct VtxPair p0 = solution[i];
        if (used_left[p0.v] || used_right[p0.w])
            return false;
        used_left[p0.v] = true;
        used_right[p0.w] = true;
        if (g0.label[p0.v] != g1.label[p0.w])
            return false;
        for (unsigned int j=i+1; j<sol_size; j++) {
            struct VtxPair p1 = solution[j];
            if (g0.adjmat[p0.v][p1.v] != g1.adjmat[p0.w][p1.w])
                return false;
        }
    }
    return true;
}

int calc_bound(const vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain &bd : domains) {
        bound += std::min(bd.left_len, bd.right_len);
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

int select_bidomain(const vector<Bidomain>& domains, const vector<int> & left,
        int current_matching_size){
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i=0; i<domains.size(); i++) {
        const Bidomain &bd = domains[i];
        if (arguments.connected && current_matching_size>0 && !bd.is_adjacent) continue;
        int len = arguments.heuristic == min_max ?
                std::max(bd.left_len, bd.right_len) :
                bd.left_len * bd.right_len;
        if (len < min_size) {
            min_size = len;
            min_tie_breaker = find_min_value(left, bd.l, bd.left_len);
            best = i;
        } else if (len == min_size) {
            int tie_breaker = find_min_value(left, bd.l, bd.left_len);
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

// multiway is for directed and/or labelled graphs
vector<Bidomain> filter_domains(const vector<Bidomain> & d, vector<int> & left,
        vector<int> & right, const Graph & g0, const Graph & g1, int v, int w,
        bool multiway){
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain &old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int left_len = partition(left, l, old_bd.left_len, g0.adjmat[v]);
        int right_len = partition(right, r, old_bd.right_len, g1.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;
        if (left_len_noedge && right_len_noedge)
            new_d.push_back({l+left_len, r+right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (multiway && left_len && right_len) {
            auto& adjrow_v = g0.adjmat[v];
            auto& adjrow_w = g1.adjmat[w];
            auto l_begin = std::begin(left) + l;
            auto r_begin = std::begin(right) + r;
            std::sort(l_begin, l_begin+left_len, [&](int a, int b)
                    { return adjrow_v[a] < adjrow_v[b]; });
            std::sort(r_begin, r_begin+right_len, [&](int a, int b)
                    { return adjrow_w[a] < adjrow_w[b]; });
            int l_top = l + left_len;
            int r_top = r + right_len;
            while (l<l_top && r<r_top) {
                unsigned int left_label = adjrow_v[left[l]];
                unsigned int right_label = adjrow_w[right[r]];
                if (left_label < right_label) {
                    l++;
                } else if (left_label > right_label) {
                    r++;
                } else {
                    int lmin = l;
                    int rmin = r;
                    do { l++; } while (l<l_top && adjrow_v[left[l]]==left_label);
                    do { r++; } while (r<r_top && adjrow_w[right[r]]==left_label);
                    new_d.push_back({lmin, rmin, l-lmin, r-rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});
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

void remove_vtx_from_left_domain(vector<int>& left, Bidomain& bd, int v){
    int i = 0;
    while(left[bd.l + i] != v) i++;
    std::swap(left[bd.l+i], left[bd.l+bd.left_len-1]);
    bd.left_len--;
}

void remove_bidomain(vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size()-1];
    domains.pop_back();
}

//////////////////////////////// GPU /////////////////////////////////////////




void GPU_solve(const unsigned depth, const Graph & g0, const Graph & g1,
                AtomicIncumbent & global_incumbent,
                PerThreadIncumbents & per_thread_incumbents,
                vector<VtxPair> & current, vector<Bidomain> & domains,
                vector<int> & left, vector<int> & right, const unsigned int matching_size_goal,
                const Position & position, unsigned long long & my_thread_nodes){

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

    int bd_idx = select_bidomain(domains, left, current.size());
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Bidomain &bd = domains[bd_idx];

    bd.right_len--;
    const int i_end = bd.right_len + 2; /* including the null */

///////////////////////////main function////////////////////
    int v = find_min_value(left, bd.l, bd.left_len);
    remove_vtx_from_left_domain(left, domains[bd_idx], v);
    int w = -1;

    for (int i = 0 ; i < i_end /* not != */ ; i++) {
        if (i != i_end - 1) {
            int idx = index_of_next_smallest(right, bd.r, bd.right_len+1, w);
            w = right[bd.r + idx];

            // swap w to the end of its colour class
            right[bd.r + idx] = right[bd.r + bd.right_len];
            right[bd.r + bd.right_len] = w;

            auto new_domains = filter_domains(domains, left, right, g0, g1, v, w,
                    arguments.directed || arguments.edge_labelled);
            current.push_back(VtxPair(v, w));
            int bound = current.size() + calc_bound(domains);
            if(bound > global_incumbent.value) {
                if (current.size() > split_levels) {
                    vector<uchar> char_left(left.size()), char_right(right.size());
                    for(int i = 0; i < left.size(); i++) {
                        char_left[i] = (uchar) left[i];
                    }
                    for(int i = 0; i < right.size(); i++) {
                        char_right[i] = (uchar) right[i];
                    }
                    compila_argomenti_globali (new_domains, bound, current, char_left, char_right, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second);
                }
                else {
                    auto new_position = position;
                    new_position.add(depth, i + 1);
                    GPU_solve(depth + 1, g0, g1, global_incumbent, per_thread_incumbents, current, new_domains, left, right, matching_size_goal, new_position, my_thread_nodes/*main_thread_nodes*/);
                }
            }
            current.pop_back();
        }
        else {
            // Last assign is null. Keep it in the loop to simplify parallelism.
            bd.right_len++;
            if (bd.left_len == 0)
                remove_bidomain(domains, bd_idx);

            int bound = current.size() + calc_bound(domains);
            if (bound > global_incumbent.value) {
                if (current.size() > split_levels) {
                    vector<uchar> char_left(left.size()), char_right(right.size());
                    for(int i = 0; i < left.size(); i++) {
                        char_left[i] = (uchar) left[i];
                    }
                    for(int i = 0; i < right.size(); i++) {
                        char_right[i] = (uchar) right[i];
                    }
                    compila_argomenti_globali (domains, bound, current, char_left, char_right, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second);
                }
                else {
                    auto new_position = position;
                    new_position.add(depth, i + 1);
                    GPU_solve(depth + 1, g0, g1, global_incumbent, per_thread_incumbents, current, domains, left, right, matching_size_goal, new_position, my_thread_nodes /*main_thread_nodes*/);
                }
            }
        }
    }
}


std::pair<vector<VtxPair>, unsigned long long> GPU_mcs(const Graph & g0, const Graph & g1) {
    vector<int> left;  // the buffer of vertex indices for the left partitions
    vector<int> right;  // the buffer of vertex indices for the right partitions

    srand(time(nullptr));
    auto domains = vector<Bidomain> {};

    std::set<unsigned int> left_labels;
    std::set<unsigned int> right_labels;
    for (unsigned int label : g0.label) left_labels.insert(label);
    for (unsigned int label : g1.label) right_labels.insert(label);
    std::set<unsigned int> labels;  // labels that appear in both graphs
    std::set_intersection(std::begin(left_labels),
                          std::end(left_labels),
                          std::begin(right_labels),
                          std::end(right_labels),
                          std::inserter(labels, std::begin(labels)));

    // Create a bidomain for each label that appears in both graphs
    for (unsigned int label : labels) {
        int start_l = left.size();
        int start_r = right.size();

        for (int i=0; i<g0.n; i++)
            if (g0.label[i]==label)
                left.push_back(i);
        for (int i=0; i<g1.n; i++)
            if (g1.label[i]==label)
                right.push_back(i);

        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        domains.push_back({start_l, start_r, left_len, right_len, false});
    }

    AtomicIncumbent global_incumbent;
    vector<VtxPair> incumbent;
    unsigned long long global_nodes = 0;

    if (arguments.big_first) {
        for (int k=0; k<g0.n; k++) {
            unsigned int goal = g0.n - k;
            auto left_copy = left;
            auto right_copy = right;
            auto domains_copy = domains;
            vector<VtxPair> current;
            PerThreadIncumbents per_thread_incumbents;
            per_thread_incumbents.emplace(std::this_thread::get_id(), vector<VtxPair>());
            Position position;
            
            GPU_solve(0, g0, g1, global_incumbent, per_thread_incumbents, current, domains_copy, left_copy, right_copy, goal, position, global_nodes);
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
        GPU_solve(0, g0, g1, global_incumbent, per_thread_incumbents, current, domains, left, right, 1, position, global_nodes);
        
        prepara_argomenti_globali(global_incumbent, incumbent);
        
        for (auto & i : per_thread_incumbents)
            if (i.second.size() > incumbent.size())
                incumbent = i.second;
    }

    return { incumbent, global_nodes };
}



//////////////////////////////////////////////////////////////////////




void solve_nopar(const unsigned depth, const Graph & g0, const Graph & g1,
        AtomicIncumbent & global_incumbent,
        vector<VtxPair> & my_incumbent,
        vector<VtxPair> & current, vector<Bidomain> & domains,
        vector<int> & left, vector<int> & right, const unsigned int matching_size_goal,
        unsigned long long & my_thread_nodes){
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

    int bd_idx = select_bidomain(domains, left, current.size());
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Bidomain &bd = domains[bd_idx];

    bd.right_len--;

    int v = find_min_value(left, bd.l, bd.left_len);
    remove_vtx_from_left_domain(left, domains[bd_idx], v);
    int w = -1;
    const int i_end = bd.right_len + 2; /* including the null */

    for (int i = 0 ; i < i_end /* not != */ ; i++) {
        if (i != i_end - 1) {
            int idx = index_of_next_smallest(right, bd.r, bd.right_len+1, w);
            w = right[bd.r + idx];

            // swap w to the end of its colour class
            right[bd.r + idx] = right[bd.r + bd.right_len];
            right[bd.r + bd.right_len] = w;

            auto new_domains = filter_domains(domains, left, right, g0, g1, v, w,
                    arguments.directed || arguments.edge_labelled);
            current.push_back(VtxPair(v, w));
            solve_nopar(depth + 1, g0, g1, global_incumbent, my_incumbent, current, new_domains, left, right, matching_size_goal, my_thread_nodes);
            current.pop_back();
        }
        else {
            // Last assign is null. Keep it in the loop to simplify parallelism.
            bd.right_len++;
            if (bd.left_len == 0)
                remove_bidomain(domains, bd_idx);

            solve_nopar(depth + 1, g0, g1, global_incumbent, my_incumbent, current, domains, left, right, matching_size_goal, my_thread_nodes);
        }
    }
}

void solve(const unsigned depth, const Graph & g0, const Graph & g1,
                AtomicIncumbent & global_incumbent,
                PerThreadIncumbents & per_thread_incumbents,
                vector<VtxPair> & current, vector<Bidomain> & domains,
                vector<int> & left, vector<int> & right, const unsigned int matching_size_goal,
                const Position & position, HelpMe & help_me, unsigned long long & my_thread_nodes){

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

    int bd_idx = select_bidomain(domains, left, current.size());
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Bidomain &bd = domains[bd_idx];

    bd.right_len--;
    std::atomic<int> shared_i{ 0 };
    const int i_end = bd.right_len + 2; /* including the null */

    // Version of the loop used by helpers
    std::function<void (unsigned long long &)> helper_function = [&shared_i, &g0, &g1, &global_incumbent, &per_thread_incumbents, &position, &depth,
        i_end, matching_size_goal, current, domains, left, right, &help_me] (unsigned long long & help_thread_nodes) {
        int which_i_should_i_run_next = shared_i++;

        if (which_i_should_i_run_next >= i_end)
            return; /* don't waste time recomputing */

        /* recalculate to this point */
        vector<VtxPair> help_current = current;
        vector<Bidomain> help_domains = domains;
        vector<int> help_left = left, help_right = right;

        /* rerun important stuff from before the loop */
        int help_bd_idx = select_bidomain(help_domains, help_left, help_current.size());
        if (help_bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
            return;
        Bidomain &help_bd = help_domains[help_bd_idx];

        int help_v = find_min_value(help_left, help_bd.l, help_bd.left_len);
        remove_vtx_from_left_domain(help_left, help_domains[help_bd_idx], help_v);

        int help_w = -1;

        for (int i = 0 ; i < i_end /* not != */ ; i++) {
            if (i != i_end - 1) {
                int idx = index_of_next_smallest(help_right, help_bd.r, help_bd.right_len+1, help_w);
                help_w = help_right[help_bd.r + idx];

                // swap w to the end of its colour class
                help_right[help_bd.r + idx] = help_right[help_bd.r + help_bd.right_len];
                help_right[help_bd.r + help_bd.right_len] = help_w;
                // don't perform recursions over vertices already considered by other threads
                if (i == which_i_should_i_run_next) {
                    which_i_should_i_run_next = shared_i++;
                    auto new_domains = filter_domains(help_domains, help_left, help_right, g0, g1, help_v, help_w,
                            arguments.directed || arguments.edge_labelled);
                    help_current.push_back(VtxPair(help_v, help_w));
                    if (depth > split_levels) {
                        solve_nopar(depth + 1, g0, g1, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, new_domains, help_left, help_right, matching_size_goal, help_thread_nodes);
                    }
                    else {
                        auto new_position = position;
                        new_position.add(depth, indice_help_me++);
                        solve(depth + 1, g0, g1, global_incumbent, per_thread_incumbents, help_current, new_domains, help_left, help_right, matching_size_goal, new_position, help_me, help_thread_nodes);
                    }
                    help_current.pop_back();
                }
            }
            else {
                // Last assign is null. Keep it in the loop to simplify parallelism.
                help_bd.right_len++;
                if (help_bd.left_len == 0)
                    remove_bidomain(help_domains, help_bd_idx);

                if (i == which_i_should_i_run_next) {
                    which_i_should_i_run_next = shared_i++;
                    if (depth > split_levels) {
                        solve_nopar(depth + 1, g0, g1, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, help_domains, help_left, help_right, matching_size_goal, help_thread_nodes);
                    }
                    else {
                        auto new_position = position;
                        new_position.add(depth, indice_help_me++);
                        solve(depth + 1, g0, g1, global_incumbent, per_thread_incumbents, help_current, help_domains, help_left, help_right, matching_size_goal, new_position, help_me, help_thread_nodes);
                    }
                }
            }
        }
    };

    // Grab this first, before advertising that we can get help
    int which_i_should_i_run_next = shared_i++;

    // Version of the loop used by the main thread
    std::function<void (unsigned long long &)> main_function = [&] (unsigned long long & main_thread_nodes) {
        int v = find_min_value(left, bd.l, bd.left_len);
        remove_vtx_from_left_domain(left, domains[bd_idx], v);
        int w = -1;
        //std::cout << "main funtion - " << std::this_thread::get_id() << std::endl;

        for (int i = 0 ; i < i_end /* not != */ ; i++) {
            if (i != i_end - 1) {
                int idx = index_of_next_smallest(right, bd.r, bd.right_len+1, w);
                w = right[bd.r + idx];

                // swap w to the end of its colour class
                right[bd.r + idx] = right[bd.r + bd.right_len];
                right[bd.r + bd.right_len] = w;

                if (i == which_i_should_i_run_next) {
                    which_i_should_i_run_next = shared_i++;
                    auto new_domains = filter_domains(domains, left, right, g0, g1, v, w,
                            arguments.directed || arguments.edge_labelled);
                    current.push_back(VtxPair(v, w));
                    if (depth > split_levels) {
                        solve_nopar(depth + 1, g0, g1, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, new_domains, left, right, matching_size_goal, main_thread_nodes);
                    }
                    else {
                        auto new_position = position;
                        new_position.add(depth, indice_help_me++);
                        solve(depth + 1, g0, g1, global_incumbent, per_thread_incumbents, current, new_domains, left, right, matching_size_goal, new_position, help_me, main_thread_nodes);
                    }
                    current.pop_back();
                }
            }
            else {
                // Last assign is null. Keep it in the loop to simplify parallelism.
                bd.right_len++;
                if (bd.left_len == 0)
                    remove_bidomain(domains, bd_idx);

                if (i == which_i_should_i_run_next) {
                    which_i_should_i_run_next = shared_i++;
                    if (depth > split_levels) {
                        solve_nopar(depth + 1, g0, g1, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, domains, left, right, matching_size_goal, main_thread_nodes);
                    }
                    else {
                        auto new_position = position;
                        new_position.add(depth, indice_help_me++);
                        solve(depth + 1, g0, g1, global_incumbent, per_thread_incumbents, current, domains, left, right, matching_size_goal, new_position, help_me, main_thread_nodes);
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

std::pair<vector<VtxPair>, unsigned long long> mcs(const Graph & g0, const Graph & g1, HelpMe & help_me) {
    vector<int> left;  // the buffer of vertex indices for the left partitions
    vector<int> right;  // the buffer of vertex indices for the right partitions

    srand(time(nullptr));
    auto domains = vector<Bidomain> {};

    std::set<unsigned int> left_labels;
    std::set<unsigned int> right_labels;
    for (unsigned int label : g0.label) left_labels.insert(label);
    for (unsigned int label : g1.label) right_labels.insert(label);
    std::set<unsigned int> labels;  // labels that appear in both graphs
    std::set_intersection(std::begin(left_labels),
                          std::end(left_labels),
                          std::begin(right_labels),
                          std::end(right_labels),
                          std::inserter(labels, std::begin(labels)));

    // Create a bidomain for each label that appears in both graphs
    for (unsigned int label : labels) {
        int start_l = left.size();
        int start_r = right.size();

        for (int i=0; i<g0.n; i++)
            if (g0.label[i]==label)
                left.push_back(i);
        for (int i=0; i<g1.n; i++)
            if (g1.label[i]==label)
                right.push_back(i);

        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        domains.push_back({start_l, start_r, left_len, right_len, false});
    }

    AtomicIncumbent global_incumbent;
    vector<VtxPair> incumbent;
    unsigned long long global_nodes = 0;

    if (arguments.big_first) {
        for (int k=0; k<g0.n; k++) {
            unsigned int goal = g0.n - k;
            auto left_copy = left;
            auto right_copy = right;
            auto domains_copy = domains;
            vector<VtxPair> current;
            PerThreadIncumbents per_thread_incumbents;
            per_thread_incumbents.emplace(std::this_thread::get_id(), vector<VtxPair>());
            Position position;
            position.add(0, indice_help_me++);
            for (auto & t : help_me.threads)
                per_thread_incumbents.emplace(t.get_id(), vector<VtxPair>());
            solve(0, g0, g1, global_incumbent, per_thread_incumbents, current, domains_copy, left_copy, right_copy, goal, position, help_me, global_nodes);
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
        position.add(0, indice_help_me++);
        for (auto & t : help_me.threads)
            per_thread_incumbents.emplace(t.get_id(), vector<VtxPair>());
        solve(0, g0, g1, global_incumbent, per_thread_incumbents, current, domains, left, right, 1, position, help_me, global_nodes);
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

GraphData write_Graph(GraphData* g0, GraphData* g1, vector<VtxPair>& solution) {

    unsigned int sol_size = solution.size();
    
    GraphData gd(Graph (sol_size), g0, g1);
    
    vector<bool> vtx_v(g0->g.n, false), vtx_w(g1->g.n, false);
    for (unsigned int i = 0; i < sol_size; i++) {
        vtx_v[solution[i].v] = true;
        vtx_w[solution[i].w] = true;
    }
    int ii = 0, jj = 0;
    for (int i = 0; i < g0->g.n; i++) {
        if (vtx_v[i]) {
            jj = 0;
            for (int j = 0; j < g0->g.n; j++) {
                if (vtx_v[j]) {
                    gd.g.adjmat[ii][jj] = g0->g.adjmat[i][j];
                    jj++;
                }
            }
            gd.g.label[ii] = g0->g.label[i];
            gd.map_g0[ii] = i;
            ii++;
        }
    }
    ii = 0;
    for (int i = 0; i < gd.g.n; i++) {
        for (unsigned int j = 0; j < sol_size; j++) {
            if (gd.map_g0.at(i) == solution.at(j).v) {
                gd.map_g1.at(i) = solution.at(j).w;
            }
        }
    }
    return gd;
}

void recursive_print (GraphData *gd, vector<int> &sol, int map) {
	if(gd->g0 != nullptr) {
		recursive_print(gd->g0, sol, gd->map_g0.at(map));
		recursive_print(gd->g1, sol, gd->map_g1.at(map));
	}
	else { // this is one of the original graphs
		sol.at(gd->ordine) = map;
	}
}

void nuova_print (vector<vector<GraphData>> &gd) {
	int n_files = arguments.filenames.size();
	vector<int> sol(n_files);
	GraphData *root = &gd.back().at(0);
	
	for (int i = 0; i < root->g.n; i++) {
        if (root->map_g0.at(i) > 20) 
            cout << "map_g0: " << root->map_g0.at(i) << endl;
		recursive_print(root->g0, sol, root->map_g0.at(i));
        if (root->map_g1.at(i) > 20)
            cout << "map_g1: " << root->map_g1.at(i) << endl;
		recursive_print(root->g1, sol, root->map_g1.at(i));
		
		cout << sol.at(0);
		for (int j = 1; j < n_files; j++) {
			cout << " -> " << sol.at(j);
		}
		cout << endl;
	}
}

////////////////////////////// GPU ////////////////////////////////////////



void GPU_produce_solution (vector<GraphData> &graphs, vector<GraphData> &sol, int index) {

	Graph *g0 = &graphs.at(index).g;
	Graph *g1 = &graphs.at(graphs.size()-1-index).g;
	vector<int> g0_deg = calculate_degrees(*g0);
	vector<int> g1_deg = calculate_degrees(*g1);

	// As implemented here, g1_dense and g0_dense are false for all instances
	// in the Experimental Evaluation section of the paper.  Thus,
	// we always sort the vertices in descending order of degree (or total degree,
	// in the case of directed graphs.  Improvements could be made here: it would
	// be nice if the program explored exactly the same search tree if both
	// input graphs were complemented.
	vector<int> vv0(g0->n);
	std::iota(std::begin(vv0), std::end(vv0), 0);
	bool g1_dense = sum(g1_deg) > g1->n*(g1->n-1);
	std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
		return g1_dense ? (g0_deg[a]<g0_deg[b]) : (g0_deg[a]>g0_deg[b]);
	});
	vector<int> vv1(g1->n);
	std::iota(std::begin(vv1), std::end(vv1), 0);
	bool g0_dense = sum(g0_deg) > g0->n*(g0->n-1);
	std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) {
		return g0_dense ? (g1_deg[a]<g1_deg[b]) : (g1_deg[a]>g1_deg[b]);
	});

	struct Graph g0_sorted = induced_subgraph(*g0, vv0);
	struct Graph g1_sorted = induced_subgraph(*g1, vv1);
	
	n0 = g0->n;
	n1 = g1->n;
	
	for (int i = 0; i < n0; i++)
        for (int j = 0; j < n0; j++)
            adjmat0[i][j] = g0_sorted.adjmat[i][j];

    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n1; j++)
            adjmat1[i][j] = g1_sorted.adjmat[i][j];
    
    checkCudaErrors(cudaDeviceReset());
    
    move_graphs_to_gpu(&g0_sorted, &g1_sorted);
	
	std::pair<vector<VtxPair>, unsigned long long> solution = GPU_mcs(g0_sorted, g1_sorted);

	// Convert to indices from original, unsorted graphs
	for (auto& vtx_pair : solution.first) {
		vtx_pair.v = vv0[vtx_pair.v];
		vtx_pair.w = vv1[vtx_pair.w];
	}

    if (!check_sol(*g0, *g1, solution.first)) {
        cout << "WRONG SOLUTION!" << endl;
        exit (-1);
    }

	sol.at(index) = write_Graph(&graphs.at(index), &graphs.at(graphs.size()-1-index), solution.first);
}


///////////////////////////////////////////////////////////////////////////////

void produce_solution (vector<GraphData> &graphs, vector<GraphData> &sol, int index, HelpMe & help_me) {
	Graph *g0 = &graphs.at(index).g;
	Graph *g1 = &graphs.at(graphs.size()-1-index).g;
	vector<int> g0_deg = calculate_degrees(*g0);
	vector<int> g1_deg = calculate_degrees(*g1);

	// As implemented here, g1_dense and g0_dense are false for all instances
	// in the Experimental Evaluation section of the paper.  Thus,
	// we always sort the vertices in descending order of degree (or total degree,
	// in the case of directed graphs.  Improvements could be made here: it would
	// be nice if the program explored exactly the same search tree if both
	// input graphs were complemented.
	vector<int> vv0(g0->n);
	std::iota(std::begin(vv0), std::end(vv0), 0);
	
	bool g1_dense = sum(g1_deg) > g1->n*(g1->n-1);
	
	std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
		return g1_dense ? (g0_deg[a]<g0_deg[b]) : (g0_deg[a]>g0_deg[b]);
	});
	
	vector<int> vv1(g1->n);
	std::iota(std::begin(vv1), std::end(vv1), 0);
	
	bool g0_dense = sum(g0_deg) > g0->n*(g0->n-1);
	
	std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) {
		return g0_dense ? (g1_deg[a]<g1_deg[b]) : (g1_deg[a]>g1_deg[b]);
	});

	struct Graph g0_sorted = induced_subgraph(*g0, vv0);
	struct Graph g1_sorted = induced_subgraph(*g1, vv1);

	std::pair<vector<VtxPair>, unsigned long long> solution = mcs(g0_sorted, g1_sorted, help_me);

	// Convert to indices from original, unsorted graphs
	for (auto& vtx_pair : solution.first) {
		vtx_pair.v = vv0[vtx_pair.v];
		vtx_pair.w = vv1[vtx_pair.w];
	}

    if (!check_sol(*g0, *g1, solution.first)) {
        cout << "WRONG SOLUTION!" << endl;
        exit (-1);
    }
	sol.at(index) = write_Graph(&graphs.at(index), &graphs.at(graphs.size()-1-index), solution.first);
}

void sort_by_size_ascending(std::vector <struct GraphData> &gi) {
    int n_graphs = gi.size();
    for (int i = 0; i < n_graphs - 1; i++) {
        for (int j = i + 1; j < n_graphs; j++) {
            if (gi.at(i).g.n > gi.at(j).g.n) {
                GraphData g = gi.at(i);
                gi.at(i) = gi.at(j);
                gi.at(j) = g;
            }
        }
    }
}

auto floatToDuration(const float time_s) {
    using namespace std::chrono;
    using fsec = duration<float>;
    return round<nanoseconds>(fsec{time_s});
}

int main(int argc, char** argv) {
    set_default_arguments();
    argp_parse(&argp, argc, argv, 0, 0, 0);
    
    int n_files = arguments.filenames.size();
    int recursion_number = ceil(log2(n_files)) + 1;
    vector<vector<GraphData>> gi_data(recursion_number);


    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';
    for (int i=0; i<n_files; i++) {
        gi_data.at(0).push_back(readGraph(arguments.filenames.at(i), format, arguments.directed,
            arguments.edge_labelled, arguments.vertex_labelled));
        gi_data.at(0).back().ordine = i;
    }

    
	struct timespec s, finish;
	float time_elapsed;
	clock_gettime(CLOCK_MONOTONIC, &s);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    HelpMe help_me(arguments.threads - 1);
    vector<std::thread> t; 
    bool aborted = false;
    bool aborted_at_least_once = false;
    float timeout = arguments.timeout;
    
    for (int j = 0; j < recursion_number-1 && !aborted; j++) {
        std::thread timeout_thread;
        std::mutex timeout_mutex;
        std::condition_variable timeout_cv;
        abort_due_to_timeout.store(false);

        if (0 != arguments.timeout) {
            timeout_thread = std::thread([&] {
                if (j != recursion_number-2) {
                    timeout = timeout/2;
                }
                auto abort_time = steady_clock::now() + floatToDuration(timeout);
                {
                    /* Sleep until either we've reached the time limit,
                        * or we've finished all the work. */
                    std::unique_lock<std::mutex> guard(timeout_mutex);
                    while (! abort_due_to_timeout.load()) {
                        if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                            /* We've woken up, and it's due to a timeout. */
                            aborted_at_least_once = true;
                            aborted = true;
                            break;
                        }
                    }
                }
                abort_due_to_timeout.store(true);
                });
        }
        sort_by_size_ascending(gi_data.at(j));
        gi_data.at(j+1).resize(gi_data.at(j).size()/2 + gi_data.at(j).size()%2);
        
        for (unsigned int i = 0; i < (unsigned int)(gi_data.at(j).size()/2); i++) {
            if (gi_data.at(j).size()/2 > 1 && i == 0) {
                t.emplace_back(std::thread([&gi_data, i, j] { GPU_produce_solution(gi_data.at(j), gi_data.at(j+1), i); } ));
            }
            else {
                t.emplace_back(std::thread([&gi_data, i, j, &help_me] { produce_solution(gi_data.at(j), gi_data.at(j+1), i, help_me); } ));
            }
        }
        if (gi_data.at(j).size() % 2) {
            gi_data.at(j+1).at(gi_data.at(j).size()/2) = gi_data.at(j).at(gi_data.at(j).size()/2);
        }
        for (unsigned int i = 0; i < (unsigned int)(gi_data.at(j).size()/2); i++) {
            t.at(i).join();
        }
        t.clear();

        if (arguments.verbose) {
            cout << "Iteration " << j << " solutions: ";
            for (auto &data : gi_data.at(j+1)) {
                cout << data.g.n << " ";
            }
            if (aborted) cout << "TIMEOUT";
            cout << endl;
        }

        // Clean up the timeout thread
        if (timeout_thread.joinable()) {
            {
                std::unique_lock<std::mutex> guard(timeout_mutex);
                abort_due_to_timeout.store(true);
                timeout_cv.notify_all();
            }
            timeout_thread.join();
            aborted = false;
        }
    }
    
    help_me.kill_workers();
    
    clock_gettime(CLOCK_MONOTONIC, &finish);
    time_elapsed = (finish.tv_sec - s.tv_sec);
    time_elapsed += (finish.tv_nsec - s.tv_nsec) / 1000000000.0;

    nuova_print(gi_data);
    int sol_size = gi_data.back().at(0).g.n;
    cout << ">>> " << sol_size << " - " << (double)time_elapsed << endl;

    if (aborted_at_least_once)
        cout << "TIMEOUT" << endl;
    
    return 0; /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
}
