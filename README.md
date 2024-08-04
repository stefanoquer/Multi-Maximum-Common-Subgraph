# McSplit Multigraph

Please refer to the following paper for a complete description of the tools:

L. Cardone, S. Quer, "The Multi-Maximum and Quasi-Maximum Common Subgraph Problem," Computation, Publisher MDPI AG, Switzerland, Vol. 11., No. 4, 2023, pp. 69-94, https://www.mdpi.com/2079-3197/11/4/69, ISSN: 2079-3197,
DOI: 10.3390/computation11040069

The repository includes the following tools

## [Quer](https://github.com/Lorenzo-Cardone/McSplit_Multigraph/tree/Quer)
Custom implementation, little to no optimizations, only sequential, written in C.

## [Quer_Parallel](https://github.com/Lorenzo-Cardone/McSplit_Multigraph/tree/Quer_Parallel)
Same as the previous, but ported to C++ and made parallel in order to allow for multithreading.

## [MultiGraph](https://github.com/Lorenzo-Cardone/McSplit_Multigraph/tree/MultiGraph)
The default implementation, derived from the original McSplit code with added multithreading, with very little optimization applied.

## [MultiGraph_FixedSize](https://github.com/Lorenzo-Cardone/McSplit_Multigraph/tree/MultiGraph_FixedSize)
Improved implementation: dynamic sized structures are mostly updated to static ones.

## [Sorted_FixedSize](https://github.com/Lorenzo-Cardone/McSplit_Multigraph/tree/Sorted_FixedSize)
Improved implementation: on top of fixed size improvements, we have added a domain sorting heuristic.

## [Sorted_FixedSize_PrintAllMCS](https://github.com/Lorenzo-Cardone/McSplit_Multigraph/tree/Sorted_FixedSize_PrintAllMCS)
Debug version only added additional prints.

## [Tree](https://github.com/Lorenzo-Cardone/McSplit_Multigraph/tree/Tree)
Approximated multiple graphs solution. Solves graph pairs in parallel then solves the returned solutions etc.

## [Tree-GPU](https://github.com/Lorenzo-Cardone/McSplit_Multigraph/tree/Tree-GPU)
Approximated multiple graphs solution. Same as previous, but adds a CUDA-enabled device as a supporting computing unit.

## [Waterfall](https://github.com/Lorenzo-Cardone/McSplit_Multigraph/tree/Waterfall)
Approximated multiple graphs solution. Solves only one pair of graphs at a time, then repeats the process between the solution and a new graph.

# Benchmarks
All tests have been conducted on the [ARG Database](https://mivia.unisa.it/datasets/graph-database/arg-database/), which is publicly available.
