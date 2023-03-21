# Usage
Usage: mcsp [OPTION...] HEURISTIC FILENAME1 FILENAME2

Find a maximum clique in a graph in DIMACS format

| Option                      | Description                                                                       |
|-----------------------------|-----------------------------------------------------------------------------------|
| -a, --labelled              | Use edge and vertex labels                                                        |
| -b, --big-first             | First try to find an induced subgraph isomorphism, then decrement the target size |
| -c, --connected             | Solve max common CONNECTED subgraph problem                                       |
| -d, --dimacs                | Read DIMACS format                                                                |
| -i, --directed              | Use directed graphs                                                               |
| -l, --lad                   | Read LAD format                                                                   |
| -q, --quiet                 | Quiet output                                                                      |
| -t, --timeout=timeout       | Specify a timeout (seconds)                                                       |
| -T, --threads=threads       | Specify how many threads to use                                                   |
| -v, --verbose               | Verbose output                                                                    |
| -x, --vertex-labelled-only  | Use vertex labels, but not edge labels                                            |
| -?, --help                  | Give this help list                                                               |
| --usage                     | Give a short usage message                                                        |

Mandatory or optional arguments to long options are also mandatory or optional for any corresponding short options.

HEURISTIC can be min_max or min_product