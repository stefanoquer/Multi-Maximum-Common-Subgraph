# Usage
mcsp [OPTIONS] [files...]

Positionals:  
&ensp; files TEXT:FILE ...         At least 2 input files

Options:
| Options:                  | Descriptions:                                                                     |
|---------------------------|-----------------------------------------------------------------------------------|
| -h,--help                 | Print this help message and exit                                                  |
| -q,--quiet                | Quiet output                                                                      |
| -v,--verbose              | Verbose output                                                                    |
| -c,--connected            | Solve max common CONNECTED subgraph problem                                       |
| -i,--directed             | Use directed graphs                                                               |
| -a,--labelled             | Use edge and vertex labels                                                        |
| -x,--vertex-labelled-only | Use vertex labels, but not edge labels                                            |
| -b,--big-first            | First try to find an induced subgraph isomorphism, then decrement the target size |
| -t,--timeout INT          | Specify a timeout (seconds)                                                       |
| -T,--thread INT           | Specify how many threads to use                                                   |

## [Option Group: Format]
   REQUIRED 
  Options:
| Options:    | Descriptions:                |
|-------------|------------------------------|
| -d,--dimacs | Read DIMACS format           |
| -l,--ladder | Read LAD format              |
| -b,--bin    | Read binary format           |
| -e,--bin_e  | Read binary alternate format |
| -o,--ioi    | Read IOI format              |

## [Option Group: Heuristic]
   REQUIRED 
| Options:      | Descriptions:             |
|---------------|---------------------------|
| --min_max     | Use min_max heuristic     |
| --min_min     | Use min_min heuristic     |
| --min_sum     | Use min_sum heuristic     |
| --min_product | Use min_product heuristic |