###pareto.py
####Nondominated sorting for multi-objective problems
by [jdherman](https://github.com/jdherman) and [matthewjwoodruff](https://github.com/matthewjwoodruff)

An epsilon-nondominated sort.  Sorts one or more files of solutions into the Pareto-efficient (or "nondominated") set.  Solutions can contain columns other than objectives, which will be carried through, unsorted, to the output.  By default, output rows are reproduced verbatim from input.  `pareto.py` assumes that all objectives are to be minimized.

This sort assumes a desired output resolution (epsilons).  If a strict nondominated sort is required, it can be approximated by setting epsilons arbitrarily small (within reason -- floating-point division is involved here.)  The default epsilon resolution of 1e-9 will effectively result in a strict nondominated sort in many cases.

Example usage for a group of input files beginning with `my_file_`:
```
python pareto.py \
	   --input my_file_* \
	   --output my_pareto_set.txt \
	   --objectives 3 5 7 \
	   --epsilons 0.01 0.05 0.1 \
	   --delimiter=' ' \
	   --print-only-objectives \
           --header=1 \
           --blank \
           --comment="#" \
           --contribution \
           --line-number
```

* `-i, --input`: Required. List of input files to sort, separated by spaces. Input files are assumed to contain floating-point values separated by `delimiter`. Input files must all contain the same number of columns. 

* `-o, --output`: Required. Filename to output your Pareto set.

* `--objectives`: Optional. A list of columns of the input files to sort (zero-indexed), separated by spaces. If not given, all columns of the input files will be sorted.

* `--epsilons`: Optional. A list of epsilon (precision) values corresponding to each objective. If not given, all objectives will use a precision of `1e-9`. 
 
* `--delimiter`: Input file delimiter (optional). Common choices:
	* Space-delimited (default): `--delimiter=' '`
	* Comma-delimited: `--delimiter=','`
	* Tab-delimited: `--delimiter=$'\t'`

* `--print-only-objectives`: Optional. Include this flag to print only the objective values in the Pareto set. If this flag is not included, all columns of the input will be printed to the output, even if they were not sorted.

* `--header`: Optional. Number of rows to skip at the top of each input.  If the header is prefixed with a comment character (see `--comment`) this is unnecessary.

* `--comment`: Optional.  Character to look for at the beginning of a row indicating that the row contains a comment and should be ignored.  Commented rows are not preserved in output.

* `--blank`: Optional. Skip blank rows instead of erroring out.

* `--contribution`:  Optional.  Append filename of origin to each solution in the output.

* `--line-number`: Optional.  Report line number in file of origin when appending contribution.


### What is this?
For more information, please consult the following references:

* <https://en.wikipedia.org/wiki/Pareto_efficiency>

* <https://en.wikipedia.org/wiki/Multi-objective_optimization>

* <http://www.iitk.ac.in/kangal/index.shtml>: Kanpur Genetic Algorithms Laboratory Homepage

* <http://waterprogramming.wordpress.com>: Waterprogramming Research Blog

* Deb, K., M. Mohan, and S. Mishra. 2005 "Evaluating the epsilon-domination based multiobjective evolutionary algorithm for a quick computation of Pareto-optimal solutions"  *Evolutionary Computation Journal* 13 (4): 501-525.

* Deb, K., A. Pratap, S. Agarwal, and T. Meyarivan. 2002. "A fast and elitist multiobjective genetic algorithm: NSGA-II." *IEEE Transactions on Evolutionary Computation* 6 (2): 182-197.

* Srinivas, N., and K. Deb. 1994. "Multiobjective optimization using nondominated sorting in genetic algorithms." *Evolutionary Computation* 2 (2): 221-248.

### But what is it for?
Among other things, this script can serve the following purposes:

* `pareto.py` may be used to post-process the output of multiple optimization runs.  Since multiple-objective evolutionary algorithms (MOEAs) are search heuristics using random numbers, they may be run more than once and produce different output each time.  Therefore the estabilshed best practice for MOEA users is to sort together the output of as many optimization runs as possible to get a good approximation of the true Pareto-efficient set.
* `pareto.py` may be used to combine the output of more than one MOEA, for research studies comparing the performance of different MOEAs.  This is useful when computing metrics of optimization performance such as hypervolume and the epsilon-indicator.
* `pareto.py` may be used to find the Pareto-efficient solutions within a set of data that does not result from optimization.
* `pareto.py` may be used to prepare reference data for validating other implementations of epsilon-nondominated sorting routines.

### Alternatives
The `ReferenceSetMerger` from moeaframework (<http://www.moeaframework.org>) is probably faster in most cases, since it is compiled Java rather than interpreted Python.  It has fewer options.  YMMV, TANSTAAFL.

### Dependencies

* Python standard library (`sys`, `math`, and `argparse`)
* Python 2.7 or later, Python 3.2 or later (for `argparse`)

### Note for Pandas users
Pandas is an excellent library for Python data analysis.  Doing a nondominated 
sort on a Pandas Data Frame requires `itertuples(False)`, because `pareto.py` expects 
to be able to iterate over rows.  The `False` argument excludes the index.  It can
be included if desired, but the objective columns should be adjusted accordingly.

```
import pandas
import pareto

table = pandas.read_table("datafile.txt")
nondominated = pareto.eps_sort([table.itertuples(False)], [3, 4, 5], [1, 0.1, 3])
```

### License
Copyright (C) 2013 Jon Herman, Matt Woodruff, Patrick Reed and others. 
Licensed under the GNU Lesser General Public License.

pareto.py is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pareto.py is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the pareto.py.  If not, see <http://www.gnu.org/licenses/>.
