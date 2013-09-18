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
	   --print-only-objectives
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

### What is this?
For more information, please consult the following references:

* [https://en.wikipedia.org/wiki/Pareto_efficiency](Wikipedia: Pareto Efficiency)

* [https://en.wikipedia.org/wiki/Multi-objective_optimization](Wikipedia: Multi-objective Optimization)

* [http://www.iitk.ac.in/kangal/index.shtml](Kanpur Genetic Algorithms Laboratory Homepage)

* [http://waterprogramming.wordpress.com](Waterprogramming Research Blog)

* Deb, K., M. Mohan, and S. Mishra. 2005 "Evaluating the epsilon-domination based multiobjective evolutionary algorithm for a quick computation of Pareto-optimal solutions"  *Evolutionary Computation Journal* 13 (4): 501-525.

* Deb, K., A. Pratap, S. Agarwal, and T. Meyarivan. 2002. "A fast and elitist multiobjective genetic algorithm: NSGA-II." *IEEE Transactions on Evolutionary Computation* 6 (2): 182-197.

* Srinivas, N., and K. Deb. 1994. "Multiobjective optimization using nondominated sorting in genetic algorithms." *Evolutionary Computation* 2 (2): 221-248.

### Alternatives
The sort in (moeaframework)[http://www.moeaframework.org] is probably faster in most cases, since it is compiled Java rather than interpreted Python.  It has fewer options.  YMMV, TANSTAAFL.

### Dependencies

* Python standard library (`sys`, `math`, and `argparse`)
* Python 2.7 or later, Python 3.2 or later (for `argparse`)

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
