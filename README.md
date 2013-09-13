###pareto.py
####Nondominated sorting for multi-objective problems

Sorts one or more files of solutions into the Pareto-efficient (or "nondominated") set. Requires [NumPy](http://www.numpy.org/). Files can contain columns other than objectives which are not sorted.

Example usage for a group of input files beginning with `my_file_`:
```
python pareto.py \
	   --input my_file_* \
	   --output my_pareto_set.txt \
	   --objectives 3 5 7 \
	   --epsilons 0.01 0.05 0.1 \
	   --delimiter=' ' \
	   --precision=8 \
	   --print-only-objectives
```

* `-i, --input`: Required. List of input files to sort, separated by spaces. Input files are assumed to contain floating-point values separated by `delimiter`. Lines beginning with `#` are treated as comments. Input files must all contain the same number of columns.

* `-o, --output`: Required. Filename to output your Pareto set.

* `--objectives`: Optional. A list of columns of the input files to sort (zero-indexed), separated by spaces. If not given, all columns of the input files will be sorted by default.

* `--epsilons`: Optional. A list of epsilon (precision) values corresponding to each objective. If not given, all objectives will use a precision of `1e-9`. 
 
* `--delimiter`: Input file delimiter (optional). Common choices:
	* Space-delimited (default): `--delimiter=' '`
	* Comma-delimited: `--delimiter=','`
	* Tab-delimited: `--delimiter=$'\t'`

* `--precision`: Digits of precision in the output file (optional). Default is 8.

* `--print-only-objectives`: Optional. Include this flag to print only the objective values in the Pareto set. If this flag is not included, by default all columns of the input will be printed to the output, even if they were not sorted.

### License
Copyright (C) 2013 Jon Herman, Patrick Reed and others. Licensed under the GNU Lesser General Public License.

The Sensitivity Analysis Library is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The Sensitivity Analysis Library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the Sensitivity Analysis Library.  If not, see <http://www.gnu.org/licenses/>.
