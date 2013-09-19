"""
Copyright (C) 2013 Jon Herman, Matthew Woodruff, Patrick Reed, and others

This script is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This script is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this script. If not, see <http://www.gnu.org/licenses/>.
===========================================================
pareto.py

Perform epsilon-nondominated sort on input files.

Please cite the following works if publishing results obtained using this
script.

For pareto.py:

@misc{herman_woodruff_2013_pareto,
    author = {Herman, Jon and Woodruff, Matthew},
    year = {2013},
    title = {pareto.py: a $\\varepsilon-nondomination$ sorting routine},
    howpublished = {https://github.com/jdherman/pareto.py}
}

For epsilon-nondomination:
@article{deb_2005_emoea,
    author = { Deb, K. and Mohan, M. and Mishra, S},
    year = {2005},
    title = {Evaluating the $\\varepsilon$-domination based
        multiobjective evolutionary algorithm for a quick
        computation of Pareto-optimal solutions.},
    journal = {Evolutionary Computation Journal},
    volume= {13}, number = {4}, pages ={501--525}
}

For a fast nondominated sort:
@article{deb_2002_nsga2,
    title="A fast and elitist multiobjective genetic algorithm: {NSGA-II}",
    author="Deb, Kalyanmoy and Pratap, A and Agarwal, S and Meyarivan, T",
    volume="6", number="2",
    journal="{IEEE} Transactions on Evolutionary Computation",
    year="2002",
    pages="182--197"
}
"""

import sys
import math
import argparse

class SortParameterError(Exception): pass

class Archive(object):
    """ An archive of epsilon-nondominated solutions """
    def __init__(self, epsilons, oindices):
        """
        epsilons: sizes of epsilon boxes to use in the sort
        oindices: indicate which indices in a solution are objectives
        """
        self.archive = []
        self.boxes = [] # remember boxes
        self.objectives = [] # remember objectives
        self.epsilons = epsilons
        self.sortinto = self._realsortinto
        if oindices is not None:
            self.oindices = oindices
            self.nobj = len(oindices)
            self.itobj = range(self.nobj)
        else:
            self.oindices = []
            self.nobj = 0
            self.itobj = range(0)
            self.sortinto = self._initialsortinto

    def add(self, solution, sobj, sbox):
        """ add a solution to the archive, plus auxiliary information """
        self.archive.append(solution)
        self.objectives.append(sobj)
        self.boxes.append(sbox)

    def remove(self, index):
        """ remove a solution from the archive """
        self.archive.pop(index)
        self.objectives.pop(index)
        self.boxes.pop(index)

    def _initialsortinto(self, solution):
        """
        Gets called the very first time, to establish the 
        number of objectives and, if not supplied, epsilons.
        """
        self.nobj = len(solution)
        self.oindices = range(self.nobj)
        self.itobj = range(self.nobj)

        if self.epsilons is None:
            self.epsilons = [1e-9]*self.nobj
        elif len(self.epsilons) != self.nobj:
            msg = "{0} epsilons specified, but found {1} columns".format(
                    len(self.epsilons), self.nobj)
            raise SortParameterError(msg)

        self.sortinto = self._realsortinto
        self.sortinto(solution)

    def _realsortinto(self, solution):
        """
        Sort a solution into the archive.  Add it if it's nondominated
        w.r.t current solutions.
        """
        # Here's how the early loop exits in this code work:
        # break:    Stop iterating the box comparison for loop because we know
        #           the solutions are in relatively nondominated boxes.
        # continue: Start the next while loop iteration immediately (i.e.
        #           jump ahead to the comparison with the next archive member).
        # return:   The candidate solution is dominated, stop comparing it to
        #           the archive, don't add it, immediately exit the method.

        sobj = [float(solution[ii]) for ii in self.oindices]
        sbox = [math.floor(sobj[ii] / self.epsilons[ii]) for ii in self.itobj]

        asize = len(self.archive)

        ai = -1
        while ai < asize - 1:
            ai += 1
            adominate = False # archive dominates
            sdominate = False # solution dominates
            nondominate = False # neither dominates

            abox = self.boxes[ai]

            for oo in self.itobj:
                if abox[oo] < sbox[oo]:
                    adominate = True
                    if sdominate: # nondomination
                        nondominate = True
                        break # for
                elif abox[oo] > sbox[oo]:
                    sdominate = True
                    if adominate: # nondomination
                        nondominate = True
                        break # for

            if nondominate:
                continue # while
            if adominate: # candidate solution was dominated
                return
            if sdominate: # candidate solution dominated archive solution
                self.remove(ai)
                ai -= 1
                asize -= 1
                continue # while

            # solutions are in the same box
            aobj = self.objectives[ai]
            corner = [sbox[ii] * self.epsilons[ii] for ii in self.itobj]
            sdist = sum([(sobj[ii] - corner[ii]) **2 for ii in self.itobj])
            adist = sum([(aobj[ii] - corner[ii]) **2 for ii in self.itobj])
            if adist < sdist: # archive dominates
                return
            else: # solution dominates
                self.remove(ai)
                ai -= 1
                asize -= 1
                # Need a continue here if we ever reorder the while loop.
                continue # while

        # if you get here, then no archive solution has dominated this one
        self.add(solution, sobj, sbox)

def get_args(argv):
    """ Get command line arguments """
    prog = argv.pop(0)
    parser = argparse.ArgumentParser(prog=prog,
        description='Nondomination Sort for Multiple Files')
    parser.add_argument('--objectives', type=int, nargs='+', required=False,
                        help='Objective Columns (zero-indexed)')
    parser.add_argument('--epsilons', type=float, nargs='+', required=False,
                        help='Epsilons, one per objective')
    parser.add_argument('-o', '--output', type=str,
                        required=True, help='Output Filename')
    parser.add_argument('-i', '--input', type=str, required=True,
                        nargs='+', help='Input filenames')
    parser.add_argument('--delimiter', type=str, required=False,
                        default=' ', help='Input column delimiter')
    parser.add_argument('--print-only-objectives', action='store_true',
                        default=False, required=False,
                        help='Print only objectives in output')
    parser.add_argument("--blank", action="store_true",
                        help="skip blank lines")
    parser.add_argument("--comment", type=str,
                        help="skip lines starting with this character")
    parser.add_argument("--header", type=int,
                        help="number of header lines to skip")
    return parser.parse_args(argv)

class SortInputError(Exception):
    """ Information about a defective input """
    def __init__(self, msg, row, table):
        super(SortInputError, self).__init__(msg)
        self.row = row
        self.table = table

def eps_sort(tables, objectives, epsilons):
    """
    Perform an epsilon-nondominated sort
    tables: input data, must support iteration
    objectives: list of column indices in which objectives can be found,
                if None default to all columns
    epsilons: list of epsilons for the sort, if None default to 1e-9
    """
    archive = Archive(epsilons, objectives)

    # for each file in argument list ...
    for counter in range(len(tables)):
        solutions = tables[counter]
        # for each line in file (new candidate solution) ...
        rownumber = 0
        for row in solutions:
            try:
                archive.sortinto(row)
                rownumber += 1
            except IndexError:
                msg = "Not enough columns in row {0} of input {1}".format(
                                               rownumber, counter)
                raise SortInputError(msg, rownumber, counter)
            except ValueError as ve:
                msg = "{0} on row {1} of input {2}".format(
                                                ve.message, rownumber, counter)
                raise SortInputError(msg, rownumber, counter)

    return archive

def rowsof(filename, delimiter):
    """ 
    Generator function yielding rows read from a file. (Lazy input.)
    Avoids having to read the whole file at once.
    """
    with open(filename, 'r') as fp:
        try:
            while True:
                line = fp.next()
                row = line.strip().split(delimiter)
                yield row
        except StopIteration:
            pass

def filter_input(rows, **kwargs):
    """
    Generator function filtering out rows.
    Use rowsof by itself if you can, as it's faster.

    rows: Anything that you can iterate over and get rows.
          A row is also iterable, and expected to be strings.
          Could be a rowsof generator.

    Keyword arguments:
    *comment* A character that, if it appears at the beginning of a row,
              indicates that the row should be skipped
    *header*  Number of rows to skip at the beginning of the file.
    *blank*   If True, ignore blank rows.  They are an error otherwise.
    """

    comment = kwargs.get("comment", None)
    header = kwargs.get("header", 0)
    blank = kwargs.get("blank", False)

    for row in rows:
        if header > 0:
            header -= 1
            continue
        if blank and len(row) == 0:
            continue

        try:
            if comment is not None and row[0].startswith(comment):
                continue
        except AttributeError as err:
            if "startswith" in err.message:
                pass # couldn't do starswith, maybe row is floats?
            else:
                raise err

def cli(args):
    """ command-line interface, execute the comparison """

    if not any([a is None for a in [args.blank, args.header, args.comment]]):
        tables = [filter_input( rowsof(fn, args.delimiter), 
                                blank=args.blank,
                                header=args.header,
                                comment=args.comment)
                      for fn in args.input]
    tables = [rowsof(fn, args.delimiter) for fn in args.input]

    if args.epsilons is not None and args.objectives is not None:
        if len(args.epsilons) != len(args.objectives):
            msg = "{0} epsilons specified for {1} objectives".format(
                    len(args.epsilons), len(args.objectives))
            raise SortParameterError(msg)
    epsilons = args.epsilons

    try:
        archive = eps_sort(tables, args.objectives, epsilons)
    except SortInputError as sie:
        table = args.input[sie.table]
        msg = sie.message.replace("input", table)
        raise SortInputError(msg, sie.row, table)

    with open(args.output, 'w') as fp:
        if args.print_only_objectives and args.objectives is not None:
            for row in archive.archive:
                obj = [row[ii] for ii in args.objectives]
                fp.write(args.delimiter.join(obj))
                fp.write("\n")
        else:
            for row in archive.archive:
                fp.write(args.delimiter.join(row))
                fp.write("\n")

if __name__ == "__main__":
    cli(get_args(sys.argv))
