"""
Copyright (C) 2013 Matthew Woodruff and Jon Herman.

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

@misc{woodruff_herman_2013_pareto,
    author = {Woodruff, Matthew and Herman, Jon},
    year = {2013},
    title = {pareto.py: a $\\varepsilon-nondomination$ sorting routine},
    howpublished = {https://github.com/matthewjwoodruff/pareto.py}
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

def get_args(argv):
    """ Get command line arguments """
    prog = argv.pop(0)
    parser = argparse.ArgumentParser(prog=prog,
        description='Nondomination Sort for Multiple Files')
    parser.add_argument('inputs', type=argparse.FileType('r'), nargs='+',
                        help='input filenames, use - for standard input')
    parser.add_argument('-o', '--objectives', type=intrange, nargs='+',
                        help='objective columns (zero-indexed)')
    parser.add_argument('-e', '--epsilons', type=float, nargs='+',
                        help='epsilons, one per objective')
    parser.add_argument('-m', '--maximize', type=intrange, nargs='+',
                        help='objective columns to maximize')
    parser.add_argument('-M', '--maximize-all', action="store_true",
                        help='maximize all objectives')
    parser.add_argument('--output', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='output filename, default to standard output')

    delimiters = parser.add_mutually_exclusive_group()
    delimiters.add_argument('-d', '--delimiter', type=str, default=' ',
                        help='input column delimiter, default to space (" ")')
    delimiters.add_argument('--tabs', action="store_true",
                        help="use tabs as delimiter")

    parser.add_argument('--print-only-objectives', action='store_true',
                        default=False, help='print only objectives in output')
    parser.add_argument("--blank", action="store_true",
                        help="skip blank lines")
    parser.add_argument("-c", "--comment", type=str, nargs="+", default=[],
                        help="skip lines starting with this character")
    parser.add_argument("--header", type=int, default=0,
                        help="number of header lines to skip")
    parser.add_argument("--contribution", action="store_true",
                        help="append filename where solution originated")
    parser.add_argument("--line-number", action="store_true",
                        help="also append line number to solution if "\
                             "--contribution is used.")
    parser.add_argument("--reverse-column-indices", action='store_true',
                        default=False, help='Reverse the order of column '\
                        'indices.  May be useful if your objectives are '\
                        'at the end of a row of unknown length.  Make sure '\
                        '-e and -m are consistent with the order you '\
                        'specify.')


    args = parser.parse_args(argv)

    args.objectives = rerange(args.objectives)
    args.maximize = rerange(args.maximize)

    if args.reverse_column_indices:
        args.objectives = [-1 - ob for ob in args.objectives]
        args.maximize = [-1 -ob for ob in args.maximize]

    if args.tabs:
        args.delimiter = "\t"

    return args

def rerange(intranges):
    """ convert a set of intranges into a list of integers """
    if intranges is None:
        return None
    thelist = []
    for therange in intranges:
        thelist.extend(therange)
    return thelist

def intrange(arg):
    """ convert a command-line argument to a list of integers """
    acceptable_chars = [str(x) for x in range(10)]
    acceptable_chars.append("-")

    partial = []
    first = None

    msg = "Could not convert {0} to index range.".format(arg)
    err = TypeError(msg)

    for char in arg:
        if char not in acceptable_chars:
            raise err
        if char == "-":
            if len(partial) == 0:
                raise err
            elif first is None:
                first = int("".join(partial))
                partial = []
            else: # this means there's a second -, which is not ok
                raise err
        else:
            partial.append(char)

    second = None
    if first is None:
        first = int("".join(partial))
    elif len(partial) == 0:
        raise err
    else:
        second = int("".join(partial))

    if second is None:
        return [first]
    elif second - first >= 0:
        return range(first, second+1)
    else:
        return range(first, second-1, -1)

class SortParameterError(Exception): pass

class Archive(object):
    """
    An archive of epsilon-nondominated solutions.
    Allows auxiliary information to tag along for the sort
    process.

    The eps_sort function provides a much more convenient interface than
    the Archive class.
    """
    def __init__(self, epsilons):
        """
        epsilons: sizes of epsilon boxes to use in the sort.  Number
                  of objectives is inferred by the number of epsilons.
        """
        self.archive = []       # objectives
        self.tagalongs = []     # tag-along data
        self.boxes = []         # remember for efficiency
        self.epsilons = epsilons
        self.itobj = range(len(epsilons)) # infer number of objectives

    def add(self, objectives, tagalong, ebox):
        """ add a solution to the archive, plus auxiliary information """
        self.archive.append(objectives)
        self.tagalongs.append(tagalong)
        self.boxes.append(ebox)

    def remove(self, index):
        """ remove a solution from the archive """
        self.archive.pop(index)
        self.tagalongs.pop(index)
        self.boxes.pop(index)

    def sortinto(self, objectives, tagalong=None):
        """
        Sort a solution into the archive.  Add it if it's nondominated
        w.r.t current solutions.

        objectives: objectives by which to sort.  Minimization is assumed.
        tagalong:   data to preserve with the objectives.  Probably the actual
                    solution is here, the objectives having been extracted
                    and possibly transformed.  Tagalong data can be *anything*.
                    We don't inspect it, just keep a reference to it for as
                    long as the solution is in the archive, and then return
                    it in the end.
        """
        # Here's how the early loop exits in this code work:
        # break:    Stop iterating the box comparison for loop because we know
        #           the solutions are in relatively nondominated boxes.
        # continue: Start the next while loop iteration immediately (i.e.
        #           jump ahead to the comparison with the next archive member).
        # return:   The candidate solution is dominated, stop comparing it to
        #           the archive, don't add it, immediately exit the method.

        ebox = [math.floor(objectives[ii] / self.epsilons[ii])
                for ii in self.itobj]

        asize = len(self.archive)

        ai = -1 # ai: archive index
        while ai < asize - 1:
            ai += 1
            adominate = False # archive dominates
            sdominate = False # solution dominates
            nondominate = False # neither dominates

            abox = self.boxes[ai]

            for oo in self.itobj:
                if abox[oo] < ebox[oo]:
                    adominate = True
                    if sdominate: # nondomination
                        nondominate = True
                        break # for
                elif abox[oo] > ebox[oo]:
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
            aobj = self.archive[ai]
            corner = [ebox[ii] * self.epsilons[ii] for ii in self.itobj]
            sdist = sum([(objectives[ii] - corner[ii]) **2
                         for ii in self.itobj])
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
        self.add(objectives, tagalong, ebox)

class SortInputError(Exception):
    """ Information about a defective input """
    def __init__(self, msg, row, table):
        super(SortInputError, self).__init__(msg)
        self.row = row
        self.table = table

def noannotation(table):
    """ produce solutions with no annotation from a table """
    empty = []
    for row in table:
        yield (row, empty)

def eps_sort(tables, objectives=None, epsilons=None, **kwargs):
    """
    Perform an epsilon-nondominated sort
    tables: input data, must support row iteration
    objectives: list of column indices in which objectives can be found,
                if None default to all columns
    epsilons: list of epsilons for the sort, if None default to 1e-9

    Keyword arguments:
    *maximize*      columns to maximize
    *maximize_all*  maximize all columns

    Duplicates some of cli() for a programmatic interface
    """
    tables = [noannotation(table) for table in tables]
    tables = [withobjectives(annotatedrows, objectives)
              for annotatedrows in tables]

    tomaximize = kwargs.get("maximize", None)
    maximize_all = kwargs.get("maximize_all", False)

    if tomaximize is not None or maximize_all:
        if objectives is None:
            mindices = tomaximize
        elif maximize_all:
            mindices = None
        else:
            mindices = [objectives.index(i) for i in tomaximize]
        tables = [maximize(solutions, mindices) for solutions in tables]

    tagalongs = eps_sort_solutions(tables, epsilons)

    return tagalongs

def eps_sort_solutions(tables, epsilons=None):
    """
    Perform an epsilon-nondominated sort
    tables: input (objectives, row) tuples
    epsilons: epsilon values for the objectives.  Assume 1e-9 if none
    """
    # slip the first row off the first table to figure out nobj
    objectives, row = next(tables[0])
    table = [(objectives, row)]
    tables = [table] + tables

    nobj = len(objectives)
    if epsilons is None:
        epsilons = [1e-9] * len(objectives)
    elif len(epsilons) != nobj:
        msg = "{0} epsilons, but {1} objectives".format(len(epsilons), nobj)
        raise SortParameterError(msg)

    archive = Archive(epsilons)

    for table in tables:
        for objectives, row in table:
            archive.sortinto(objectives, row)

    return archive.tagalongs

def attribution(stream, tag, number=False):
    """
    extract lines from stream and augment with tag
    """
    if number:
        linenumber = 0
        for line in stream:
            linenumber += 1
            line = line.strip()
            yield (line, [tag, str(linenumber)])
    else:
        for line in stream:
            line = line.strip()
            yield (line, [tag])

def noattribution(stream):
    """
    extract lines from stream and augment with null attribution
    """
    empty = []
    for line in stream:
        line = line.strip()
        yield (line, empty)

def filter_lines(annotatedlines, **kwargs):
    """
    remove commented, blank, and header lines
    """
    comment = kwargs.get("comment", [])
    header = kwargs.get("header", 0)
    blank = kwargs.get("blank", False)

    for line, annot in annotatedlines:

        # skip header lines
        if header > 0:
            header -= 1
            continue

        # skip comment lines
        iscomment = False
        for commentchar in comment:
            iscomment = iscomment or line.startswith(commentchar)
        if iscomment:
            continue

        # skip blank lines
        if blank and len(line) == 0:
            continue

        yield (line, annot)

def rowsof(annotatedlines, delimiter):
    """ split lines using delimiter, yielding annotated rows """
    for line, annot in annotatedlines:
        yield (line.split(delimiter), annot)

def withobjectives(annotatedrows, oindices):
    """ extract objectives and convert to float """
    if oindices is not None:
        for row, annot in annotatedrows:
            objectives = []
            for oo in oindices:
                objectives.append(float(row[oo]))
            row.extend(annot)
            yield objectives, row
    else:
        for row, annot in annotatedrows:
            objectives = [float(x) for x in row]
            row.extend(annot)
            yield objectives, row

def maximize(solutions, mindices=None):
    """
    mindices: which objectives to maximize.  If None, maximize all.

    These indices are indices into the list of objectives, not
    into the input row.  So if the objectives are 2, 3, 13, and 9,
    in that order, and you want to maximize column 2, specify 0
    to this function, and if you want to maximize column 13,
    specify 2 to this function.
    """
    if mindices is None:
        for objectives, row in solutions:
            objectives = [-x for x in objectives]
            yield objectives, row
    else:
        for objectives, row in solutions:
            for ii in mindices:
                objectives[ii] = 0 - objectives[ii]
            yield objectives, row

def cli(args):
    """ command-line interface, execute the comparison """
    if args.contribution:
        tables = [attribution(fp, fp.name, args.line_number)
                  for fp in args.inputs]
    else:
        tables = [noattribution(fp) for fp in args.inputs]

    if args.header > 0 or len(args.comment) > 0 or args.blank:
        tables = [filter_lines(annotatedlines, comment=args.comment,
                              header=args.header, blank=args.blank)
                  for annotatedlines in tables]

    tables = [rowsof(annotatedlines, args.delimiter)
              for annotatedlines in tables]

    tables = [withobjectives(annotatedrows, args.objectives)
              for annotatedrows in tables]

    if args.maximize is not None or args.maximize_all:
        if args.objectives is None:
            mindices = args.maximize
        elif args.maximize_all:
            mindices = None
        else:
            mindices = [args.objectives.index(i) for i in args.maximize]
        tables = [maximize(solutions, mindices) for solutions in tables]

    tagalongs = eps_sort_solutions(tables, args.epsilons)

    if args.print_only_objectives and args.objectives is not None:
        for row in tagalongs:
            obj = [row[ii] for ii in args.objectives]
            args.output.write(args.delimiter.join(obj))
            args.output.write("\n")
    else:
        for row in tagalongs:
            args.output.write(args.delimiter.join(row))
            args.output.write("\n")

    args.output.close()

if __name__ == "__main__":
    cli(get_args(sys.argv))
