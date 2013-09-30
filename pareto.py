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

# Alright, I think we can agree that everything can be made to work 
# really well *if we specify objectives*
# I'm at the point where I want to rip out _initialsortinto and 
# purge the Archive class of crap.
# Now, if we don't specify objectives, what breaks?
# The archive needs somebody to tell it what the objectives are.
# That's right and proper.  So eps_sort does the inference, and
# in doing so it takes the first row off of the first table, so
# it puts that row in its own table.
# OK, so now it's a problem for eps_sort, not the archive.
# And the problem is, we've already monkeyed with the row by the
# time it gets there.  Maybe I need to push the "secret extra
# table" back to the cli function?  Because this is only ever
# a problem for command-line invocation.  Otherwise if you've 
# monkeyed with your table before calling eps_sort, you're 
# naturally obliged to specify the objectives your damn self.
# So when does this happen, that I lift a row from the first 
# input table?  If I do it before I call filter_inputs, it messes
# up the line numbers for the whole first table, plus I might
# sample the table by mistake.  So I need to have the filtered
# row instead.  But this is the row I've monkeyed with already!
# Now how do I figure out how long it was before I pulled my 
# shenanigans?  Maybe filter_inputs need to yield the pure 
# row plus monkey business, and let cli combine them.
# It breaks generator chaining, but it makes everything else 
# easier.

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
        self.sortinto = self._initialsortinto
        self.epsilons = epsilons
        self.oindices = oindices
        self.nobj = 0
        self.itobj = range(self.nobj)

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
        number of objectives and what the epsilons are.
        """
        # should this method look for some kind of a guard that indicates 
        # the end of the "natural" row?
        # something that has a string representation of "empty string?"
        # grrr, you get an extra separator if you pull this sort of
        # nonsense.
        # Either that or you need to do an output filter to strip it
        # back out.
        # what if an optional manifest gets put at the end of the first row?
        # Really don't like the idea of tying those two together so closely.
        # The guard is better, it's basically free, there only needs to be
        # one of it, and it only needs to be stripped out if you put it in.
        # still, it's a hack and I don't like the way it interferes with the
        # purity of the archive
        if self.oindices is None:
            self.nobj = len(solution)
            self.oindices = range(self.nobj)
        else:
            self.nobj = len(self.oindices)

        self.itobj = range(self.nobj)

        if self.epsilons is None:
            self.epsilons = [1e-9]*self.nobj
        elif len(self.epsilons) != self.nobj:
            msg = "{0} epsilons specified, but {1} objectives".format(
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

class MixedArchive(Archive):
    """ an archive with maximization objectives """
    def __init__(self, epsilons, oindices, maximize, maximize_all):
        """
        epsilons: sizes of epsilon boxes to use in the sort
        oindices: indicate which indices in a solution are objectives
        maximize: column indices to maximize
        maximize_all: maximize all columns
        """
        super(self, MixedArchive).__init__(epsilons, oindices)
        self.maximize = maximize
        self.maximize_all = maximize_all

    def _initalsortinto(self, solution):
        """
        Establish number of objectives and which are maximized
        """
        if self.oindices is None:
            self.nobj = len(solution)
            self.oindices = range(self.nobj)
        else:
            self.nobj = len(self.oindices)
        if self.maximize_all:
            self.maximize = self.oindices
        if self.maximize is None:
            self.maximize = []
        super(self, MixedArchive)._initialsortinto(solution)

    def _realsortinto(self, solution):
        """
        cache the 

        


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
    else:
        return range(first, second+1)

class SortInputError(Exception):
    """ Information about a defective input """
    def __init__(self, msg, row, table):
        super(SortInputError, self).__init__(msg)
        self.row = row
        self.table = table

def eps_sort(tables, objectives, epsilons, **kwargs):
    """
    Perform an epsilon-nondominated sort
    tables: input data, must support row iteration
    objectives: list of column indices in which objectives can be found,
                if None default to all columns
    epsilons: list of epsilons for the sort, if None default to 1e-9

    Keyword arguments:
    *maximize*      invert the sense of these columns
    *maximize_all*  invert the sense of all columns
    """
    maximize = kwargs.get("maximize", None)
    maximize_all = kwargs.get("maximize_all", False)

    archive = Archive(epsilons, objectives)

    maximize = kwargs.get("maximize", None)False
    minimize_objectives = objectives

    # chain flip onto any other generators to augment rows with negated values
    if kwargs.get("maximize_all", False):
        tables = [flip(table, maximize_all=True) for table in tables]
        # all of the objectives are now after the end of the row
        # but we have a bootstrapping problem here -- we need to 
        # push a solution into a phony archive to figure out how 
        # big it is, so that we can create the real archive

        # so here's what we're going to do: snarf one solution out of the 
        # first table, push it into a phony archive, and make a new table
        # that's just that solution

        # it doesn't break contribution, since contribution has already
        # happened at this point

        # Wait, what if somebody wants to do contribution and maximize all
        # objectives?  Oy.  Need to be aggressive about catching valueerrors?
        # hmmm. Actually, what happens if you want to do contribution 
        # without specifying objectives, period?  bollocks.  That's a bug.

        # all this is to say that the lazy-objectives folks who specify 
        # nothing make this whole project about 10x harder!

    elif maximize is not None:
        tables = [flip(table, columns=maximize) for table in tables]
        # easier case, we know which columns get flipped, so we change 
        # the indices in the list of the objectives


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
                                                str(ve), rownumber, counter)
                raise SortInputError(msg, rownumber, counter)

    return archive

def rowsof(stream, delimiter):
    """
    Generator function yielding rows read from a stream. (Lazy input.)
    Avoids having to read the whole file at once.
    """
    try:
        while True:
            line = next(stream)
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
    *comment*       A character that, if it appears at the beginning of a row,
                    indicates that the row should be skipped
    *header*        Number of rows to skip at the beginning of the file.
    *blank*         If True, ignore blank rows.  They are an error otherwise.
    *contribution*  A tag to append to each row indicating where it came from.
                    Can be anything printable.
    *number*        Include line number in contribution if True.  Defaults to
                    False.  Counts from 1 (because usual use case is files,
                    where lines are conventionally numbered from 1.)
    """

    comment = kwargs.get("comment", None)
    header = kwargs.get("header", 0)
    if header is None:
        header = 0
    blank = kwargs.get("blank", False)
    contribution = kwargs.get("contribution", None)
    number = kwargs.get("number", False)

    counter = 1

    for row in rows:
        counter += 1
        if header > 0:
            header -= 1
            continue
        if blank:
            if len(row) == 0:
                continue
            elif len(row) == 1 and len(row[0]) == 0:
                continue

        try:
            if comment is not None:
                iscomment = False
                for commentchar in comment:
                    if row[0].startswith(commentchar):
                        iscomment = True
                if iscomment:
                    continue
        except AttributeError as err:
            if "startswith" in str(err):
                # couldn't do starswith, maybe row is floats?
                pass
            else:
                raise

        if contribution is not None:
            # uh-oh, what if they didn't specify objectives?
            row.append(str(contribution))
            if number:
                row.append(str(counter))

        yield row

def use_filter(args):
    """ return True if we need to use filtered input """
    if args.header is not None:
        return True
    if args.comment is not None:
        return True
    if args.blank:
        return True
    if args.contribution:
        return True
    return False

def flip(rows, **kwargs):
    """
    Append inverted values to each row in the input

    Keyword arguments:
    *columns*       which columns to invert, invert all if not specified
    *maximize_all*  invert all columns if True
    """
    counter = -1
    row = []
    try:
        columns = kwargs.get("columns", None)
        if columns is not None:
            import copy
            for row in rows:
                counter += 1
                maxrow = copy.copy(row)
                for cc in columns:
                    maxrow.append(0.0 - float(maxrow[cc]))
                yield maxrow
        elif kwargs.get("maximize_all", True):
            import copy
            for row in rows:
                counter += 1
                maxrow = copy.copy(row)
                maxrow.extend([0.0 - float(val) for val in row])
                yield maxrow
        else: # why are you using this function again?
            for row in rows:
                yield row
    except ValueError as ve:
        msg = "On row {0} of input: {1}, encountered error {2}".format(
                    counter, row, ve)
        raise SortInputError(msg)

def get_args(argv):
    """ Get command line arguments """
    prog = argv.pop(0)
    parser = argparse.ArgumentParser(prog=prog,
        description='Nondomination Sort for Multiple Files')
    parser.add_argument('inputs', type=argparse.FileType('r'), nargs='+',
                        help='input filenames, use - for standard input')
    parser.add_argument('-o', '--objectives', type=intrange, nargs='+',
                        help='objective columns (zero-indexed)')
    parser.add_argument('-m', '--maximize', type=intrange, nargs='+',
                        help='objective columns to maximize (zero-indexed)')
    parser.add_argument('-M', '--maximize-all', action="store_true",
                        help='maximize all objectives')
    parser.add_argument('-e', '--epsilons', type=float, nargs='+',
                        help='epsilons, one per objective')
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
    parser.add_argument("-c", "--comment", type=str, nargs="+",
                        help="skip lines starting with this character")
    parser.add_argument("--header", type=int,
                        help="number of header lines to skip")
    parser.add_argument("--contribution", action="store_true",
                        help="append filename where solution originated")
    parser.add_argument("--line-number", action="store_true",
                        help="also append line number to solution if "\
                             "--contribution is used.")
    args = parser.parse_args(argv)

    if args.objectives is not None:
        objectives = []
        for indexrange in args.objectives:
            objectives.extend(indexrange)
        args.objectives = objectives

    if args.maximize is not None:
        cols = []
        for indexrange in args.maximize:
            cols.extend(indexrange)
        args.maximize = cols

    if args.maximize_all: # override args.maximize if maximizing all
        args.maximize = args.objectives
        args.maximize_all = False

    if args.tabs:
        args.delimiter = "\t"

    return args

def cli(args):
    """ command-line interface, execute the comparison """

    tables = [rowsof(fp, args.delimiter) for fp in args.inputs]

    if use_filter(args):
        if args.contribution:
            tags = [i.name for i in args.inputs]
        else:
            tags = [None] * len(args.inputs)
        tables = [filter_input(table, blank=args.blank, header=args.header,
                               comment=args.comment, contribution=tag,
                               number=args.line_number)
                  for table, tag in zip(tables, tags)]

    if args.epsilons is not None and args.objectives is not None:
        if len(args.epsilons) != len(args.objectives):
            msg = "{0} epsilons specified for {1} objectives".format(
                    len(args.epsilons), len(args.objectives))
            raise SortParameterError(msg)
    epsilons = args.epsilons

    try:
        archive = eps_sort(tables, args.objectives, epsilons, 
                           maximize = args.maximize, 
                           maximize_all = args.maximize_all)
    except SortInputError as sie:
        table = args.inputs[sie.table].name
        msg = str(sie).replace("input", table)
        raise SortInputError(msg, sie.row, table)

    if args.print_only_objectives and args.objectives is not None:
        for row in archive.archive:
            obj = [row[ii] for ii in args.objectives]
            args.output.write(args.delimiter.join(obj))
            args.output.write("\n")
    else:
        for row in archive.archive:
            args.output.write(args.delimiter.join(row))
            args.output.write("\n")

    args.output.close()

if __name__ == "__main__":
    cli(get_args(sys.argv))
