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
import numpy as np
import math
import argparse

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
        self.oindices = oindices
        self.nobj = len(oindices)
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

    def sortinto(self, solution):
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

        sobj = [solution[ii] for ii in self.oindices]
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
    parser.add_argument('--precision', type=int, required=False, default=8,
                        help='Output floating-point precision')
    return parser.parse_args(argv)

class SortInputError(Exception):
    """ Information about a defective input """
    def __init__(self, msg, row, table):
        super(SortError, self).__init__(msg)
        self.row = row
        self.table = table

def eps_sort(tables, objectives, epsilons):
    archive = Archive(epsilons, objectives)

    # for each file in argument list ...
    for counter in range(len(tables)):
        solutions = tables[counter]
        # for each line in file (new candidate solution) ...
        for i in range(0, solutions[:,0].size):
            try:
                candidateSolution = solutions[i,:]
                archive.sortinto(candidateSolution)
            except IndexError:
                msg = "Not enough columns in row {0} of input {1}".format(
                                                                    i, counter)
                raise SortInputError(msg, i, counter)
            except ValueError as ve:
                msg = "{0} on row {1} of input {2}".format(
                                                        ve.message, i, counter)
                raise SortInputError(msg, i, counter)

    return archive

def cli(args):
    """ command-line interface, execute the comparison """
    tables = [np.loadtxt(fn, delimiter = args.delimiter) for fn in args.input]
    if args.objectives is None:
        objectives = range(tables[0][0,:].size)
    else:
        objectives = args.objectives

    if args.epsilons:
        if len(args.epsilons) != len(objectives):
            msg =  "Error: Number of epsilon values must match "\
                   "number of objectives.\n"
            sys.stdout.write(msg)
            exit()
        epsilons = args.epsilons
    else:
        # if epsilons not defined, use 1e-9 for all objectives
        epsilons = [1e-9]*len(objectives) 

    try:
        archive = eps_sort(tables, objectives, epsilons)
    except SortInputError as sie:
        table = args.input[sie.table]
        raise SortInputError(sie.message, sie.row, table)

    # Convert ParetoSet to nparray and print it
    ParetoSet = np.array(archive.archive)
    if args.print_only_objectives:
        ParetoSet = ParetoSet[:, objectives]
    np.savetxt(args.output, ParetoSet, delimiter=args.delimiter, fmt='%.' + str(args.precision) + 'e')

if __name__ == "__main__":
    cli(get_args(sys.argv))
