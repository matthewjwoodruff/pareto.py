from __future__ import division
import sys
import numpy as np
import math
import argparse
from sys import exit

class Archive(object):
    def __init__(self, epsilons, oindices):
        self.archive = []
        self.boxes = [] # remember boxes
        self.objectives = [] # remember objectives
        self.epsilons = epsilons
        self.oindices = oindices
        self.nobj = len(oindices)
        self.itobj = range(self.nobj)

    def add(self, solution, sobj, sbox):
        self.archive.append(solution)
        self.objectives.append(sobj)
        self.boxes.append(sbox)

    def remove(self, index):
        self.archive.pop(index)
        self.objectives.pop(index)
        self.boxes.pop(index)

    def sortinto(self, solution):
        """
        Sort a solution into the archive.  Add it if it's nondominated
        w.r.t current solutions.
        """
        # this code has a lot of early exits, keep control flow in mind
        # break -- stop iterating and do not evaluate any more of the loop
        # continue -- start the next loop iteration immediately
        # return -- immediately stop executing the function

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
                asize -=1
                continue # while

        # if you get here, then no archive solution has dominated this one
        self.add(solution, sobj, sbox)

# Get command line arguments and check for errors
parser = argparse.ArgumentParser(description='Nondomination Sort for Multiple Files')
parser.add_argument('--objectives', type=int, nargs='+', required=False, help='Objective Columns (zero-indexed)') 
parser.add_argument('--epsilons', type=float, nargs='+', required=False, help='Epsilons, one per objective')
parser.add_argument('-o', '--output', type=str, required=True, help='Output Filename')
parser.add_argument('-i', '--input', type=str, required=True, nargs='+', help='Input filenames')
parser.add_argument('--delimiter', type=str, required=False, default=' ', help='Input column delimiter')
parser.add_argument('--print-only-objectives', action='store_true', default=False, required=False, help='Print only objectives in output')
parser.add_argument('--precision', type=int, required=False, default=8, help='Output floating-point precision')
args = parser.parse_args()

# Get the first input file to check size
tempinput = np.loadtxt(args.input[0], delimiter=args.delimiter)

# If objectives not given, use all columns
if not args.objectives:
  args.objectives = range(tempinput[0,:].size)
# If given objectives are out of bounds, exit with error
elif (min(args.objectives) < 0 or max(args.objectives) > tempinput[0,:].size):
  msg = "Error: One or more objective values exceed input matrix bounds\n"
  sys.stdout.write(msg)
  exit()
if args.epsilons:
  if len(args.epsilons) != len(args.objectives): # make sure number of epsilons matches number of objectives
    msg =  "Error: Number of epsilon values must match number of objectives.\n"
    sys.stdout.write(msg)
    exit()
else:
  args.epsilons = [1e-9]*len(args.objectives) # if epsilons not defined, use 1e-9 for all objectives

archive = Archive(args.epsilons, args.objectives)

# for each file in argument list ...
for filename in args.input:
  solutions = np.loadtxt(filename, delimiter=args.delimiter)

  # for each line in file (new candidate solution) ...
  for i in xrange(0, solutions[:,0].size):
    candidateSolution = solutions[i,:]
    archive.sortinto(candidateSolution)

# Convert ParetoSet to nparray and print it
ParetoSet = np.array(archive.archive)
if args.print_only_objectives:
  ParetoSet = ParetoSet[:, args.objectives]
np.savetxt(args.output, ParetoSet, delimiter=args.delimiter, fmt='%.' + str(args.precision) + 'e')
