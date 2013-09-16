from __future__ import division
import sys
import numpy as np
import math
import argparse
from sys import exit

# Epsilon box comparison function
def compare(solution1, solution2):

  dominate1 = False
  dominate2 = False

  for i in xrange(0, len(args.objectives)):
    epsilon = args.epsilons[i]
    index1 = math.floor(solution1[args.objectives[i]]/epsilon)
    index2 = math.floor(solution2[args.objectives[i]]/epsilon)
#    if secretcounter == magicnumber:
#      sys.stdout.write("({0}, {1}) ".format(index1, index2))

    if index1 < index2:
      dominate1 = True
#      if secretcounter == magicnumber:
#        sys.stdout.write("l ")
      if dominate2:
#        if secretcounter == magicnumber:
#          sys.stdout.write("* ")
        return 0
    elif index1 > index2:
#      if secretcounter == magicnumber:
#        sys.stdout.write("g ")
      dominate2 = True
      if dominate1:
#        if secretcounter == magicnumber:
#          sys.stdout.write("* ")
        return 0

  # If one clearly dominates the other, return
  if dominate1:
#    if secretcounter == magicnumber:
#      sys.stdout.write("L ")
    return -1
  elif dominate2:
#    if secretcounter == magicnumber:
#      sys.stdout.write("G ")
    return 1

#  if secretcounter == magicnumber:
#    sys.stdout.write("= ")
  # If neither dominates the other in any objective, they are in the same epsilon box
  if not dominate1 and not dominate2:
    dist1 = 0.0
    dist2 = 0.0

    for i in xrange(0, len(args.objectives)):
      epsilon = args.epsilons[i]
      index1 = math.floor(solution1[args.objectives[i]]/epsilon)
      index2 = math.floor(solution2[args.objectives[i]]/epsilon)
      dist1 += math.pow(solution1[args.objectives[i]] - index1*epsilon, 2.0);
      dist2 += math.pow(solution2[args.objectives[i]] - index2*epsilon, 2.0);

    if (dist1 < dist2): # compare squared distances
#      if secretcounter == magicnumber:
#        sys.stdout.write("{0} - {1}".format(dist1, dist2))
#        sys.stdout.write("L ")
      return -1
    else:
#      if secretcounter == magicnumber:
#        sys.stdout.write("{0} - {1}".format(dist1, dist2))
#        sys.stdout.write("G ")
      return 1

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
  print "Error: One or more objective values exceed input matrix bounds"
  exit()
if args.epsilons:
  if len(args.epsilons) != len(args.objectives): # make sure number of epsilons matches number of objectives
    print "Error: Number of epsilon values must match number of objectives."
    exit()
else:
  args.epsilons = [1e-9]*len(args.objectives) # if epsilons not defined, use 1e-9 for all objectives

# Define ParetoSet (a list of nparrays) and build it from input files
ParetoSet = []; 

secretcounter = -1
magicnumber = 70

# for each file in argument list ...
for filename in args.input:
  solutions = np.loadtxt(filename, delimiter=args.delimiter)

  # for each line in file (new candidate solution) ...
  idominated = 0
  for i in xrange(0, solutions[:,0].size):
    secretcounter += 1
#    sys.stdout.write("\narchive size {0}, I dominated {1}\n".format(
#            len(ParetoSet), idominated))
    idominated = 0
    candidateSolution = solutions[i,:]
#    sys.stdout.write("{0}\n".format(" ".join([str(x) for x in candidateSolution])))
    addCandidateSolution = True
    n = len(ParetoSet)

    if n == 0:
#      sys.stdout.write("add {0}\n".format(" ".join([str(x) for x in candidateSolution])))
      ParetoSet.append(candidateSolution)
    else:
      p = 0
      while p < n:
        existingSolution = ParetoSet[p]
        flag = compare(existingSolution, candidateSolution)
        if flag > 0:
          del ParetoSet[p]
          n = n - 1
          idominated += 1
        elif flag < 0:
          addCandidateSolution = False
          break
        else:
          p = p + 1

      if addCandidateSolution:
#        sys.stdout.write("add {0}\n".format(" ".join([str(x) for x in candidateSolution])))
        ParetoSet.append(candidateSolution)

# Convert ParetoSet to nparray and print it
ParetoSet = np.array(ParetoSet)
if args.print_only_objectives:
  ParetoSet = ParetoSet[:, args.objectives]
np.savetxt(args.output, ParetoSet, delimiter=args.delimiter, fmt='%.' + str(args.precision) + 'e')
