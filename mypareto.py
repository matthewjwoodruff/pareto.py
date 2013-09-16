import sys
import pandas
import numpy
import argparse

def squared_distance(distances):
    """ squared euclidean distance metric """
    return sum([x**2 for x in distances])

def compare(x1, x2):
    """ replacement for cmp """
    if x1 == x2:
        return 0
    diff = x1 - x2
    return diff / abs(diff)

def nondom(r1, r2):
    """ nondominance comparison """
    return [compare(x1, x2) for x1, x2 in zip(r1, r2)]

class Comparator(object):
    def __init__(self, cols, eps):
        self.cols = cols
        self.eps = eps

    def boxes(self, row):
        """ box for each objective """
        return [numpy.floor(x/e) for x, e in zip(row, self.eps)]

    def remainder(self, box, row):
        """ distance to box corner for each objective """
        return [x - b*e for b, x, e in zip(box, row, self.eps)]

    def objectives(self, row):
        """ objectives of row """
        return [row[i] for i in self.cols]

    def compare(self, r1, r2):
        """ eps-nondominance comparison """
        global secretcounter
        global magicnumber

        o1 = self.objectives(r1)
        o2 = self.objectives(r2)

        b1 = self.boxes(o1)
        b2 = self.boxes(o2)

        dom = nondom(b1, b2)

        if all([x == 0 for x in dom]): # same box
            rem1 = self.remainder(b1, o1)
            rem2 = self.remainder(b1, o2)
            sd1 = squared_distance(rem1)
            sd2 = squared_distance(rem2)
#            if secretcounter == magicnumber:
#                sys.stdout.write("{0} - {1}".format(sd1, sd2))
            return compare(sd1, sd2)
        if not any([x > 0 for x in dom]):
            return -1 # r1 dominates r2
        if not any([x < 0 for x in dom]):
            return 1 # r1 dominated by r2
        return 0 # r1 and r2 do not dominate each other

def eps_sort(tables, objectives, epsilons):
    global secretcounter
    global magicnumber
    comp = Comparator(objectives, epsilons)
    archive = []
    secretcounter = -1
    magicnumber = 70

    for table in tables:
        idominated = 0
        for solution in table.iterrows():
            if secretcounter % 100 == 99:
                sys.stdout.write("\n")
            sys.stdout.write(".")
#            sys.stdout.write("\narchive size {0}, I dominated {1}\n".format(len(archive), idominated))
            idominated = 0
            secretcounter += 1
            solution = solution[1]
#            sys.stdout.write("{0}\n".format(" ".join([str(x) for x in solution])))
            obj = comp.objectives(solution)
            box = comp.boxes(obj)
#            sys.stdout.write("versus {0} / {1}\n".format(
#                ["{0:.3f}".format(x) for x in obj],
#                ["{0:.0f}".format(b) for b in box]
#            ))
            dominated_by_archive = False
            dominated_by_solution = []
            for ii in range(len(archive)):
                if ii % 10 == 0:
                    sys.stdout.write(",")
                comparison = comp.compare(solution, archive[ii])
#                if secretcounter == magicnumber and False:
#                    obj = comp.objectives(archive[ii])
#                    box = comp.boxes(obj)
#                    sys.stdout.write("       {0} / {1}\n".format(
#                        ["{0:.3f}".format(x) for x in obj],
#                        ["{0:.0f}".format(b) for b in box]
#                    ))
#                    obj = comp.objectives(solution)
#                    box = comp.boxes(obj)
#                    sys.stdout.write("versus {0} / {1}\n".format(
#                        ["{0:.3f}".format(x) for x in obj],
#                        ["{0:.0f}".format(b) for b in box]
#                    ))
#                    sys.stdout.write("{0} ".format(comparison))
                if comparison < 0: # solution dominates archive[ii]
                    dominated_by_solution.append(ii)
                    idominated += 1
                if comparison > 0: # archive dominates solution
                    dominated_by_archive = True
                    break
#            if secretcounter == magicnumber:
#                sys.stdout.write("\n")
            # purge dominated solutions from the archive
            archive = [archive[i] for i in range(len(archive))
                       if i not in dominated_by_solution]
            # add the solution if it was not dominated
            if not dominated_by_archive:
#                sys.stdout.write("add {0}\n".format(" ".join([str(x) for x in solution])))
                archive.append(solution)

#            sys.stdout.write("----\n")
            for xs in archive:
                obj = comp.objectives(xs)
                box = comp.boxes(obj)
#                sys.stdout.write( "{0} / {1}\n".format(
#                    ["{0:.3f}".format(x) for x in obj],
#                    ["{0:.0f}".format(b) for b in box]))
                    
    return archive

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=argparse.FileType("r"),
                        help="input files")
    parser.add_argument("-f", "--first-objective", type=int, default=0,
                        help="first objective column, default is zero") 
    parser.add_argument("-n", "--number-objectives", type=int, default=1,
                        help="how many objective columns (default is one)")
    parser.add_argument("-e", "--epsilons", nargs="+", type=float,
                        help="list of epsilons", default=[])
    return parser.parse_args()

def cli():
    args = get_args()
    fo = args.first_objective
    archive = eps_sort([pandas.read_table(t, sep=" ", header=None)
                        for t in args.inputs],
                        range(fo, fo + args.number_objectives), args.epsilons)
    for xs in archive:
        sys.stdout.write( "{0}\n".format(" ".join(["{0:.20g}".format(x) 
                                          for x in xs])))
    

if __name__ == "__main__":
    cli()

# vim:ts=4:sw=4:expandtab:ai:colorcolumn=70
