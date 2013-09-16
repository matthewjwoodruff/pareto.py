""" eps boxes of input file """

import sys
from numpy import floor

fp = open(sys.argv[1], 'r')

eps = [0.15, 30.0, 6.0, 0.03, 30.0, 3000.0, 150.0, 0.3, 3.0, 0.3]
eps = [0.1,0.1,0.1]

for line in fp:
    row = line.strip().split(" ")
    converted = [float(x) for x in row]
    boxes = [floor(x/e) for x,e in zip(converted, eps)]
    newrow = [str(box) for box in boxes]
    sys.stdout.write("{0}\n".format(" ".join(newrow)))
