from __future__ import print_function
from numpy import *
from matplotlib.pyplot import *
import sys, os
import pickle

# get filename on command line
if len(sys.argv) != 2:
    print("Usage: python {} FILE".format(sys.argv[0]))
    sys.exit(2)
datafilename = sys.argv[1]

# check whether data file exists
if not os.path.exists(datafilename):
    print("ERROR: {} doesn't exist.".format(datafilename))
    sys.exit(1)

print("Reading data from {}.".format(datafilename))
datafile = open(datafilename, 'r')
ts, Es = pickle.load(datafile)
datafile.close()

ts = array(ts)
Es = array(Es)

plot(ts, Es)
show()
