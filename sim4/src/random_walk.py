from __future__ import print_function, division
from numpy import *
from matplotlib.pyplot import *
import scipy.stats
import sys, os
import pickle 

from inspect import currentframe, getframeinfo

# Additional bibs
import argparse
import time

# Own bibs
from random_numbers import *

# Commandline options
parser = argparse.ArgumentParser()

parser.add_argument( "--Xlcg", type=double, help="initial value for the linear congruental generator\nif not set: time.time() is used" )

args = parser.parse_args()

# initialize LCG
if args.Xlcg:
	init_LCG(args.Xlcg)
else:
	init_LCG(do2int(time.time()))

# 1D-random walk
N = 1000
x = random_walk(N=N)

# Plot
plot(range(0,N+1),x,label=r'random walk')
xlabel('time')
ylabel('position')
#show()

# Box-Mueller
print(BM())
show_BM_hist()
