from __future__ import print_function, division
from numpy import *
from matplotlib.pyplot import *

# Additional bibs
import argparse
import time

# Own bibs
from random_numbers import *

# Commandline options
parser = argparse.ArgumentParser()

parser.add_argument( "--Xlcg", type=double, help="Initial value for the linear congruental generator\nif not set: time.time() is used" )

args = parser.parse_args()

# initialize LCG
if args.Xlcg:
	init_LCG(args.Xlcg)
else:
	init_LCG(do2int(time.time()))

# 1D-random walk
x = random_walk(N=1000)

# Box-Mueller
show_BM_hist(N=10000)

# Velocity distribution
show_vel_hist(N=1000)
