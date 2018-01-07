from __future__ import print_function, division
from numpy import *
from matplotlib.pyplot import *
import scipy.stats
import sys, os
import pickle 

from inspect import currentframe, getframeinfo

# Additional bibs
import argparse

# Own files
from random_numbers import *

# Argparse command line options
parser = argparse.ArgumentParser()

parser.add_argument( '--ID', type=str, help='Simulation ID' )
parser.add_argument( '--dtmax', type=int, help='Maximum of step width (in terms of ts[1]-ts[0])' )
parser.add_argument( '--linreg', type=int, nargs='+', help='--linreg Min Max: area in which it is linear' )

args = parser.parse_args()

### READING DATA-FILE #############################################################################

if not args.ID:
    print('ERROR: simuation ID missing.\nUse --ID <ID>')
    sys.exit(2)
simulation_id = args.ID

datafilename = '../dat/{}.dat'.format(simulation_id)

print("Reading data from {}.".format(datafilename))
datafile = open(datafilename, 'rb')
step, t, x, v, ts, Tms, xs, vs = pickle.load(datafile)
datafile.close()

xs = array(xs)
N = shape(xs)[0]
d = shape(xs)[1]
if args.dtmax:
    dtmax = args.dtmax
else:
    dtmax = N-1
if args.linreg:
    linreg = args.linreg
else:
    linreg = [0,dtmax]

### CALCULATE MSD #################################################################################

def calc_msd( ndt, x=xs):
    '''
    calculates the msd for time steps of dt = ts[ndt] - ts[0]
    '''
    n = shape(x)[0]/ndt
    k = 1
    msd = 0.
    while k < n:
        msd += sum( ( x[ndt*k,:,0] - x[ndt*(k-1),:,0] )**2 )
        k += 1
    return msd / (k-1)  #TODO errorbars
    
msd_vec = zeros(dtmax)
dt = zeros(dtmax)
for k in range(1,dtmax+1):
    msd_vec[k-1] = calc_msd( k, xs )
    dt[k-1] = ts[k]-ts[0]

# Plot
plot(dt,msd_vec,'o')

### FIT ON DATA ###################################################################################

p = polyfit( dt[linreg[0]:linreg[1]], msd_vec[linreg[0]:linreg[1]], 1 )
D = p[0] * 0.5 / float(d)
print('D = {}'.format(D))

# Plot
x = linspace(linreg[0],linreg[1],1000)
plot( x, p[0]*x+p[1], '-r')
xlabel("Duration $\Delta t$")
ylabel("MSD of $x$")
savefig('../dat/{}_MSD.png'.format(simulation_id))
close()

### VACF ##########################################################################################

vs=array(vs)
#TODO expand for VACF

    


