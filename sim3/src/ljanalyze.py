from __future__ import print_function
from numpy import *
from matplotlib.pyplot import *
import sys, os
import pickle
import argparse

# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--M", type=int, help="window size")
parser.add_argument("--datafile", type=string, help="datafilename")

args = parser.parse_args()

if args.M:
	MM = args.M
else:
	MM = 10

# get filename on command line
datafilename = args.datafile

# check whether data file exists
if not os.path.exists(datafilename):
    print("ERROR: {} doesn't exist.".format(datafilename))
    sys.exit(1)

print("Reading data from {}.".format(datafilename))
datafile = open(datafilename, 'r')
ts, Es, Ts, Ps, x, v = pickle.load(datafile)
datafile.close()

def compute_running_average(O,M):
	m = M/2
	N = len(O)
	Oa = zeros(N)						# new empty array for average
	Oa[:m] = nan						# set values to NaN 
	Oa[-m:] = nan
	for k in xrange(m,N-m):				# loop over all observables
		for l in xrange(k-m,k+m+1):		# loop over window values
			Oa[k] += O(l)
	return Oa/(M+1.)

# Plot Runing Averages
Ea[:,0] = compute_running_average(Es[:,0],MM)
Ea[:,1] = compute_running_average(Es[:,1],MM)
Ea[:,2] = compute_running_average(Es[:,2],MM)
Ta = compute_running_average(Ts,MM)
Pa = compute_running_average(Ps,MM)
# Energies
plot(ts, Ea[:,0],'g-',label=r'E_{pot}')
plot(ts, Ea[:,1],'r-',label=r'E_{kin}')
plot(ts, Ea[:,2],'b-',label=r'E_{tot}')
xlabel("Time t [s]")
ylabel("Energy E")
legend()
savefig('../dat/avEnergies.png')
close()
# Total Energy
plot(ts, Ea[:,2],'b-',label=r'E_{tot}')
xlabel("Time t [s]")
ylabel("Energy E")
legend()
savefig('../dat/avTotal_Energy.png')
close()
# Temperature
plot(ts, Ta,'g-',label=r'Temperature T')
xlabel("Time t [s]")
ylabel("Temperature T")
#legend()
savefig('../dat/avTemperature.png')
close()
# Pressure
plot(ts, Pa,'g-',label=r'Pressure P')
xlabel("Time t [s]")
ylabel("Pressure P")
#legend()
savefig('../dat/avPressure.png')
close()
