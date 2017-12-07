from __future__ import print_function
from numpy import *
from matplotlib.pyplot import *
import sys, os
import pickle
import argparse

# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--M", type=int, help="window size | default: M=10")
parser.add_argument("--datafile", type=str, help="datafilename | default: '../dat/ljsim.dat'")
parser.add_argument("--teq", type=double, help="equilibration time [s] | calculates averages after equilibration time")

args = parser.parse_args()

# get M from command line
if args.M:
	MM = args.M
else:
	MM = 10
print(MM)

# get filename on command line
if args.datafile:
	datafilename = args.datafile
else:
	datafilename = '../dat/ljsim.dat'

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
			Oa[k] += O[l]
	return Oa/(M+1.)
	
def compute_mean_value(O,keq):
	N = len(O)
	Om = 0					# mean value
	for k in xrange(keq,N):	# loop over teq -> tmax
		Om += O[k]
	return Om/(N-keq)

# Plot Runing Averages
print("Calculating avrage values...")
ts = array(ts)
Es = array(Es)
Ts = array(Ts)
Ps = array(Ps)
Ea = zeros(shape(Es))
Ea[:,0] = compute_running_average(Es[:,0],MM)
Ea[:,1] = compute_running_average(Es[:,1],MM)
Ea[:,2] = compute_running_average(Es[:,2],MM)
Ta = compute_running_average(Ts,MM)
Pa = compute_running_average(Ps,MM)
print("Finished.")
print("Plotting...")
# Energies
plot(ts, Ea[:,0],'g-',label=r'E_{pot}')
plot(ts, Ea[:,1],'r-',label=r'E_{kin}')
plot(ts, Ea[:,2],'b-',label=r'E_{tot}')
xlabel("Time t [s]")
ylabel("Energy E")
legend()
savefig('../dat/avEnergies_M{}.png'.format(MM))
close()
# Total Energy
plot(ts, Ea[:,2],'b-',label=r'E_{tot}')
xlabel("Time t [s]")
ylabel("Energy E")
legend()
savefig('../dat/avTotal_Energy_M{}.png'.format(MM))
close()
# Temperature
plot(ts, Ta,'g-',label=r'Temperature T')
xlabel("Time t [s]")
ylabel("Temperature T")
#legend()
savefig('../dat/avTemperature_M{}.png'.format(MM))
close()
# Pressure
plot(ts, Pa,'g-',label=r'Pressure P')
xlabel("Time t [s]")
ylabel("Pressure P")
#legend()
savefig('../dat/avPressure_M{}.png'.format(MM))
close()
print("Finished.")

# Compute mean values
if args.teq:
	print("Computing mean values...")
	# compute index of teq
	keq = 0
	while ts[keq] < args.teq:
		keq += 1
	# Compute mean values
	E_potm = compute_mean_value(Es[:,0],keq)
	E_kinm = compute_mean_value(Es[:,1],keq)
	E_totm = compute_mean_value(Es[:,2],keq)
	Tm = compute_mean_value(Ts,keq)
	Pm = compute_mean_value(Ps,keq)
	print("Mean Values for t >= {}:\n\tE_pot={}\n\tE_kin={}\n\tE_tot={}\n\tT={}\n\tP={}".format(ts[keq],E_potm, E_kinm, E_totm, Tm, Pm))
	print("Finished.")
