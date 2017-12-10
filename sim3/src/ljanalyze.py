from __future__ import print_function
from numpy import *
from matplotlib.pyplot import *
import sys, os
import pickle
import argparse

# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--M", type=int, default=10, help="window size | default: M=10")
parser.add_argument("--datafile", type=str, default='../dat/ljsim.dat', help="datafilename | default: '../dat/ljsim.dat'")
parser.add_argument("--teq", type=double, help="equilibration time | calculates averages after equilibration time")
parser.add_argument("--tlim", type=double, nargs='+', help="type in the time limits you want to plot")

args = parser.parse_args()

MM = args.M
datafilename = args.datafile
bins = np.linspace(0.8, 5.0, 100)
# check whether data file exists
if not os.path.exists(datafilename):
    print("ERROR: {} doesn't exist.".format(datafilename))
    sys.exit(1)

# Getting additional information from datafilename    
ending = datafilename[12:-4]

print("Reading data from {}.".format(datafilename))
datafile = open(datafilename, 'r')
ts, Es, Ts, Ps, x, v, hs = pickle.load(datafile)
datafile.close()

def compute_running_average(O,M):
	m = M/2
	N = len(O)
	Oa = zeros(len(O))				# new empty array for average
	Oa[:m] = nan					# set values to NaN 
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

def compute_mean_rdf(O,keq):
    N = O.shape[0]
    Om = empty(O.shape[1])
    for k in xrange(keq,N):
        Om += O[k,:]
    return Om/(N-keq)

# Plot Runing Averages
print("Calculating average values...")
ts = array(ts)
Es = array(Es)
Ts = array(Ts)
Ps = array(Ps)
hs = array(hs)
Ea = zeros(shape(Es))
Ea[:,0] = compute_running_average(Es[:,0],MM)
Ea[:,1] = compute_running_average(Es[:,1],MM)
Ea[:,2] = compute_running_average(Es[:,2],MM)
Ta = compute_running_average(Ts,MM)
Pa = compute_running_average(Ps,MM)
print("Finished.")
print("Plotting...")
# Energies
ylabel("Energy E")
if args.tlim:
	t1 = 0
	while ts[t1] < args.tlim[0]:
		t1 += 1
	t2 = 0
	while ts[t2] < args.tlim[1]:
		t2 += 1
	xlim(args.tlim)
	plot(ts[t1:t2], Ea[t1:t2,0],'g-',label=r'E_{pot}')
	plot(ts[t1:t2], Ea[t1:t2,1],'r-',label=r'E_{kin}')
	plot(ts[t1:t2], Ea[t1:t2,2],'b-',label=r'E_{tot}')
else:
	plot(ts, Ea[:,0],'g-',label=r'E_{pot}')
	plot(ts, Ea[:,1],'r-',label=r'E_{kin}')
	plot(ts, Ea[:,2],'b-',label=r'E_{tot}')
xlabel("Time t")
legend()
savefig('../dat/avEnergies{}_M{}.png'.format(ending,MM))
close()
# Total Energy
if args.tlim:
	xlim(args.tlim)
	plot(ts[t1:t2], Ea[t1:t2,2],'b-',label=r'E_{tot}')
else:
	plot(ts, Ea[:,2],'b-',label=r'E_{tot}')
xlabel("Time t")
ylabel("Energy E")
legend()
savefig('../dat/avTotal_Energy{}_M{}.png'.format(ending,MM))
close()
# Temperature
if args.tlim:
	xlim(args.tlim)
	plot(ts[t1:t2], Ta[t1:t2],'g-',label=r'Temperature T')
else:
	plot(ts, Ta,'g-',label=r'Temperature T')
xlabel("Time t")
ylabel("Temperature T")
#legend()
savefig('../dat/avTemperature{}_M{}.png'.format(ending,MM))
close()
# Pressure
if args.tlim:
	xlim(args.tlim)
	plot(ts[t1:t2], Pa[t1:t2],'g-',label=r'Pressure P')
else:
	plot(ts, Pa,'g-',label=r'Pressure P')
xlabel("Time t")
ylabel("Pressure P")
#legend()
savefig('../dat/avPressure{}_M{}.png'.format(ending,MM))
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
        hm = compute_mean_rdf(hs,keq)
        plot(bins,hm, 'g-', label=r'Radial Distribution Function')
        xlabel("Radial Distance r")
        ylabel("Probability")
        savefig('../dat/meanRDF.png')
        close()
        
	print("Mean Values for t >= {}:\n\tE_pot={}\n\tE_kin={}\n\tE_tot={}\n\tT={}\n\tP={}".format(ts[keq],E_potm, E_kinm, E_totm, Tm, Pm))
	print("Finished.")
