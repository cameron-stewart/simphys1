from __future__ import print_function
from numpy import *
from matplotlib.pyplot import *
import sys, os
import pickle
from lj import *
import argparse
import scipy.misc as binomial

def compute_histogram(r):
        h, bins = histogram(r, bins=100, range=(0.8, 5.),density=True)
        return  h*vol/(4*pi*bins[:-1]*bins[:-1])

# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--cont", type=double, help="continue calculation with for cont further time")
parser.add_argument("--time", type=double, help="How long do you want to run the simulation? | default time=10")
parser.add_argument("--tstat", type=double, help="Uses a thermostat with a given temperature")
parser.add_argument("--warm", type=double, help="Use the warming up | pass force")
parser.add_argument("--ctstat", type=double, help="continue for the simulation with tstat")
parser.add_argument("--cwarm", type=double, help="continue for the simulation with warm")

args = parser.parse_args()

# SYSTEM CONSTANTS
# density
density = 0.316
# timestep
dt = 0.01
# number of particles per side for cubic setup
n = 10

# SIMULATION CONSTANTS
# skin size
skin = 0.4
# number of steps to do before the next measurement
measurement_stride = 100
# cutoff length
rcut = 2.5
# potential shift
shift = -0.016316891136

# Filename addings
if args.ctstat:
	tstat = "_T" + str(args.ctstat)
elif args.tstat:
	tstat = "_T" + str(args.tstat) 
else:
	tstat = ""
if args.cwarm:
	warm = "_F" + str(args.cwarm)
elif args.warm:
	warm = "_F" + str(args.warm)
else:
	warm = ""
		
# VTF filename 
vtffilename = "../dat/ljsim{}{}.vtf".format(tstat,warm)
# DATA filename 
datafilename = "../dat/ljsim{}{}.dat".format(tstat,warm)        

# COMPUTED CONSTANTS
# total number of particles
N = n*n*n
# volume of the system
volume = N/density
# side length of the system
L = volume**(1./3.)
# volume for RDF
vol = 4*pi*(5)**3/3

r = empty((binomial.comb(N,2),1))

# Import previous data
if args.cont:
	# open datafile
	datafile = open(datafilename,'r')
	ts, Es, Ts, Ps, x, v, hs = pickle.load(datafile)
	datafile.close()
	
	# length of run
	t = ts[-1]
	tmax = t+args.cont
	step = 0

else:
	# length of run
	if args.time:
		tmax = args.time
	else:
		tmax = 10.0

	print("Starting simulation...")
	# particle positions on cubic lattice
	t = 0.0
	step = 0
	
	# Initialize particle position
	if args.warm: 	
		# The warming up random positions
		x = L*random.random((3,N))
	else:
		x = empty((3,N))
		l = L/n
		count = 0
		for i in range(n):
			for j in range(n):
				for k in range(n):
				    x[:,count] = [i*l, j*l, k*l]
				    count += 1
		        
	# random particle velocities
	v = 0.1*(2.0*random.random((3,N))-1.0)

	# variables to cumulate data
	ts = []
	Es = []
	Ts = []
	Ps = []
        hs = []

print("density={}, L={}, N={}".format(density, L, N))

# check whether vtf file already exists
print("Creating {}...".format(vtffilename))
# create a new file and write the structure
if args.cont:
	vtffile = open(vtffilename, 'a')
else:
	vtffile = open(vtffilename, 'w')
	# write the structure of the system into the file: 
	# N particles ("atoms") with a radius of 0.5
	vtffile.write('atom 0:{} radius 0.5\n'.format(N-1))
	vtffile.write('pbc {} {} {}\n'.format(L, L, L))


# write out that a new timestep starts
vtffile.write('timestep\n')
# write out the coordinates of the particles
for i in range(N):
    vtffile.write("{} {} {}\n".format(x[0,i], x[1,i], x[2,i]))

def step_vv(x, v, f, dt, xup):
    global rcut, skin

    # update positions
    x += v*dt + 0.5*f * dt*dt

    # compute maximal position update
    # vectorial
    dx = x - xup
    # square
    dx *= dx
    # sum up 
    dx = dx.sum(axis=0)
    # test whether the neighbor list needs to be rebuilt
    if max(dx) > (0.5*skin)**2:
        rebuild_neighbor_lists(x, rcut+skin)
        xup = x.copy()
    
    # half update of the velocity
    v += 0.5*f * dt
        
    # compute new forces
    # we assume that m=1 for all particles
    f = compute_forces(x)
    if args.warm:
        force_capping(f,args.warm)

    # second half update of the velocity
    v += 0.5*f * dt

    return x, v, f, xup

# main loop
set_globals(L, N, rcut, shift)
rebuild_neighbor_lists(x, rcut+skin)
xup = x.copy()
f = compute_forces(x)
if args.warm:
    force_capping(f,args.warm)

print("Simulating until tmax={}...".format(tmax))
while t < tmax:
    x, v, f, xup = step_vv(x, v, f, dt, xup)
    t += dt
    step += 1

    if step % measurement_stride == 0:
        E_pot, E_kin, E_tot = compute_energy(x, v)
        T = 2*E_kin/(3*N)
        P = compute_pressure(E_kin, x)
        print("t={}:\n\tE_pot={}\n\tE_kin={}\n\tE_tot={}\n\tT={}\n\tP={}".format(t, E_pot, E_kin, E_tot, T, P))
        compute_distances(x, r)
        h = compute_histogram(r)

        ts.append(t)
        Es.append([E_pot,E_kin,E_tot])
        Ts.append(T)
        Ps.append(P)
        hs.append(h)

        # write out that a new timestep starts
        vtffile.write('timestep\n')
        # write out the coordinates of the particles
        for i in range(N):
            vtffile.write("{} {} {}\n".format(x[0,i], x[1,i], x[2,i]))
            
        # velocity rescaling thermostat
        if args.tstat:
        	velocity_rescaling(args.tstat,T,v)
        	
        # rescale fcap
        if args.warm:
        	args.warm *= 1.1
                
# close vtf file
print("Closing {}.".format(vtffilename))
vtffile.close()

# write out simulation data
print("Writing simulation data to {}.".format(datafilename))
datafile = open(datafilename, 'w')
pickle.dump([ts, Es, Ts, Ps, x, v], datafile)
datafile.close()

print("Finished simulation.")

print("Plotting...")
ts = array(ts)
Es = array(Es)
Ts = array(Ts)
Ps = array(Ps)
# Energies
plot(ts, Es[:,0],'g-',label=r'E_{pot}')
plot(ts, Es[:,1],'r-',label=r'E_{kin}')
plot(ts, Es[:,2],'b-',label=r'E_{tot}')
xlabel("Time t")
ylabel("Energy E")
legend()
savefig('../dat/Energies{}{}.png'.format(tstat,warm))
close()
# Total Energy
plot(ts, Es[:,2],'b-',label=r'E_{tot}')
xlabel("Time t")
ylabel("Energy E")
legend()
savefig('../dat/Total_Energy{}{}.png'.format(tstat,warm))
close()
# Temperature
plot(ts, Ts,'g-',label=r'Temperature T')
xlabel("Time t")
ylabel("Temperature T")
#legend()
savefig('../dat/Temperature{}{}.png'.format(tstat,warm))
close()
# Pressure
plot(ts, Ps,'g-',label=r'Pressure P')
xlabel("Time t")
ylabel("Pressure P")
#legend()
savefig('../dat/Pressure{}{}.png'.format(tstat,warm))
close()

print("Finished.")
