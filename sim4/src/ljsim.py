from __future__ import print_function, division
from numpy import *
from matplotlib.pyplot import *
import scipy.stats
import sys, os
import pickle 

from inspect import currentframe, getframeinfo

# Additional bibs
import argparse

# Own bibs
from random_numbers import *

# Argparse command line options
parser = argparse.ArgumentParser()

parser.add_argument( '--ID', type=str, help='simulation ID' )

args = parser.parse_args()

# SYSTEM CONSTANTS
# timestep
dt = 0.01
# length of run
trun = 1000.0
# desired temperature
T = 0.3
# total number of particles
N = 1
# side length of the system
L = 10.

# SIMULATION CONSTANTS
# warmup?
freq_andersen = 0.1
tau_berendsen = 3*dt
gamma_langevin = 0.3

# number of steps to do before the next measurement
measurement_stride = 100
# cutoff length

# OPEN SIMULATION FILES
if not args.ID:
    print('ERROR: simuation ID missing.\nUse --ID <ID>')
    sys.exit(2)
simulation_id = args.ID

log_vel= False
vtffilename = '../dat/{}.vtf'.format(simulation_id)
velfilename = '../dat/{}.vel'.format(simulation_id)
datafilename = '../dat/{}.dat'.format(simulation_id)

def compute_temperature(v):
    _, N = v.shape
    Tm = (v*v).sum()/(3*N)
    return Tm


def step_vv(x, v, f, dt, xup):
    global rcut, skin

    # update positions
    x += v*dt + 0.5*f * dt*dt

    # half update of the velocity
    v += 0.5*f * dt
        
    # for this excercise no forces from other particles
    f = zeros_like(x)

    # second half update of the velocity
    v += 0.5*f * dt

    return x, v, f, xup


# SET UP SYSTEM OR LOAD IT
# check whether data file already exists
if os.path.exists(datafilename):
    print("Reading data from {}.".format(datafilename))
    datafile = open(datafilename, 'rb')
    step, t, x, v, ts, Tms = pickle.load(datafile)
    datafile.close()
    print("Restarting simulation at t={}...".format(t))
else:
    print("Starting simulation...")
    t = 0.0
    step = 0
    x = np.array([[0.1], [2], [0]])
    v = np.array([[0.1], [0], [1]])

 # variables to cumulate data
    ts = []
    Tms = []

# check whether vtf file already exists
if os.path.exists(vtffilename):
    print("Opening {} to append new timesteps...".format(vtffilename))
    vtffile = open(vtffilename, 'a')

else:
    print("Creating {}...".format(vtffilename))
    # create a new file and write the structure
    vtffile = open(vtffilename, 'a')

    # write the structure of the system into the file: 
    # N particles ("atoms") with a radius of 0.5
    vtffile.write('atom 0:{} radius 0.5\n'.format(N-1))
    vtffile.write('pbc {} {} {}\n'.format(L, L, L))
    
    # write out that a new timestep starts
    vtffile.write('timestep\n')
    # write out the coordinates of the particles
    for i in range(N):
        vtffile.write("{} {} {}\n".format(x[0,i], x[1,i], x[2,i]))

# check whether velocity file already exists
if os.path.exists(velfilename):
    print("Opening {} to append new timesteps...".format(velfilename))
    velfile = open(velfilename, 'a')

else:
    print("Creating {}...".format(velfilename))
    # create a new file
    velfile = open(velfilename, 'a')

# main loop
xup = x.copy()
f = zeros_like(x)

tmax = t+trun



print("Simulating until tmax={}...".format(tmax))

while t < tmax:
    if  simulation_id.startswith('berendsen'):
        # change this
        x, v, f, xup = step_vv(x, v, f, dt, xup)

    elif simulation_id.startswith('langevin'):
        # change this
        x, v, f, xup = step_vv(x, v, f, dt, xup)

    elif simulation_id.startswith('andersen'):
        # change this
        x, v, f, xup = step_vv(x, v, f, dt, xup)
    else:
        raise Exception("Please implement the integrator for " + simulation_id)

    t += dt
    step += 1

    if step % measurement_stride == 0:
        Tm = compute_temperature(v)
        print("t={}, T_m={}".format(t, Tm))

        ts.append(t)
        Tms.append(Tm)

        for i in range(N):
            velfile.write("{}\t{}\t{}\n".format(v[0,i], v[1,i], v[2,i]))
        
        # write out that a new timestep starts
        vtffile.write('timestep\n')
        # write out the coordinates of the particles
        for i in range(N):
            vtffile.write("{} {} {}\n".format(x[0,i], x[1,i], x[2,i]))
            


# at the end of the simulation, write out the final state
print("Writing simulation data to {}.".format(datafilename))
datafile = open(datafilename, 'wb')
pickle.dump([step, t, x, v, ts, Tms ], datafile)
datafile.close()

# close vtf file
print("Closing {}.".format(vtffilename))
vtffile.close()
velfile.close()

print("Finished simulation.")
