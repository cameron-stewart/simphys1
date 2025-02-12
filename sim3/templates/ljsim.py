from __future__ import print_function
from numpy import *
from matplotlib.pyplot import *
import sys, os
import pickle 
from lj import *

# SYSTEM CONSTANTS
# density
density = 0.316
# timestep
dt = 0.01
# length of run
tmax = 10.0
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
# VTF filename 
vtffilename = "ljsim.vtf"
# DATA filename 
datafilename = "ljsim.dat"

# COMPUTED CONSTANTS
# total number of particles
N = n*n*n
# volume of the system
volume = N/density
# side length of the system
L = volume**(1./3.)

print("Starting simulation...")
# particle positions on cubic lattice
t = 0.0
step = 0

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

print("density={}, L={}, N={}".format(density, L, N))

# check whether vtf file already exists
print("Creating {}...".format(vtffilename))
# create a new file and write the structure
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
    f = compute_forces(x)
    # we assume that m=1 for all particles

    # second half update of the velocity
    v += 0.5*f * dt

    return x, v, f, xup

# main loop
set_globals(L, N, rcut, shift)
rebuild_neighbor_lists(x, rcut+skin)
xup = x.copy()
f = compute_forces(x)
print("Simulating until tmax={}...".format(tmax))
while t < tmax:
    x, v, f, xup = step_vv(x, v, f, dt, xup)
    t += dt
    step += 1

    if step % measurement_stride == 0:
        E_pot, E_kin, E_tot = compute_energy(x, v)
        print("t={}, E={}".format(t, E_tot))

        ts.append(t)
        Es.append(E_tot)

        # write out that a new timestep starts
        vtffile.write('timestep\n')
        # write out the coordinates of the particles
        for i in range(N):
            vtffile.write("{} {} {}\n".format(x[0,i], x[1,i], x[2,i]))

# close vtf file
print("Closing {}.".format(vtffilename))
vtffile.close()

# write out simulation data
print("Writing simulation data to {}.".format(datafilename))
datafile = open(datafilename, 'w')
pickle.dump([ts, Es], datafile)
datafile.close()

print("Finished simulation.")

print("Plotting...")
ts = array(ts)
Es = array(Es)
plot(ts, Es)
show()
print("Finished.")
