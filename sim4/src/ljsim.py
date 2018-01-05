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
parser.add_argument( '--gamma', type=float, default=0.3, help='Friction coefficient\nDefault: 0.3' )
parser.add_argument( '--nu', type=float, default=0.1, help='Frequency of stochastic collisions\nDefault: 0.1' )
parser.add_argument( '--tau', type=float, default=3.0, help='Constant for berendsen in terms of dt\nDefault: 3.0s' )
parser.add_argument( '--tstat', type=float, default=1.0, help='Desired temperature for thermostat\nDefault: 1.0' )
parser.add_argument( '--restart', action='store_false', default=True, help='Do you want to restart the simulation instead of continuing?' )

args = parser.parse_args()

### INITIALIZATION ###############################################################################

# SET GLOBALS
global gamma;   gamma   = args.gamma
global nu;      nu      = args.nu
global tau;     tau     = args.tau
global T_des;   T_des   = args.tstat

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

### FUNCTIONS ####################################################################################

def compute_temperature(v):
    _, N = v.shape
    Tm = (v*v).sum()/(3*N)
    return Tm

def rand_vel(T):
    '''
    Creates 3D vector of normal distributed random numbers with sigma=sqrt(T)
    Because the expectation value of T for the normal distribution is sigma**2
    '''
    return sqrt(T)*random.randn(3)


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
    
def step_vv_langevin(x, v, f, dt, xup):
    '''
    velocity verlet for langevin thermostat
    '''
    global rcut, skin, gamma, T_des

    # update positions
    x += v*dt*(1-0.5*gamma*dt) + 0.5*f*dt*dt

    # half update of the velocity
    v += -0.5*gamma*dt + 0.5*f*dt
       
    # for this excercise no forces from other particles
    f = sqrt( 24*T_des*gamma/dt ) * ( random.random(x.shape)-0.5 )

    # second half update of the velocity
    v += 0.5*f*dt
    v /= 1+0.5*gamma*dt

    return x, v, f, xup
    
def step_vv_andersen(x, v, f, dt, xup):
    '''
    velocity verlet for anderson thermostat
    '''
    global rcut, skin, nu, T_des

    # update positions
    x += v*dt + 0.5*f * dt*dt

    # half update of the velocity
    v += 0.5*f * dt
        
    # for this excercise no forces from other particles
    f = zeros_like(x)

    # second half update of the velocity
    v += 0.5*f * dt
    
    # random velocity replacing:
    for k in range(0,v.shape[1]): # (only one particle to test, but extendable or more)
        if random.random() < nu*dt: 
            v[:,k] = rand_vel(T_des)

    return x, v, f, xup
   
def step_vv_berendsen(x, v, f, dt, xup):
    '''
    velocity verlet for berendsen thermostat
    '''
    global rcut, skin, tau, T_des
    
    # update positions
    x += v*dt + 0.5*f * dt*dt

    # half update of the velocity
    v += 0.5*f * dt
        
    # for this excercise no forces from other particles
    f = zeros_like(x)

    # second half update of the velocity
    v += 0.5*f * dt
    
    # velocity rescaling
    T_act = compute_temperature(v)
    v *= sqrt( 1 + ( T_des/T_act - 1 )/tau )

    return x, v, f, xup  
   
### PREPARING FOR SIMULATION #####################################################################
    
# SET UP SYSTEM OR LOAD IT
# check whether data file already exists
if os.path.exists(datafilename) and args.restart:
    print("Reading data from {}.".format(datafilename))
    datafile = open(datafilename, 'rb')
    step, t, x, v, ts, Tms, vs = pickle.load(datafile)
    datafile.close()
    print("Restarting simulation at t={}...".format(t))
else:
    print("Starting simulation...")
    t = 0.0
    step = 0
    x = np.array([[0.1], [2.], [0.]])
    v = np.array([[0.1], [0.], [1.]])

    # variables to cumulate data
    ts = []
    Tms = []
    vs = []

# check whether vtf file already exists
if os.path.exists(vtffilename) and args.restart:
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
if os.path.exists(velfilename) and args.restart:
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

### SIMULATION ###################################################################################

print("Simulating until tmax={}...".format(tmax))

while t < tmax:
    if  simulation_id.startswith('berendsen'):
        # Using berendsen thermostat
        x, v, f, xup = step_vv_berendsen(x, v, f, dt, xup)

    elif simulation_id.startswith('langevin'):
        # Using langevin thermostat 
        x, v, f, xup = step_vv_langevin(x, v, f, dt, xup)

    elif simulation_id.startswith('andersen'):
        # Using andersen thermostat
        x, v, f, xup = step_vv_andersen(x, v, f, dt, xup)
    else:
        raise Exception("Please implement the integrator for " + simulation_id)

    t += dt
    step += 1

    if step % measurement_stride == 0:
        Tm = compute_temperature(v)
        #print("t={}, T_m={}".format(t, Tm))

        ts.append(t)
        Tms.append(Tm)
        vs.append(linalg.norm(v,axis=0))

        for i in range(N):
            velfile.write("{}\t{}\t{}\n".format(v[0,i], v[1,i], v[2,i]))
        
        # write out that a new timestep starts
        vtffile.write('timestep\n')
        # write out the coordinates of the particles
        for i in range(N):
            vtffile.write("{} {} {}\n".format(x[0,i], x[1,i], x[2,i]))
            
### MEASUREMENT ##################################################################################

# at the end of the simulation, write out the final state
print("Writing simulation data to {}.".format(datafilename))
datafile = open(datafilename, 'wb')
pickle.dump([step, t, x, v, ts, Tms, vs ], datafile)
datafile.close()

# close vtf file
print("Closing {}.".format(vtffilename))
vtffile.close()
velfile.close()

print("Finished simulation.")

# PLOTS
print('Plotting.')

# Temperature over time
plot(ts,Tms)
xlabel("Time t")
ylabel("Temperature Tm")
savefig('../dat/{}_Tm.png'.format(simulation_id))
close()

# Velocity distribution
r = linspace(0, 10, 1000)
gauss_fun = gauss_3d(r, T_des**2)       # function from random_numbers.py
plot(r, gauss_fun, label=r'gauss       : $T_d$$_e$$_s={}$'.format(T_des))
hist(array(vs).flatten(), linspace(0,10,100), normed=1, label=r'simulation: $T_d$$_e$$_s={}$'.format(T_des))
xlabel('velocity v')
ylabel('frequenzy of v in time')
legend()
savefig('../dat/{}_vv.png'.format(simulation_id))
close()

print('Finished')
