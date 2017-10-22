# Simulate the solar system in 2d
import numpy as np
import matplotlib.pyplot as plt
import pickle, gzip
import numpy.linalg as la

# Read in initial planetary positions, velocities, and masses
datafile = gzip.open('Solar_system.pkl.gz')
name, x_init, v_init, m, g = pickle.load(datafile)
datafile.close()

# Compute array shape
N,M = x_init.shape

# Define functions
def compute_forces(x):
    f = np.zeros(x_init.shape)
    for i in range(M):
        for j in range(M):
            if i != j:
                r_ij = x[:,i] - x[:,j]
                f[:,i] += -g*m[i]*m[j]*r_ij/(la.norm(r_ij)**3)
    return f

# Define a simple euler integration scheme
def step_euler(x, v, dt):
    f = compute_forces(x)
    x += v*dt
    v += f*dt/m    
    return x, v

# Define a symplectic Euler algorithm 
def step_eulersym(x, v, dt):
    f = compute_forces(x)
    v += f*dt/m    
    x += v*dt
    return x, v

# Define a velocity Verlet integrator
def step_vv(x, v, a, dt):
    x += v*dt+a*dt**2/2
    v += a*dt/2
    a = compute_forces(x)/m
    v += a*dt/2
    return x, v, a

# Define the main simulation loop, takes the timestep as an argument
def sim_loop(dt, tmax, flag):
    traj = []
    t = 0
    x = x_init.copy()
    v = v_init.copy()
    a = compute_forces(x)/m
    traj.append(x.copy())
    while t <= tmax:
        if flag==0:
            x, v = step_euler(x, v, dt)
        if flag==1:
            x, v = step_eulersym(x, v, dt)
        if flag==2:
            x, v, a = step_vv(x, v, a, dt)
        t += dt
        traj.append(x.copy())
    return np.array(traj)

# Set time constants
step = 0.01
tmax = 1.0

# Simulate and plot simple Euler integrator
nptraj = sim_loop(step, tmax, 0)
moon_traj = nptraj[:,:,2]-nptraj[:,:,1]
plt.subplot(131)
plt.plot(moon_traj[:,0],moon_traj[:,1], label = "Simple Euler")

# Simulate and plot symplectic Euler integrator
nptraj = sim_loop(step, tmax, 1)
moon_traj = nptraj[:,:,2]-nptraj[:,:,1]
plt.subplot(132)
plt.plot(moon_traj[:,0],moon_traj[:,1], label = "Symplectic Euler")

# Simulate and plot velocity Verlet integrator
nptraj = sim_loop(step, tmax, 2)
moon_traj = nptraj[:,:,2]-nptraj[:,:,1]
plt.subplot(133)
plt.plot(moon_traj[:,0],moon_traj[:,1], label = "Velocity Verlet")


plt.legend()
plt.show()
