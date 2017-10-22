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
    return x, v

# Define the main simulation loop, takes the timestep as an argument
def sim_loop(dt, flag):
    traj = []
    t = 0
    x = x_init.copy()
    v = v_init.copy()
    a = compute_forces(x)/m
    traj.append(x.copy())
    while t <= 1.0:
        if flag==-1:
            x, v = step_euler(x, v, dt)
        if flag==0:
            x, v = step_eulersym(x, v, dt)
        if flag==1:
            x, v = step_vv(x, v, a, dt)
        t += dt
        traj.append(x.copy())
    return np.array(traj)

# Run the simulation
step = 0.01
for i in range(1):
    nptraj = sim_loop(step, 1)
    plt.plot(nptraj[:,0,2]-nptraj[:,0,1],nptraj[:,1,2]-nptraj[:,1,1], label = str(i))
plt.title("Trajectory of the Moon for Different Time Steps: Simple Euler")
plt.legend()
plt.show()
