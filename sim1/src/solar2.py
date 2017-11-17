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

# Define the main simulation loop, takes the timestep as an argument
def sim_loop(dt):
    traj = []
    t = 0
    x = x_init.copy()
    v = v_init.copy()
    traj.append(x.copy())
    while t <= 1.0:
        x, v = step_euler(x, v, dt)
        t += dt
        traj.append(x.copy())
    return np.array(traj)

# Run the simulation
for i in range(4):
    step = 0.0001+i*0.0003
    nptraj = sim_loop(step)
    plt.plot(nptraj[:,0,2]-nptraj[:,0,1],nptraj[:,1,2]-nptraj[:,1,1], label = str(step))
plt.title("Trajectory of the Moon for Different Time Steps: Simple Euler")
plt.legend()
plt.xlabel("x position [AU]")
plt.ylabel("y position [AU]")

plt.savefig('../fig/solar2.png', dpi=300, bbox_inches='tight')
