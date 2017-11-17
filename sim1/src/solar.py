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

def step_euler(x, v, dt):
    f = compute_forces(x)
    x += v*dt
    v += f*dt/m    
    return x, v

def sim_loop(dt):
    traj = []
    t = 0
    x = x_init
    v = v_init
    traj.append(x.copy())
    
    while t <= 1.0:
        x, v = step_euler(x, v, dt)
        t += dt
        traj.append(x.copy())
    return np.array(traj)

nptraj = sim_loop(0.0001)

for i in range(M):
    plt.plot(nptraj[:,0,i],nptraj[:,1,i], label = name[i])
plt.title("Trajectory of Solar System for 1 Year, dt = 0.0001")
plt.legend(loc="lower right", ncol=2)
plt.xlabel("x position [AU]")
plt.ylabel("y position [AU]")
plt.savefig('../fig/solar1.png', dpi=300, bbox_inches='tight')
