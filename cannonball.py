# Simulate a cannonball
import numpy as np
import matplotlib.pyplot as plt

# Constants
m = 2.0
g = 9.81
dt = 0.1

# Define functions
def compute_forces(x):
    f = np.array([0.0, -m*g])
    return f

def step_euler(x, v, dt):
    f = compute_forces(x)
    x += v*dt
    v += f*dt/m    
    return x, v

# Initialize variables
t = 0.0
x = np.array([0.0, 0.0])
v = np.array([50.0, 50.0])
traj = []

while x[1] >= 0.0:
    x, v = step_euler(x, v, dt)
    t += dt
    traj.append(x.copy())

# Plot
traj = np.array(traj)
plt.plot(traj[:,0], traj[:,1], '-')
plt.show()

#test comment
