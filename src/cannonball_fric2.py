# Simulate a cannonball
import numpy as np
import matplotlib.pyplot as plt

# Constants
m = 2.0
g = 9.81
dt = 0.1

# Define functions
def compute_forces(x, v, y, vw):
    f_fric = -y*(v - np.array([vw, 0.0]))
    f = np.array([0.0, -m*g])+f_fric
    return f

def step_euler(x, v, y, vw):
    f = compute_forces(x, v, y, vw)
    x += v*dt
    v += f*dt/m    
    return x, v

def sim_loop(y, vw):
    traj = []
    t = 0
    x = np.array([0.0, 0.0])
    v = np.array([50.0, 50.0])
    traj.append(x.copy())
    
    while x[1] >= 0.0:
        x, v = step_euler(x, v, y, vw)
        t += dt
        traj.append(x.copy())
    return traj

vw = 0.0
while vw >= -200.0:
    traj = sim_loop(0.1, vw)
    traj = np.array(traj)
    plt.plot(traj[:,0], traj[:,1], label = r'$v_\mathrm{w}$ = ' + str(vw))
    vw -= 40.0



plt.axis([0,600,0,140])
plt.title('Trajectory of a Cannonball with Various Wind Speeds')
plt.xlabel('x position')
plt.ylabel('y position')
plt.legend()
plt.savefig('../fig/cannonball3.png', dpi=300, bbox_inches='tight')
