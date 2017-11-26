from numpy import *
from matplotlib.pyplot import *
from ljlib import *

def compute_forces(x):
    """Compute and return the forces acting onto the particles,
    depending on the positions x."""
    global epsilon, sigma
    _, N = x.shape
    f = zeros_like(x)
    for i in range(1,N):
        for j in range(i):
            # distance vector
            rij = x[:,j] - x[:,i]
            rij = image(rij)
            fij = compute_cutoff_force(rij)
            f[:,i] -= fij
            f[:,j] += fij
    return f

def compute_energy(x, v):
    """Compute and return the total energy of the system with the
    particles at positions x."""
    _, N = x.shape
    E_pot = 0.0
    E_kin = 0.0
    # sum up potential energies
    for i in range(1,N):
        for j in range(i):
            # distance vector
            rij = x[:,j] - x[:,i]
            rij = image(rij)
            E_pot += compute_cutoff_potential(rij)
    # sum up kinetic energy
    for i in range(N):
        E_kin += 0.5 * dot(v[:,i],v[:,i])
    return E_pot + E_kin
    
def step_vv(x, v, f, dt):
    # update positions
    x += v*dt + 0.5*f * dt*dt
    # half update of the velocity
    v += 0.5*f * dt
        
    # compute new forces
    f = compute_forces(x)
    # we assume that m=1 for all particles

    # second half update of the velocity
    v += 0.5*f * dt

    return x, v, f

def image(rij):
    """Return the minimum image of particle j for particle i"""
    rij -= L*rint(rij/L)
    return rij

# constants
dt = 0.01
tmax = 20.0
L = 10.0

# running variables
t = 0.0

# particle positions
x = zeros((3,2))
x[:,0] = [3.9, 3, 4]
x[:,1] = [6.1, 5, 6]


# particle velocities
v = zeros((3,2))
v[:,0] = [-2.0, -2.0,-2.0]
v[:,1] = [2.0, 2.0, 2.0]

f = compute_forces(x)

# variables to cumulate data
traj = []
Es = []

# number of particles
N = x.shape[1]

# open the trajectory file
vtffile = open('../dat/ljbillards.vtf', 'w')
# write the structure of the system into the file:
# N particles ("atoms") with a radius of 0.5
vtffile.write('atom 0:{} radius 0.5\n'.format(N-1))
vtffile.write('pbc 10.0 10.0 10.0\n')

# main loop
while t < tmax:
    x, v, f = step_vv(x, v, f, dt)
    t += dt

    traj.append(x.copy()-L*floor(x.copy()/L))
    
    Es.append(compute_energy(x, v))
    
    # write out that a new timestep starts
    vtffile.write('timestep\n')
    # write out the coordinates of the particles
    for i in range(N):
        vtffile.write("{} {} {}\n".format(x[0,i], x[1,i], x[2,i]))

vtffile.close()

traj = array(traj)

# plot the trajectory
figure()
for i in range(N):
    plot(traj[:,0,i], traj[:,1,i], 'o', label='{}'.format(i))
axes().set_aspect('equal')
legend()

# plot the total energy
figure()
xlabel("Time step")
ylabel("Total energy")
plot(Es)

show()
