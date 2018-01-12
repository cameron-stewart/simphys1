# Some things I learned in the tutorial

# random LCG as class
class Random(object):
    a = 1664525
    c = 1013904223
    m = 2**32
    def __init__(self, seed = 0):
        self.seed = seed
    def __call__(self):
        self.seed = (self.a*self.seed+self.c) % self.m
        return self.seed / m - 0.5
        

# The right velocity verlet (langevin)
# First position
x += v*dt*(1.-0.5*gamma*dt) + 0.5*f*dt*dt

# Then velocity:
v += -0.5*gamma*dt*v + 0.5*f*dt
v /= 1.+0.5*gamma*dt
# Correct velocity after 1/2

# Force (may depend on velocity (not in the task))    
f = sqrt( 24.*T_des*gamma/dt ) * ( random.random(x.shape) - 0.5 )

# Velocity 2/2
v += 0.5*f*dt / ( 1.+0.5*gamma*dt )


# andersen (vectorial "if request")
chosen = np.random.random(N) < freq*dt
# >> array([True, False, False, True, ...], dtype=boolean)
v[:, chosen] = np.random.normal( ... )
# Changes all values for which it is true.

# diffusion msd
np.linalg.norm( x[dt:dt*N:dt] - x[:dt*(N-1):dt] , axis=1)

# fit:
# option weight to pass errorbars. 

# autocorrelation:
# do not norm like I did
