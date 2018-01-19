import numpy as np
from mon_car import *
import matplotlib.pyplot as plt

N = 100000                  # number of trial moves
dxs = [0.1, 1., 10., 100.]  # list of dx
dim = 1                     # dimension of x
states = []                 # list of the states
accept = []                 # acceptance rates
for k in range(4):
    x0 = 2*dxs[k]*np.random.random(1) - dxs[k]
    trial_move_set_dx( dxs[k] )     
    state, acc = metropolis( N, runge, trial_move_runge, x0 )
    states.append( state )
    accept.append( acc )

# plot metropolis numbers
bins = np.linspace(-5,5,100)
for k in range(4):
    plt.hist( np.array( states[k] ), bins, normed=1, label=r'$\Delta x={}$\t, acc={}%'.format(dxs[k], int(accept[k]*1000)/10) )
    
# plot runge function
x = np.linspace(-5,5,100)
y = runge(x) / runge_int(5.,-5.)
plt.plot( x, y, label=r'runge-function')
plt.xlim([-5,7])
plt.legend()
plt.show()
