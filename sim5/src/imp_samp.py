import numpy as np
import matplotlib.pyplot as plt
import argparse
from mon_car import *
from metropolis import *

# Initial state
parser = argparse.ArgumentParser()
parser.add_argument( '--x0', type=float, default=10*np.random.random() - 5, help='Initial state' )
args = parser.parse_args()

# Importance sampling
N = 100000                      # number of trial moves
dxs = [0.1, 1., 10., 100.]      # list of dx
states = []                     # list of the states
accept = []                     # acceptance rates
x0 = args.x0                    # initial state (random)
print('initial state:   x0 = {}'.format(x0))

for dx in dxs:
    trial_move_set_dx( dx )     
    state, acc = metropolis( N, runge, trial_move_runge, x0 )
    states.append( state )
    accept.append( acc )

# plot metropolis numbers
bins = np.linspace(-5,5,100)
x = np.linspace(-5,5,100)
y = runge(x) / runge_int(5.,-5.)
splt = [221, 222, 223, 224]

plt.figure(1)
for k in range(4):
    plt.subplot(splt[k])
    plt.plot( x, y, label=r'runge (normed)', linewidth=2, color='blue' )
    plt.hist( np.array( states[k] ), bins, normed=1, label=r'$\Delta x={}$\t, acc={}%'.format(dxs[k], int(accept[k]*1000)/10), color='red' )
    plt.xlim([-5,5])
    plt.title(r'$\Delta x={}$, $\alpha={}$%'.format(dxs[k], int(accept[k]*1000)/10) )
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.5), fancybox=True)
    plt.ylabel('frequency of $x$')
    plt.xlabel('state $x$')
    
# plot runge function
plt.tight_layout()
plt.savefig('../fig/metroplot.png',bbox_inches='tight')
