"""
Plots the LJ potential and force in order to test these functions

Author: Cameron Stewart
"""

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import ljlib as lj

d = np.linspace(0.85, 2.5, 1000)
r = np.zeros((1000,3))
r[:,0] = d
v = np.zeros((1000,1))
f = np.zeros((1000,3))

for i in range(1000):
    v[i] = lj.compute_lj_potential(r[i])
    f[i] = lj.compute_lj_force(r[i])
    
plt.plot(d,v, label="LJ Potential")
plt.plot(d,f[:,0], label="LJ Force")
plt.axis([0.85,2.5,-5,10])
plt.plot(d,np.zeros((1000,1)), 'k')
plt.xlabel("Distance [A.U.]")
plt.ylabel("Magnitude [A.U.]")
plt.grid(True)
plt.legend()
plt.savefig('../fig/ljplot.png')
#plt.show()
