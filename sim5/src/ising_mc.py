import numpy as np
from ising_lib import *
import matplotlib.pyplot as plt

def ising_mc(T, L):
    e = []
    m = []
    state = initial(L)
    for i in range(10000):
        e.append(energy(state, L))
        m.append(mag(state, L))
        montecarlo(T, state, L)
    e = np.array(e)
    m = np.array(m)
    return e/(L*L), m

L = 4
T = np.arange(1.0, 5.1, 0.1)
e0 = []
m0 = []
for t in T:
    e1, m1 = ising_mc(t, L)
    e0.append(e1)
    m0.append(m1)

e0 = np.array(e0)
m0 = np.array(m0)

np.savetxt('../dat/mc_en_series.txt', e0)
np.savetxt('../dat/mc_mag_series.txt', m0)

