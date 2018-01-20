import numpy as np
from ising_lib import *
import matplotlib.pyplot as plt


def exact_sum(T, states, L):
    en = 0.
    m = 0.
    z = 0.
    for i in range(len(states)):
        ei = energy(states[i], L)
        mi = mag(states[i], L)
        en += ei*np.exp(-ei/T)
        m += mi*np.exp(-ei/T) 
        z += np.exp(-ei/T)
    return en/(L*L*z), m/z

L = 4
states = []
generate(states)
T = np.arange(1.0, 5.1, 0.1)
e0 = []
m0 = []
for t in T:
    e1, m1 = exact_sum(t, states, L)
    e0.append(e1)
    m0.append(m1)

e0 = np.array(e0)
m0 = np.array(m0)

np.savetxt('../dat/exact_en.txt', e0)
np.savetxt('../dat/exact_mag.txt', m0)
