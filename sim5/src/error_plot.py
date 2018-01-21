import numpy as np
import matplotlib.pyplot as plt
from ising_lib import *

e0 = np.loadtxt('../dat/exact_en.txt')
e1 = np.loadtxt('../dat/mc_en_series.txt')
m0 = np.loadtxt('../dat/exact_mag.txt')
m1 = np.loadtxt('../dat/mc_mag_series.txt')
T = np.arange(1.0, 5.1, 0.1)

e2 = np.empty(len(T))
e2_err = np.empty(len(T))
m2 = np.empty(len(T))
m2_err = np.empty(len(T))

for i in range(len(T)):
    e_temp = error_analyze(e1[i,:])
    e2[i] = e_temp[0]
    e2_err[i] = e_temp[1]
    m_temp = error_analyze(m1[i,:])
    m2[i] = m_temp[0]
    m2_err[i] = m_temp[1]

plt.plot(T, e0, label="Exact")
plt.errorbar(T, e2, yerr=e2_err, fmt='o', label='Monte Carlo')
plt.legend(loc='best')
plt.xlabel('Temperature T')
plt.ylabel('Energy Density')
plt.savefig('../fig/errorplot_e.png', bbox_inches='tight')
plt.figure()
plt.plot(T, m0, label='Exact')
plt.errorbar(T, m2, yerr=m2_err, fmt='o', label ='Monte Carlo')
plt.legend(loc='best')
plt.xlabel('Temperature T')
plt.ylabel('Magnetization Density')
plt.savefig('../fig/errorplot_m.png', bbox_inches='tight')
