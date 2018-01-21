import numpy as np
import matplotlib.pyplot as plt

e0 = np.loadtxt('../dat/exact_en.txt')
e1 = np.loadtxt('../dat/mc_en_series.txt')
m0 = np.loadtxt('../dat/exact_mag.txt')
m1 = np.loadtxt('../dat/mc_mag_series.txt')
T = np.arange(1.0, 5.1, 0.1)

e2 = np.empty(len(T))
m2 = np.empty(len(T))
for i in range(len(T)):
    e2[i] = np.mean(e1[i,:])
    m2[i] = np.mean(m1[i,:])

plt.plot(T, e0, label='Exact')
plt.plot(T, e2, 'o', label='Monte Carlo')
plt.legend(loc='best')
plt.xlabel('Temperature T')
plt.ylabel('Energy Density e')
plt.tight_layout()
plt.savefig('../fig/isingplot_e.png', bbox_inches='tight')
plt.figure()
plt.plot(T, m0, label='Exact')
plt.plot(T, m2, 'o', label='Monte Carlo')
plt.legend(loc='best')
plt.xlabel('Temperature T')
plt.ylabel('Magnetization Density m')
plt.tight_layout()
plt.savefig('../fig/isingplot_m.png', bbox_inches='tight')
