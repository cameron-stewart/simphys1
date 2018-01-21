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

plt.plot(T, e0)
plt.plot(T, e2)
plt.figure()
plt.plot(T, m0)
plt.plot(T, m2)
plt.show()

