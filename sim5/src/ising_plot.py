import numpy as np
import matplotlib.pyplot as plt

e0 = np.loadtxt('../dat/exact_en.txt')
e1 = np.loadtxt('../dat/mc_en.txt')
m0 = np.loadtxt('../dat/exact_mag.txt')
m1 = np.loadtxt('../dat/mc_mag.txt')
T = np.arange(1.0, 5.1, 0.1)

plt.plot(T, e0)
plt.plot(T, e1)
plt.figure()
plt.plot(T, m0)
plt.plot(T, m1)
plt.show()

