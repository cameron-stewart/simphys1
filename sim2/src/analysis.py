import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b):
    return a*x*x+b*x

N = np.power(np.array([5,6,7,8,9,10,11,12]),3)
# Don't judge me
t = np.array([0.188, 0.340, 0.592, 1.124, 1.988, 3.376, 5.772, 9.300])

x = np.linspace(0,12*12*12,100)

popt, pcov = curve_fit(func, N, t)

plt.plot(N,t, 'o', label="data")
plt.plot(x, func(x,popt[0], popt[1]), label="quadratic fit")
plt.legend(loc=0)
plt.xlabel("$N = n^3$")
plt.ylabel("Time [s]")
plt.savefig("../fig/fit.png")
