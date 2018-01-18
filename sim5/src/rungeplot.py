import numpy as np
from mon_car import *
import matplotlib.pyplot as plt

x = np.linspace(-5,5,100)
y = runge(x)

plt.plot(x,y)
plt.savefig('../fig/rungeplot.png',bbox_inches='tight')
