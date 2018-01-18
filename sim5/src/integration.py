import numpy as np
from mon_car import *
import matplotlib.pyplot as plt

I_est = np.empty(19)
std_err = np.empty(19)
act_err = np.empty(19)
N_range = np.empty(19)

for i in range(2,21):  
    I_est[i-2], std_err[i-2] = simple_sampling( runge, -5., 5., 2**i )
    act_err[i-2] = np.absolute(I_est[i-2] - runge_int(5, -5))
    N_range[i-2] = 2**i

plt.loglog(N_range, std_err, label='Statistical Error', basex=2)
plt.loglog(N_range, act_err, 'r^', label='Actual Error', basex=2)
plt.legend()
plt.xlabel('Number of Samples N')
plt.ylabel('Error')
plt.savefig('../fig/simple_err.png', bbox_inches='tight')
