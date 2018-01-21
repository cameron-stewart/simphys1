import numpy as np
from ising_lib import *
import matplotlib.pyplot as plt


N = int(1e6)
rho = 0.99005
e = bivariate_gaussian(N, rho)
print(error_analyze(e))
