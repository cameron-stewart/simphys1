#!/usr/bin/python
import numpy as np

# 1.1 calculate the entropy
# partition function
omega1 = 31 # omega = exponent of partition function 10**x
omega2 = 28
omega = omega1+omega2
print 'omega = {}'.format(omega)

# entropy s: # maybe as num*kb
kb = 1.38*10**-23 #J/K
s1 = kb*omega1*np.log(10)
s2 = kb*omega2*np.log(10)
s = kb*omega*np.log(10)
print 's1 = {}\ns2 = {}\ns = {}\n'.format(s1,s2,s)

# isothermal:
mol = 6.022*10**23
fV = 1.01
p0 = 101325. #Pa
V0 = 1. #m^3
T = 298. #K
fO = fV**(p0*V0/(T))
N = p0*V0/(kb*T)/mol
print 'fO = {}^1/kB\nN = {}\n'.format(fO,N)

# isochoric:
N2 = mol
T0 = 400 #K
dU = 100000 #J
fT = 1.+2.*dU/(3.*kb*N2*T0)
power = (3.*N2/2.)
print 'fT = {}\nfO = {}^{}'.format(fT,fT,power)
