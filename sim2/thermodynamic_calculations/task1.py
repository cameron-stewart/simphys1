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
