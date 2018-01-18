import numpy as np

def runge( x ):
    return ( 1. + x**2 )**-1
    
def runge_int( x2, x1 ):
    return np.arctan(x2) - np.arctan(x1) 

def simple_sampling( f, a, b, N):
    x = (b-a)*np.random.random( (N) ) + a
    fx = (b-a)*f(x)
    return np.mean(fx), np.std(fx)/np.sqrt(N)
    
