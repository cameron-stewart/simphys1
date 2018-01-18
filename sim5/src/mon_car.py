import numpy as np

def simple_sampling( f, a, b, N):
    x = (b-a)*np.random.random( (N) ) + a
    fx = np.sum( f(x), axis=0 )
    return fx * (b-a)/N

def runge( x ):
    return ( 1. + x**2 )**-1
    

    
