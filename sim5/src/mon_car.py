import numpy as np

def runge( x ):
    return ( 1. + x**2 )**-1
    
def runge_int( x2, x1 ):
    return np.arctan(x2) - np.arctan(x1) 

def simple_sampling( f, a, b, N ):
    x = (b-a)*np.random.random( (N) ) + a
    fx = (b-a)*f(x)
    return np.mean(fx), np.std(fx)/np.sqrt(N)

def metropolis( N, P, trial_move, phi0 ):
    '''
    returns the list of new states and the acceptance rate
    '''
    accept = 0.                 # acceptance rate
    states = [ phi0 ]           # list of all states (first is phi0)
    for k in range( 0, N-1 ):
        newstate = trial_move( states[k] )
        if np.random.random() < np.minimum( 1, P(newstate)/P(states[k]) ):
            states.append( newstate) 
            accept += 1
        else:
            states.append( states[k] )
    return states, accept/(N-1.)

def trial_move_runge( x ):
    global tm_dx
    return x + 2*tm_dx*np.random.random(np.shape(x)) - tm_dx

def trial_move_set_dx( dx ):
    global tm_dx; tm_dx = dx
