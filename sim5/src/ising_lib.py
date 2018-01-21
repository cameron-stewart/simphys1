import numpy as np
import matplotlib.pyplot as plt

def initial(L):
    # Generates a random initial state of side length L
    state = np.zeros([L,L], dtype=int)

    for i in range(L):
        for j in range(L):
            if np.random.randint(0,2):
                state[i,j]=1
            else:
                state[i,j]=-1
    return state

def energy(state, L):
    # Returns the energy of a given state of side length L
    E = 0.
    for i in range(L):
        for j in range(L):
            E += -state[i,j]*(state[(i-1)%L,j] + state[(i+1)%L,j] + state[i,(j+1)%L] 
                    + state[i,(j-1)%L])
    return E/2

def mag(state, L):
    # Returns the absolute value of the magnetization per site for a state of side length L
    mu = 0.
    for i in range(L):
        for j in range(L):
            mu += state[i,j]
    return np.abs(mu/(L*L))

def montecarlo(T, state, L):
    # Performs a single Monte Carlo step on a state with side length L using the 
    # Metropolis-Hastings algorithm
    i = np.random.randint(0,L)
    j = np.random.randint(0,L)
    delta = 2*state[i,j]*(state[(i-1)%L,j] + state[(i+1)%L,j] + state[i,(j+1)%L] 
                + state[i,(j-1)%L])
    if delta < 0:
        state[i,j] *= -1
    elif np.random.rand() < np.exp(-delta/T):
        state[i,j] *= -1

def generate(states):
    # Generates all possible states with a side length L = 4. Don't judge me for this
    # monstrosity.
    state = np.zeros([4,4])
    for i in range(2**16):
        binstr = "{0:b}".format(i)
        for j in range(16-len(binstr)):
            binstr = '0' + binstr
        for k in range(len(binstr)):
            if binstr[k]=='1':
                state[k%4,k/4]=binstr[k]
            else:
                state[k%4,k/4]=-1
        states.append(np.copy(state))

def bivariate_gaussian(N, rho):
    # generates N samples from a bivariate gaussian
    e = np.random.normal(size=N)
    for i in range(1, N):
        e[i] = rho*e[i-1] + np.sqrt(1-rho**2)*e[i]
    return e

def error_analyze(e):
    # Given a time series of observables, returns the mean, size of bins, autocorrelation time,
    # effective number of samples, and measured error of the mean using binning analysis
    k_final = 0
    tau_final = 0
    N_eff = 0
    error = 0
    N = len(e)
    
    for k in range(1,N/10):
        if (N%k==0):
            Nb = N/k
            varB = 0 
            Ob = np.zeros(Nb)
            for n in range(Nb):
                for i in range(k):
                    Ob[n] += e[n*k+i]/k
            Ob_mean = np.mean(Ob)
            for n in range(Nb):
                varB += (Ob[n] - Ob_mean)**2/(Nb-1)
            tau = k*varB/(2*np.var(e))
            if tau > tau_final:
                tau_final = tau
                k_final = k
                N_eff = N/(2*tau)   
                error = np.sqrt(varB/Nb)

    return Ob_mean, error, k_final, tau_final, N_eff

