from numpy import *
from matplotlib.pyplot import *
import sys, os

def init_LCG(X):
	'''
	Initialize the random nuber generator with Xlgc.
	'''
	global Xlcg
	Xlcg = X
	
def LCG( m=2**32, a=1664525, c=1013904223 ):
	'''
	Returns a normal distributed random number in [0,m-1] (natural number).
	'''
	global Xlcg
	Xlcg = ( a*Xlcg + c ) % m
	return Xlcg

def normal_LCG( m=2**32, a=1664525, c=1013904223 ):
	'''
	Returns normal ditributed random numbers in [0,1]
	'''
	return LCG( m, a, c ) / ( m - 1 )
	
def do2int(x):
	'''
	takes a double and converts it into an integer by changing the place of the decimal dot behind the last relevant place.
	EXAMPLE: 123.456 -> int(123456)
	'''
	while x % 1:
		x *= 10
	return int(x)
	
def random_walk( N=1000 ):
	'''
	returns a random walk for N steps with an deviation in [-0.5,0.5]
	'''
	x = zeros(N+1)
	for k in range(0,N):
		x[k+1] = x[k] + normal_LCG() - 0.5
	
	# Plot
	plot(range(0,N+1),x,label=r'random walk')
	xlabel('time')
	ylabel('position')
	show()
	
	return x
	
def calc_BM( u1, u2 ):
	'''
	converts uniform random numbers into normal distributet random numbers
	'''
	return sqrt(-2.*log(u1)) * array([cos(2.*pi*u2), sin(2.*pi*u2)])
	
def BM( N, mu=2.0, sigma=5.0 ):
	'''
	returns a numpy array of N normal ditributed random numbers
	'''
	n = int(ceil(N/2.))
	bm = zeros((n,2))
	for k in range(0,n):
		bm[k,:] = sigma*calc_BM(random.random(), random.random()) + mu
	return bm.flatten()
	
def gauss( x, mu=2.0, sigma=5.0 ):
	'''
	returns the gaussian-function
	'''
	return exp(-0.5*((x-mu)/sigma)**2) / (sqrt(2.0*pi)*sigma)

def show_BM_hist( N=10000, mu=2.0, sigma=5.0, limits=[-15,20], bars=150 ):
	'''
	draws gaussian and hist of random numbers with BM
	'''
	bins=linspace(limits[0], limits[1], bars)
	
	BM_numbers = BM(N, mu, sigma)
	
	x = linspace(limits[0], limits[1], 1000)
	gauss_fun = gauss(x, mu, sigma)
	
	plot(x, gauss_fun, label=r'gauss: $\mu={} \sigma={}$'.format(mu,sigma))
	hist(BM_numbers, bins, normed=1, label=r'BM    : $\mu={} \sigma={}$'.format(mu,sigma))
	legend()
	show()
	
def rand_vel_vec( N=1000, sigma=1. ):
	'''
	returns an numpy array with shape (N,3) of normal distributed random numbers
	'''
	vel_vec = zeros((N,3))
	for k in range(0,3):
		vel_vec[:,k] = BM(N, mu=0., sigma=sigma)
	return vel_vec
	
def gauss_3d( r, sigma=1. ):
	'''
	returns the 3D gaussian-function
	'''
	return exp(-0.5*(r/sigma)**2) * r**2 * sqrt(2./pi) * sigma**-6
	
def show_vel_hist( N=1000, sigma=1., limits=[-0,5], bars=100 ):
	'''
	draws gaussian and hist of the velocity distribution
	'''
	bins=linspace(limits[0], limits[1], bars)
	
	vel_abs = linalg.norm(rand_vel_vec(N, sigma), axis=1)
	
	r = linspace(limits[0], limits[1], 1000)
	gauss_fun = gauss_3d(r, sigma)
	
	plot(r, gauss_fun, label=r'gauss: $\sigma={}$'.format(sigma))
	hist(vel_abs, bins, normed=1, label=r'BM    : $\sigma={}$'.format(sigma))
	legend()
	show()
