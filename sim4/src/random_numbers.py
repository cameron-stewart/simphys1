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
	
def calc_BM( u1, u2 ):
	n = sqrt(-2*log(u1)) * [cos(2*pi*u2),sin(2*pi*u2)]
