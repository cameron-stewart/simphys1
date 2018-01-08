from __future__ import print_function, division
from numpy import *
from matplotlib.pyplot import *
import scipy.stats
import sys, os
import pickle 

from inspect import currentframe, getframeinfo

# Additional bibs
import argparse

# Own files
from random_numbers import *

# Argparse command line options
parser = argparse.ArgumentParser()

parser.add_argument( '--ID', type=str, help='Simulation ID' )
parser.add_argument( '--dtmax', type=int, help='Maximum of step width (in terms of ts[1]-ts[0])' )
parser.add_argument( '--linreg', type=int, nargs='+', help='--linreg Min Max: area in which it is linear' )

args = parser.parse_args()

### READING DATA-FILE ############################################################################

if not args.ID:
    print('ERROR: simuation ID missing.\nUse --ID <ID>')
    sys.exit(2)
simulation_id = args.ID

datafilename = '../dat/{}.dat'.format(simulation_id)

print("Reading data from {}.".format(datafilename))
datafile = open(datafilename, 'rb')
step, t, x, v, ts, Tms, xs, vs = pickle.load(datafile)
datafile.close()

xs = array(xs)[:,:,0]
vs = array(vs)[:,:,0]
N = shape(xs)[0]
d = shape(xs)[1]
if args.dtmax:
    dtmax = args.dtmax
else:
    dtmax = N-1
if args.linreg:
    linreg = args.linreg
else:
    linreg = [0,dtmax]

### CALCULATE MSD ################################################################################

def calc_msd( ndt, x ):
    '''
    calculates the msd for time steps of dt = ts[ndt] - ts[0]
    '''
    n = shape(x)[0]/ndt
    k = 1
    msd = 0.
    while k < n:
        msd += sum( ( x[ndt*k,:] - x[ndt*(k-1),:] )**2 )
        k += 1
    return msd / (k-1.) 
    
def calc_err( ndt, msd, x ):
    '''
    calculates the error of msd_vec
    '''
    n = shape(x)[0]/ndt
    k = 1
    err = 0.
    while k < n:
        err += ( sum( ( x[ndt*k,:] - x[ndt*(k-1),:] )**2 ) - msd )**2
        k += 1
    if k > 2:
        return sqrt( err / ( (k-1.)*(k-2.) ) )
    else:
        return 0.
    
print('Calculating MSD...')
msd_vec = zeros(dtmax)
dt = zeros(dtmax)
err_vec = zeros(dtmax)
for k in range(1,dtmax+1):
    msd_vec[k-1] = calc_msd( k, xs )
    dt[k-1] = ts[k]-ts[0]
    err_vec[k-1] = calc_err( k, msd_vec[k-1], xs )

### FIT ON DATA ##################################################################################

print('Fit onto MSD...')
p = polyfit( dt[linreg[0]:linreg[1]], msd_vec[linreg[0]:linreg[1]], 1 )
Dx = p[0] * 0.5 / float(d)

# Plot
print('Plot MSD...')
errorbar(dt,msd_vec,yerr=err_vec,fmt='*')
x = linspace(linreg[0],linreg[1],2)
plot( x, p[0]*x+p[1], '-r')
xlabel("Duration $\Delta t$")
ylabel("MSD of $x$")
savefig('../dat/{}_MSD.png'.format(simulation_id))
close()

print('Finished MSD.')

### VACF #########################################################################################

def autocor(v):
    '''
    returns the autocorrelation function of vestorial observable v
    It uses scalar produkt of v with itself
    The result is normed to autocor(v)[0] = 1
    '''
    ftv = real( fft.ifft( conj( fft.fft(v, axis=0) ) * fft.fft(v, axis=0), axis=0 ) )
    while size(ftv.shape) > 1:  # only one loop in this case 
         ftv = sum( ftv, axis=1 ) 
    return ftv[:ftv.shape[0]//2] / ftv[0]

 
# Calculate VACF  
print('Calculate VACF...')
vacf = autocor(vs)

print('Plot VACF...')
plot(arange(0,vacf.shape[0]),vacf)
xlabel("VACF")
ylabel("Time t")
savefig('../dat/{}_VACF.png'.format(simulation_id))
close()

# Integrate VACF
print('Integrate VACF...')
Dv = trapz(vacf)

print('Finished VACF.')

print('\nThe diffusion coefficients are:')
print('Dx = {}'.format(Dx))
print('Dv = {}'.format(Dv))

    


