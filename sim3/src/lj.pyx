import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern void c_set_globals(double _L, int _N, double _rcut, double _shift)
cdef extern void c_compute_forces(double *x, double *f)
cdef extern double c_compute_energy (double *x, double *v, double *E_pot, double *E_kin)
cdef extern double c_rebuild_neighbor_lists(double *x, double vlsize)
cdef extern double c_compute_pressure(double Ekin, double* x)
cdef extern void c_velocity_rescaling(double T0, double T, double* v)
cdef extern double c_force_capping(double* f, double fcap)

def set_globals(double L, int N, double rcut, double shift):
    c_set_globals(L, N, rcut, shift)

def compute_forces(np.ndarray[double, ndim=2, mode='c'] x not None):
    cdef np.ndarray[double, ndim=2, mode='c'] f = np.zeros_like(x)
    N = x.shape[1]
    c_compute_forces (&x[0,0], &f[0,0])
    return f

def compute_energy(np.ndarray[double, ndim=2, mode='c'] x not None, 
                   np.ndarray[double, ndim=2, mode='c'] v not None):
    N = x.shape[1]
    cdef double E_tot, E_kin, E_pot
    E_tot = c_compute_energy(&x[0,0], &v[0,0], &E_pot, &E_kin)
    return E_pot, E_kin, E_tot

def rebuild_neighbor_lists(np.ndarray[double, ndim=2, mode='c'] x not None,
                           double vlsize):
    N = x.shape[1]
    c_rebuild_neighbor_lists(&x[0,0], vlsize)
    
def compute_pressure(double E_kin, 
                  np.ndarray[double, ndim=2, mode='c'] x not None):
    cdef double P
    N = x.shape[1]
    P = c_compute_pressure(E_kin, &x[0,0])
    return P

def velocity_rescaling(T0, T, np.ndarray[double, ndim=2, mode='c'] v not None):
    N = v.shape[1]
    c_velocity_rescaling(T0, T, &v[0,0])
    
def force_capping(np.ndarray[double, ndim=2, mode='c'] f not None, fcap):
    N = f.shape[1]
    fcap = c_force_capping(&f[0,0], fcap)
