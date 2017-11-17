"""
A library containing functions to compute the LJ potential and force

Author: Cameron Stewart
"""

import numpy as np
import numpy.linalg as la

def compute_lj_potential(r_ij):
    """Compute LJ potential on particle j from particle i"""
    d = la.norm(r_ij)
    return 4*(d**-12 - d**-6)

def compute_lj_force(r_ij):
    """Compute LJ force on particle j due to particle i"""
    d = la.norm(r_ij)
    f = 48*(d**-13 - d**-7)
    return f*r_ij/d
