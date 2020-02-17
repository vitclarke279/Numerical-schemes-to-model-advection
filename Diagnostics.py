# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 16:34:19 2019

@author: vlopa
"""
from scipy.integrate import quad
import numpy as np

def L2ErrorNorm(phi, phiExact):
    """Calculates the L2 error norm (RMS error) of phi in comparison to
    phiExact, ignoring the boundaries"""
    
    # calculate the error and the error norms
    phiError = phi - phiExact
    L2 = np.sqrt(sum(phiError**2)/sum(phiExact**2))

    return L2

def mass(phi, dx):
    return np.sum(phi) * dx
    

    
