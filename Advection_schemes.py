# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 16:26:30 2019

@author: vlopa
"""

import numpy as np

import numpy.linalg as la

def CTCS(phiOld, c, nt):
    "Advection of profile in phiOld using CTCS non-dimentional Courant"
    "number, c"
    
    nx=len(phiOld)
    
    # New time-step array for phi
    phi = phiOld.copy()
    phiNew = phiOld.copy()
    
    
    #FTCS for the first time-step, looping over space
    for j in range(nx):
        phi[j]=phiOld[j] - 0.5*c*(phiOld[(j+1)%nx] - phiOld[(j-1)%nx])
 
    # CTCS for all time steps
    for n in range (1,nt):
        # Code for CTCS at each time-step 
        # Loop over space
        for j in range(nx):
            phiNew[j] = phiOld[j] - c*(phi[(j+1)%nx] - phi[(j-1)%nx])
        # Update phi for the next time-step
        phiOld = phi.copy()
        phi = phiNew.copy()
        
    return phi

def BTCS(phi, c, nt):
    nx = len(phi)
    
    M = np.zeros([nx,nx])
    for j in range(nx):
        M[j,j] = 1
        M[(j-1)%nx,j] = 0.5*c
        M[(j+1)%nx,j] = -0.5*c

    # Solution for nt time steps
    for it in range(nt):
        phi = np.linalg.solve(M, phi) 
    
    return phi

def CNCS (phi, c, nt):
    
    #Array for the RHS of the matrix equation
    RHS = phi.copy()
    
    nx = len(phi)
    
    M = np.zeros([nx,nx])
    for j in range (nx):
        M[j,j] = 1
        M[(j-1)%nx,j] = 0.25*c
        M[(j+1)%nx,j] = -0.25*c
        
    #Solution for nt time steps
    for it in range(nt):
        for j in range (nx):
            RHS[j] = phi[j] - 0.25*c*(phi[(j+1)%nx] - phi[(j-1)%nx])
        
        phi = np.linalg.solve(M, RHS)
        
    return phi
