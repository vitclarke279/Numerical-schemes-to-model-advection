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
    for j in range (1,nx-1):
        phi[j]=phiOld[j] - 0.5*c*(phiOld[j+1] - phiOld[j-1])
    #apply periodic boundary conditions
    phi[0] = phiOld[0] -0.5*c*(phiOld[1] - phiOld[nx-1])
    phi[nx-1] = phi[0]
    
 
    # CTCS for all time steps
    for n in range (1,nt):
        # Code for CTCS at each time-step 
        # Loop over space
        for j in range (1, nx-1):
            phiNew[j] = phiOld[j] - c*(phi[j+1] - phi[j-1])
        # Apply periodic boundary conditions
        phiNew[0] = phiOld[0] - c*(phi[1] - phi[nx-1])
        phiNew[nx-1] = phiNew[0]
        # Update phi for the next time-step
        phiOld = phi.copy()
        phi = phiNew.copy()
        
    return phi
            
        
    
    
