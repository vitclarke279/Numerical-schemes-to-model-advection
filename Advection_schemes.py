# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 16:26:30 2019

@author: vlopa
"""

import numpy as np
import math as m
from Diagnostics import *

import numpy.linalg as la

def CTCS(phiOld, c, nt, dx):
    "Advection of profile in phiOld using CTCS non-dimentional Courant"
    "number, c"
    
    nx=len(phiOld)
    
    # New time-step array for phi
    phi = phiOld.copy()
    phiNew = phiOld.copy()
    
    totalMass = np.zeros(nt+1)
    totalMass[0] = mass(phi, dx)
    
    #FTCS for the first time-step, looping over space
    for j in range(nx):
        phi[j]=phiOld[j] - 0.5*c*(phiOld[(j+1)%nx] - phiOld[(j-1)%nx])
        
    totalMass[1] = mass(phi, dx)
        
    # CTCS for all time steps
    for n in range (1,nt):
        # Code for CTCS at each time-step 
        # Loop over space
        for j in range(nx):
            phiNew[j] = phiOld[j] - c*(phi[(j+1)%nx] - phi[(j-1)%nx])
        # Update phi for the next time-step
        phiOld = phi.copy()
        phi = phiNew.copy()
        totalMass[n+1] = mass(phi, dx)
    
    #L2 error norm
    #L2ErrorNorm(phi, phiExact)
    
    return phi, totalMass

def BTCS(phi, c, nt, dx):
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

def CNCS (phi, c, nt, dx):
    
    #Array for the RHS of the matrix equation
    RHS = phi.copy()
    
    nx = len(phi)
    
    totalMass = np.zeros(nt+1)
    totalMass[0] =mass(phi, dx)
    
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
        totalMass[it+1] = mass(phi, dx)
        
    return phi, totalMass



def Semi_Lagrangian (phi, c, nt, dx):
    
    nx = len(phi)
    
    # New timestep array for phi
    phiNew = phi.copy()
    
    totalMass = np.zeros(nt+1)
    totalMass[0] = mass(phi, dx)
    
    for n in range (nt):
        
        for j in range(nx):
            k = m.floor(j-c)
            beta = j-k-c
            phiNew[j] = -((1/6)*beta*(1-beta)*(2-beta)*phi[(k-1)%nx]) \
                      + (0.5*(1+beta)*(1-beta)*(2-beta)*phi[k%nx])\
                      + (0.5*(1+beta)*beta*(2-beta)*phi[(k+1)%nx])\
                      - ((1/6)*(1+beta)*beta*(1-beta)*phi[(k+2)%nx])
        phi = phiNew.copy()
        totalMass[n+1] = mass(phi, dx)
    
    return phi, totalMass

def Semi_Lagrangian_mass (phi, c, nt, dx):
    
    nx = len(phi)
    
    # New timestep array for phi
    phiNew = phi.copy()
    
    totalMass = np.zeros(nt+1)
    totalMass[0] = mass(phi, dx)
    
    for n in range (nt):
        
        for j in range(nx):
            k = m.floor(j-c)
            beta = j-k-c
            phiNew[j] = -((1/6)*beta*(1-beta)*(2-beta)*phi[(k-1)%nx]) \
                      + (0.5*(1+beta)*(1-beta)*(2-beta)*phi[k%nx])\
                      + (0.5*(1+beta)*beta*(2-beta)*phi[(k+1)%nx])\
                      - ((1/6)*(1+beta)*beta*(1-beta)*phi[(k+2)%nx])
        phi = phiNew.copy()
        totalMass[n+1] = mass(phi, dx)
    
    return totalMass











