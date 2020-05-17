# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 16:34:19 2019

@author: vlopa

Script with diagnostic functins to aid the comparison of the models.
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
    """Clculates the total mass of phi"""
    return np.sum(phi) * dx
    
def L2_E_N(nx, c, t):  
    "Function to calculate the L2 error norm of CTCS, CNCS and Semilagrangian" 
    "advection schemes. Function converg loops over different values of nt to"
    "show how the l2 error norm changes as a function of dx and dt."

    #Parameters
    u = 1.0
    #Setting empty arrays for dx,dt and nt   
    dx = np.zeros(len(nx))
    dt = np.zeros(len(nx))
    nt = np.zeros(len(nx), dtype='int')
    #Setting empty arrays for the error norms for all schemes being tested
    L2ErrorNormCTCS = np.zeros(len(nx))
    L2ErrorNormCNCS = np.zeros(len(nx))
    L2ErrorNormSemiLag = np.zeros(len(nx))
        #Loops over all of the space steps
    for i in range (len(nx)):
        dx[i] = 1.0/ (nx[i])
        dt[i] = c[i]*dx[i]/u
        nt[i] = t/dt[i]
        
        # Spatial variable going from zero to one exclusive
        x = np.arange(0, 1, dx[i])
        
        phiOld = initialBell(x)
        
        #Exact solution of the advection
        phiExact = initialBell(x - u*t)
        
        #Advection using CTCS
        [phiCTCS, CTCSmass] = CTCS(phiOld.copy(), c[i], nt[i], dx[i])
        
        L2ErrorNormCTCS[i] = L2ErrorNorm(phiCTCS, phiExact)
        
        
        #Advection using CNCS
        [phiCNCS, CNCSmass] = CNCS(phiOld.copy(), c[i], nt[i], dx[i])
        
        L2ErrorNormCNCS[i] = L2ErrorNorm(phiCNCS, phiExact)
        

        #Advection using Semi-Lagrangian 
        [phiSemiLagrangian, SLmass] = Semi_Lagrangian(phiOld.copy(), c[i], nt[i], dx[i])
        
        L2ErrorNormSemiLag[i] = L2ErrorNorm(phiSemiLagrangian, phiExact)

        
    #Plotting Error norms of advection schemes   
    plt.figure()
    plt.plot(dx, L2ErrorNormCTCS, label='CTCS error', color='blue', marker='*')
    plt.plot(dx, L2ErrorNormCNCS, label='CNCS error', color='green', marker='*')
    plt.plot(dx, L2ErrorNormSemiLag, label='SL error', color='orange', marker='*')
    plt.axhline(0, linestyle=':', color='black')
    #plt.legend(bbox_to_anchor=(1.1,1))
    #plt.title('Graph of L2 Error against the Courant number', fontsize=15)
    plt.xlabel('dx')
    plt.ylabel('$l_{2}$')
    plt.show
    
    #Plotting error norms for the Semi-Lagrangian scheme
    plt.figure()
    plt.plot(dx, L2ErrorNormSemiLag, label='SL error', color='orange', marker='*')
    plt.axhline(0, linestyle=':', color='black')
    plt.xlabel('dx')
    plt.ylabel('$l_{2}$')
    plt.show
    
    #Plotting error norms for the CTCS scheme
    plt.figure()
    plt.plot(dx, L2ErrorNormCTCS, label='CTCS error', color='blue', marker='*')
    plt.axhline(0, linestyle=':', color='black')
    plt.xlabel('dx')
    plt.ylabel('$l_{2}$')
    plt.show
  
    
    y2=4*(dx**2) 
    y3=4*(dx**3)     
    
    #Order of accuracy with dx
    plt.figure(figsize=(10,5))
    plt.loglog(dx, L2ErrorNormCTCS, label='CTCS error', color='blue', marker='*')
    plt.loglog(dx, L2ErrorNormCNCS, label='CNCS error', color='green', marker='*')
    plt.loglog(dx, L2ErrorNormSemiLag, label='SL error', color='orange', marker='*')
    slopeCTCS, interseptCTCS = np.polyfit(np.log(dx), np.log(L2ErrorNormCTCS), 1)
    slopeCNCS, interseptCNCS = np.polyfit(np.log(dx), np.log(L2ErrorNormCNCS), 1)
    slopeSL, interseptSL = np.polyfit(np.log(dx), np.log(L2ErrorNormSemiLag), 1)
    print('CTCS--Slope of graph(order of convergence)=',slopeCTCS)
    print('CNCS--Slope of graph(order of convergence)=',slopeCNCS)
    print('SL--Slope of graph(order of convergence)=',slopeSL)
    plt.axhline(0, linestyle=':', color='black')
    plt.loglog(dx, y2, linestyle=':', color='pink')
    plt.loglog(dx, y3, linestyle=':', color='red')
    plt.title('Log Log graph of L2 Error against dx', fontsize=12)
    plt.xlabel('$dx$')
    plt.ylabel('L2 Error', fontsize=12)
    
    
    
    #Order of accuracy with c
    plt.figure(figsize=(10,5))
    plt.loglog(c, L2ErrorNormCTCS, label='CTCS error', color='blue', marker='*')
    plt.loglog(c, L2ErrorNormCNCS, label='CNCS error', color='green', marker='*')
    plt.loglog(c, L2ErrorNormSemiLag, label='SL error', color='orange', marker='*')
    slopeCTCS, interseptCTCS = np.polyfit(np.log(c), np.log(L2ErrorNormCTCS), 1)
    slopeCNCS, interseptCNCS = np.polyfit(np.log(c), np.log(L2ErrorNormCNCS), 1)
    slopeSL, interseptSL = np.polyfit(np.log(c), np.log(L2ErrorNormSemiLag), 1)
    print('CTCS--Slope of graph(order of convergence)=',slopeCTCS)
    print('CNCS--Slope of graph(order of convergence)=',slopeCNCS)
    print('SL--Slope of graph(order of convergence)=',slopeSL)
    plt.axhline(0, linestyle=':', color='black')
    plt.title('Log Log graph of L2 Error against c', fontsize=12)
    plt.xlabel('$c$')
    plt.ylabel('L2 Error', fontsize=12)  
    
