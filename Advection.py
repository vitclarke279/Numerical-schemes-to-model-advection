# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 16:30:51 2019

@author: vlopa

Script running advection models.
"""

import matplotlib.pyplot as plt

import numpy as np

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from Advection_schemes import *
from Diagnostics import *
from Initial_conditions import *

nx = 100
c = 0.5
t = 3

def main(nx, c, t):
    "Advect a bell curve on a domain between x = xmin and x = xmax over nx "
    "spatial steps with advection coefficient c, time step dt for nt time steps"
    
    # Parameters
    u = 1.0
    
    #derived parameters
    dx = 1.0/nx
    dt = c*dx/u
    nt = int(t/dt)
    
    print ("Courant Number =" , c)
    print ("dx =" , dx, "dt =", dt, "nt = ",nt)
    print ("end time =" ,t)
    
    # Spatial variable going from zero to one inclusive
    x = np.arange(0,1,dx)
    print ('x=', x)
    
    time=np.arange(0,t+dt,dt)
    
    # Three time levels of the dependant variable, phi
    phiOld = initialBell(x)
    
    #Advection using CTCS, CNCS and Semi-Langrangian
    [phiCTCS, CTCSmass] = CTCS(phiOld.copy(), c, nt, dx)
    [phiCNCS, CNCSmass] = CNCS(phiOld.copy(), c, nt, dx)
    [phiSL, SLmass] = Semi_Lagrangian(phiOld.copy(), c, nt, dx)
    
    
    #phiExact
    phiExact = initialBell(x - u*t)
    
    
    #Calculate and print out error norms
    CTCS_L2_error = L2ErrorNorm(phiCTCS, phiExact)
    print ('L2 Error Norm for CTCS:',CTCS_L2_error)
    CNCS_L2_error = L2ErrorNorm(phiCNCS, phiExact)
    print ('L2 Error Norm for CNCS:', CNCS_L2_error)
    SL_L2_error = L2ErrorNorm(phiSL, phiExact)
    print ('L2 Error Norm for Semi-Lagrangian:' , SL_L2_error)
    
    #Boundedness of the schemes
    print('Exact solution phi max=',max(phiExact))
    print('Exact solution phi min=',min(phiExact))
    print('CTCS phi max=',max(phiCTCS))
    print('CTCS phi min=', min(phiCTCS))
    print('CNCS phi max=',max(phiCNCS))
    print('CNCS phi min=', min(phiCNCS))
    print('SL phi max=',max(phiSL))
    print('SL phi min=', min(phiSL))    
   
    #Plotting all advection scheme
    plt.figure()
    plt.figure(figsize=(10,5))
    font = {'size' :20}
    plt.rc('font', **font)
    #Plotting starting phi
    #plt.plot(x, phiOld, label= 'Initial', color= 'black', linestyle= ':')
    #Plotting exact solution
    plt.plot(x, phiExact, label='Analytic' , color='black', linewidth=2)
    #plotting advection schemes
    plt.plot (x, phiCTCS, label='CTCS', color='orange')
    plt.plot (x, phiCNCS, label='CNCS', color= 'green')
    plt.plot(x, phiSL, label='Semi Lagrangian', color='blue')
    plt.ylim(-0.1,1.1)
    plt.xlim(0,1)
    plt.legend(loc='best', fontsize=20)
    plt.title(f'c={c},nx={nx}', fontsize=22)
    plt.xlabel ('$x$', fontsize=20)
    plt.ylabel('$\phi$', fontsize=20)
    plt.show
    
    #plotting CTCS scheme
    plt.figure()
    font = {'size' :20}
    plt.rc('font', **font)
    #Plotting starting phi
    plt.plot(x, phiOld, label= 'Initial', color= 'black', linestyle= ':')
    #Plotting exact solution
    plt.plot(x, phiExact, label='Analytic' , color='gray', linestyle= '--', linewidth=2)
    #Plotting scheme
    plt.plot (x, phiCTCS, label='CTCS', color='blue')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim(-0.1,1.1)
    plt.xlim(0,1)
    plt.title(f'c={c},nx={nx}', fontsize=15)
    plt.xlabel ('$x$')
    plt.ylabel('$\phi$')
    plt.show
    
    
    #Plotting Crank-Nicolson scheme
    plt.figure()
    font = {'size' :20}
    plt.rc('font', **font)
    #Plotting starting phi
    plt.plot(x, phiOld, label= 'Initial', color= 'black', linestyle= ':')
    #Plotting exact solution
    plt.plot(x, phiExact, label='Analytic' , color='gray', linestyle= '--', linewidth=2)
    #Plotting scheme
    plt.plot (x, phiCNCS, label='CNCS', color= 'green')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim(-0.1,1.1)
    plt.xlim(0,1)
    plt.title(f'c={c},nx={nx}', fontsize=15)
    plt.xlabel ('$x$')
    plt.ylabel('$\phi$')
    plt.show
    
    
    #Plotting Semi-LAgrangian scheme with cubic interpolation
    plt.figure()
    font = {'size' :20}
    plt.rc('font', **font)
    #Plotting starting phi
    plt.plot(x, phiOld, label= 'Initial', color= 'black', linestyle= ':')
    #Plotting exact solution
    plt.plot(x, phiExact, label='Analytic' , color='gray', linestyle= '--', linewidth=2)
    #Plotting scheme
    plt.plot(x, phiSL, label='Semi Lagrangian', color='orange')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim(-0.1,1.1)
    plt.xlim(0,1)
    plt.title(f'c={c},nx={nx}', fontsize=15)
    plt.xlabel ('$x$')
    plt.ylabel('$\phi$')
    plt.show
    
    
    #Plotting conservation of phi for all advection schemes used
    plt.figure()
    plt.plot(time, CNCSmass, label='CNCS mass of $\phi$', color='green')
    plt.plot(time, CTCSmass, label='CTCS mass of $\phi$', color='blue')
    plt.plot(time, SLmass, label='SL mass of $\phi$', color='orange')
    plt.grid(b=None, which='major', axis='both')
    plt.legend(loc='best', fontsize=12)
    #plt.title(f'c={c},nx={nx}', fontsize=15)
    plt.title("Conservation of the schemes")
    plt.ylim(0.249,0.251)
    plt.xlabel('$s$')
    plt.ylabel('mass of $\phi$')
    plt.show
    
    
    #Plotting conservation of phi for CTCS 
    plt.figure()
    plt.plot(time, CTCSmass, label='CTCS mass of $\phi$', color='blue')
    plt.axhline(0, linestyle=':', color='black')
    plt.grid(b=None, which='major', axis='both')
    plt.title(f'c={c},nx={nx}', fontsize=15)
    plt.ylim(0.249,0.251)
    plt.xlabel('$s$')
    plt.ylabel('mass of $\phi$')
    plt.show
    
    
    #Plotting conservation of phi for CNCS
    plt.figure()
    plt.plot(time, CNCSmass, label='CNCS mass of $\phi$', color='green')
    plt.axhline(0, linestyle=':', color='black')
    plt.grid(b=None, which='major', axis='both')
    plt.title(f'c={c},nx={nx}', fontsize=15)
    plt.ylim(0.249,0.251)
    plt.xlabel('$s$')
    plt.ylabel('mass of $\phi$')
    plt.show
    
    
    #Plotting conservation of phi for Semi Lagrangian
    plt.figure()
    plt.plot(time, SLmass, label='SL mass of $\phi$', color='orange')
    plt.axhline(0, linestyle=':', color='black')
    plt.grid(b=None, which='major', axis='both')
    plt.title(f'c={c},nx={nx}', fontsize=15)
    plt.ylim(0.249,0.251)
    plt.xlabel('$s$')
    plt.ylabel('mass of $\phi$')
    plt.show
    

  
    

    
    
    
    
    
    
    
    
    
    
    

    
