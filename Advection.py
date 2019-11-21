# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 16:30:51 2019

@author: vlopa
"""

import matplotlib.pyplot as plt

import numpy as np

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from Advection_schemes import CTCS
#from Advection_schemes import CNCS
#from Diagnostics import
from Initial_conditions import initialBell



def main():
    "Advect a bell curve on a domain between x = xmin and x = xmax over nx "
    "spatial steps with advection coefficient c, time step dt for nt time steps"
    
    # Parameters
    nx = 100
    c = 0.9
    u = 1.0
    nt=111
    
    
    # Derived parameters
    dx = 1.0/nx
    dt = c*dx/u
    t = nt*dt
    print ("Courant Number =" , c)
    print ("dx =" , dx, "dt =", dt, "nt = ",nt)
    print ("end time =" ,t)
    
    # Spatial variable going from zero to one inclusive
    x = np.linspace (0.0, 1.0, nx+1)
    print ('x=', x)
    
    # Three time levels of the dependant variable, phi
    phiOld = initialBell(x)
    
    #Advection using CTCS, CNCS and Semi-Langrangian
    phiCTCS = CTCS(phiOld.copy(), c, nt)
    #phiCNCS = CNCS()
    #phiSL = Semi_Langrangian
    
    #Calculate and print out error norms
    #
    #
    
    
    font = {'size' :20}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, initialBell(x - u*t), label='Analytic' , color='black', linestyle= '--', linewidth=2)
    plt.plot (x, phiCTCS, label='CTCS', color='blue')
    #plt.plot(x, phiBTCS, label='BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim(0,1.2)
    plt.legend(bbox_to_anchor=(1.1, 1))
    plt.xlabel ('$x$')
    
main()    
    
    
    

    