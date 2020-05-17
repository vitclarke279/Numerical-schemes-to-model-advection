# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 12:53:25 2020

@author: vlopa

Main script running all the advection models with different parameters and 
producing graphs to to perform a comparison.
"""
import matplotlib.pyplot as plt

import numpy as np

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from Advection_schemes import *
from Diagnostics import *
from Initial_conditions import *
from Advection import*
from L2_Error_Norm import*


#nx = np.array([500,250,166,125,100,50,34,25])
#c = np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5])
nx = np.array([100])
c = np.array([10.5])
u = 1.0
t = 3

#L2_E_N(nx, c, t)


nx1 = nx[0]
c1 = c[0]

main(nx1, c1, t)

nx2 = nx[1]
c2 = c[1]

#main(nx2, c2, t)

nx3 = nx[2]
c3 = c[2]

#main(nx3, c3, t)

nx4 = nx[3]
c4 = c[3]

#main(nx4, c4, t)

nx5 = nx[4]
c5 = c[4]

#main(nx5, c5, t)

nx6 = nx[5]
c6 = c[5]

#main(nx6, c6, t)

#nx7 = nx[6]
#c7 = c[6]

#main(nx7, c7, t)

#nx8 = nx[7]
#c8 = c[7]

#main(nx8, c8, t)