# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 12:53:25 2020

@author: vlopa
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


nx = np.array([100,100,100,100])
c = np.array([0.2,0.5,1.0,1.2])
u = 1.0
t = 3

L2_E_N(nx, c, t)


nx1 = nx[0]
c1 = c[0]

main(nx1, c1, t)

nx2 = nx[1]
c2 = c[1]

main(nx2, c2, t)

nx3 = nx[2]
c3 = c[2]

main(nx3, c3, t)

nx4 = nx[3]
c4 = c[3]

main(nx4, c4, t)

