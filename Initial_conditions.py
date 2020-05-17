# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 16:35:55 2019

@author: vlopa

Function with the initial conditions entered into the advection models.
"""
import numpy as np

def initialBell(x):
    "Normal distribution of phi as initial condition"
    return np.where(x%1. < 0.5, np.power (np.sin(2*x*np.pi), 2), 0)
