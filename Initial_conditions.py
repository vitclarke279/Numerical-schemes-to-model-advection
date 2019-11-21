# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 16:35:55 2019

@author: vlopa
"""
import numpy as np

def initialBell(x):
    return np.where(x%1. < 0.5, np.power (np.sin(2*x*np.pi), 2), 0)
