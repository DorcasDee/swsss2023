# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 16:14:23 2023

@author: USER
"""
# Question: find x such that 0 = f(x)

from scipy.optimize import fsolve
import numpy as np
# right hand side
def f(x):
    return np.exp(x) - 4*x
# initial guess
x0 = 2.0
# solve statement assings solution to x
x = fsolve(f, x0)

print(x)

# has two roots, different values for a positive and negative initial guess