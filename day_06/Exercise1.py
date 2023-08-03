# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:37:17 2023

@author: USER
"""

#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 10, 10)
y = x**3 + 2*x**2 - x

plt.plot(x, y)

# first derivative
dydx = 3*x**2 + 4*x - 1
plt.plot(x, dydx)

d2ydx2 = 6*x + 4
plt.plot(x, d2ydx2)

print(x)
# print(y)
# print(dydx)
# print(d2ydx2)

#%%
def eqn_fn(x):
    y = x**3 + 2*x**2 - x
    return y

def first_derivative(x):
    y = x**3 + 2*x**2 - x
    return y

# eqn_fn(x[2]) - eqn_fn(x[0]) / 2 * (x[1] - x[0])

dfdx = []
for i in np.arange(1,len(x)-1):
    num = eqn_fn(x[i+1]) - eqn_fn(x[i-1])
    den = 2 * (x[i+1] - x[i])
    
    dfdx.append(num/den)
dfdx.insert(0, dfdx[0])
dfdx.insert(-1, dfdx[-1])

print(dfdx)
plt.plot(x, dydx)
plt.plot(x, dfdx)
#%%
d2fdx2 = []
for i in np.arange(1, len(x)-1):
    num = eqn_fn(x[i+1]) + eqn_fn(x[i-1]) - 2 * eqn_fn(x[i])
    den = (x[i+1] - x[i])**2
    d2fdx2.append(num/den)

d2fdx2.insert(0, d2fdx2[0])
d2fdx2.insert(-1, d2fdx2[-1])

plt.plot(x, d2ydx2)
plt.plot(x, d2fdx2)
    