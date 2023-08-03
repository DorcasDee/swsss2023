# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 09:43:37 2023

@author: USER
"""

import numpy as np
import matplotlib.pyplot as plt

x0 = 3.0
# step size
h = 0.2
#horizon
T = 2

# Numbrt of steps needed
N = int(np.ceil(T/h))


def f(x):
    return -2*x

def analytical_sol(t):
    return 3*np.exp(-2*t)

xs = np.zeros(N+1)
ts = np.zeros(N+1)

xs[0] = x0
ts[0] = 0.0

for i in range(N):
    xs[i+1] = xs[i] + h*f(xs[i]) # f(xs[i]) = dxs[i]/dt
    ts[i+1] = ts[i] + h
    
fig, ax = plt.subplots()
ax.plot(ts, xs, 'k.', label = 'Explicit Euler')

trange = np.linspace(0, T, 100)
ax.plot(trange, analytical_sol(trange), linewidth=2, color = 'red', label = 'Analytical Euler')

legend = ax.legend(loc='upper center', fontsize='x-large')