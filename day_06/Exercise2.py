# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 15:26:11 2023

@author: USER
"""

# FIRST ATTEMPT

import numpy as np
import matplotlib.pyplot as plt

r = 6370 # radius, kilometres
n_0 = 1 * 10**19 # density
nPts = 100 # number of points
m = 28 * 1.67 * 10**(-27) # mass, kilograms
k = 1.38 * 10**(-23) # boltzmann constant
T = np.linspace(200, 1000, nPts) # temperature array
alt = np.linspace(100, 500, nPts) # altitude array
g = (3.99 * 10**14)/((1000 * (r + alt))**2) # gravity array
H = (k * T)/(m * g) # scale height array

# for i in np.arange(len(alt)):
#     Temp = T[i] + T
#     H = (k * T)/(m * g) # scale height array
dz = (alt[1]-alt[0]) * 1000
n = []
n.insert(0, n_0)
n.insert(1, T[0]/T[1] * n_0 * np.exp(- dz/H[0]))
print(n)

for i in np.arange(2, len(alt)):
    # starting from n[2]
    # n[0] and n[1] set above
    dens = T[i-1]/T[i] * n[i-1] * np.exp(- dz/H[i-1])
    n.append(dens)

fig, ax = plt.subplots(1, figsize=(8, 6))
ax.plot(alt, n)
plt.xscale('log')
ax.set_title('Plot of Altitude against Density')
ax.set_xlabel('Density')
ax.set_ylabel('Altitude')

# height = (k * T[0])
print('H[0]', H[0])
print('H[-1]', H[-1])
print('dz', dz)
print('n[1]', n[1])
print('n[-1]', n[-1])
print('g[0]', g[0])

#%%

# SECOND ATTEMPT



import numpy as np
import matplotlib.pyplot as plt

r = 6370 # radius, kilometres
n_0 = 1 * 10**19 # density
nPts = 100 # number of points
m = 28 * 1.67 * 10**(-27) # mass, kilograms
k = 1.38 * 10**(-23) # boltzmann constant
T = np.linspace(200, 1000, nPts) # temperature array
alt = np.linspace(100, 500, nPts) # altitude array
g = (3.99 * 10**14)/((1000 * (r + alt))**2) # gravity array
H = (k * T)/(m * g) # scale height array


dz = (alt[1]-alt[0]) * 1000 # converts from km to m
n = []
n.insert(0, n_0)
# Formula for n_1 is =  T[0]/T[1] * n_0 * np.exp(- dz/H[0]))
print(n)

for i in np.arange(1, len(alt)):
    # starting from n[1], n[0] set above
    # T[i-1] and n[i-1] corresponds to T_0 and n_0 respectively
    dens = T[i-1]/T[i] * n[i-1] * np.exp(- dz/H[i-1])
    n.append(dens)

fig, ax = plt.subplots(1, figsize=(8, 6))
plt.semilogx(n, alt)
ax.set_title('Plot of Altitude against Nitrogen Density')
ax.set_xlabel('Nitrogen Density [m$^-3)$]')
ax.set_ylabel('Altitude [km]')

# height = (k * T[0])
print('H[0]', H[0])
print('H[-1]', H[-1])
print('dz', dz)
print('n[1]', n[1])
print('n[-1]', n[-1])
print('g[0]', g[0])

#%%
print('I ate {} for breakfast'.format('beans'))
#%%
# Add oxygen 2 and oxygen 1

import numpy as np
import matplotlib.pyplot as plt

nPts = 100 # number of points
m = 1.67 * 10**(-27) # atomic mass unit (amu)
alt = np.linspace(100, 500, nPts) # altitude array

mN2 = 28*m # nitrogen mass, kilograms
mO2 = 32*m # oxygen2
mOx = 16*m # oxygen
N2_n_0 = 1 * 10**19 # density of nitrogen
O2_n_0 = 0.3 * N2_n_0 # density of oxygen
Ox_n_0 = 1 * 10**18 # density of oxygen


def dens_of_elements(m, n_0):
    '''
    reads in mass and initial density of the element
    
    returns the density array
    
    contains constants r, k, T, g, H, dz
    ''' 
    r = 6370 # radius, kilometres
    k = 1.38 * 10**(-23) # boltzmann constant
    T = np.linspace(200, 1000, nPts) # temperature array
    g = (3.99 * 10**14)/((1000 * (r + alt))**2) # gravity array
    H = (k * T)/(m * g) # scale height array
    dz = (alt[1]-alt[0]) * 1000 # converts from km to m
    
    density = []
    density.insert(0, n_0) # insert first value of the densities
    
    for i in np.arange(1, len(alt)):
        # starting from density[1], density[0] set above
        dens = T[i-1]/T[i] * density[i-1] * np.exp(- dz/H[i-1])
        density.append(dens)
    return density

def plot_of_density(density, name_of_element):
    ''' 
    creates a plot of density returned from previous function
    '''
    fig, ax = plt.subplots(1, figsize=(8, 6))
    ax.semilogx(density, alt)
    ax.set_title('Plot of Altitude against {} Density'.format(name_of_element))
    ax.set_xlabel('{} Density [m$^-3)$]'.format(name_of_element))
    ax.set_ylabel('Altitude [km]')
    
N2_dens = dens_of_elements(mN2, N2_n_0)
plot_of_density(N2_dens, 'Nitrogen')

O2_dens = dens_of_elements(mO2, O2_n_0)
plot_of_density(O2_dens, 'Oxygen 2')

Ox_dens = dens_of_elements(mOx, Ox_n_0)
plot_of_density(Ox_dens, 'Oxygen')

#%%
# Add oxygen 2 and oxygen 1

import numpy as np
import matplotlib.pyplot as plt

nPts = 100 # number of points
m = 1.67 * 10**(-27) # atomic mass unit (amu)
alt = np.linspace(100, 500, nPts) # altitude array

mN2 = 28*m # nitrogen mass, kilograms
mO2 = 32*m # oxygen2
mOx = 16*m # oxygen
N2_n_0 = 1 * 10**19 # density of nitrogen
O2_n_0 = 0.3 * N2_n_0 # density of oxygen
Ox_n_0 = 1 * 10**18 # density of oxygen


T = np.linspace(200, 1000, nPts) # temperature array
def dens_of_elements(m, n_0, T):
    '''
    reads in mass and initial density of the element
    
    returns the density array
    
    contains constants r, k, T, g, H, dz
    ''' 
    r = 6370 # radius, kilometres
    k = 1.38 * 10**(-23) # boltzmann constant
    g = (3.99 * 10**14)/((1000 * (r + alt))**2) # gravity array
    H = (k * T)/(m * g) # scale height array
    dz = (alt[1]-alt[0]) * 1000 # converts from km to m
    
    density = []
    density.insert(0, n_0) # insert first value of the densities
    
    for i in np.arange(1, len(alt)):
        # starting from density[1], density[0] set above
        dens = T[i-1]/T[i] * density[i-1] * np.exp(- dz/H[i-1])
        density.append(dens)
    return density

def plot_of_density(dens_array, name_of_elements):
    ''' 
    creates a plot of density returned from previous function
    '''
    fig, ax = plt.subplots(len(name_of_elements), 1, figsize=(8, 6), sharex = False)
    for i in range(len(name_of_elements)):
        
        ax[i].semilogx(dens_array[i], alt)
        ax[i].set_title('Plot of Altitude against {} Density'.format(name_of_elements[i]))
        # ax[i].set_xlabel('{} Density [m$^-3)$]'.format(name_of_elements[i]))
        ax[i].set_ylabel('Altitude [km]')

N2_dens = dens_of_elements(mN2, N2_n_0)
O2_dens = dens_of_elements(mO2, O2_n_0)
Ox_dens = dens_of_elements(mOx, Ox_n_0)

dens_array = [N2_dens, O2_dens, Ox_dens]
name_of_elements = ['Nitrogen', 'Oxygen 2', 'Oxygen']

plot_of_density(dens_array, name_of_elements)