import numpy as np
import matplotlib.pyplot as plt

def calc_gravity(alt):
    ''' 
    calculates gravity at changing altitudes
    '''
    r = 6370 # radius, kilometres
    g = (3.99 * 10**14)/((1000 * (r + alt))**2) # gravity array
    return g

def calc_scaleheight(m, alt, T):
    '''
    calculates scale height from inputs mass, altitude and temperature
    '''
    k = 1.38 * 10**(-23) # boltzmann constant
    H = (k * T)/(m * calc_gravity(alt)) # scale height array
    return H
    
def build_hydro(m, n_0, T, alt):
    '''
    input:  mass and initial density of the element, temperature and altitude
    
    routput: density array
    ''' 
    dz = (alt[1]-alt[0]) * 1000 # converts from km to m
    
    density = []
    density.insert(0, n_0) # insert first value of the densities
    
    for i in np.arange(1, len(alt)):
        # starting from density[1], density[0] set above
        H = calc_scaleheight(m, alt, T)
        dens = T[i-1]/T[i] * density[i-1] * np.exp(- dz/H[i-1])
        density.append(dens)
    return density


def plot_of_density(dens_array, name_of_elements, alt):
    ''' 
    creates a plot of density returned from previous function
    '''
    fig, ax = plt.subplots(len(name_of_elements), 1, figsize=(8, 6), sharex = False)
    for i in range(len(name_of_elements)):
        
        ax[i].semilogx(dens_array[i], alt)
        ax[i].set_title('Plot of Altitude against {} Density'.format(name_of_elements[i]))
        # ax[i].set_xlabel('{} Density [m$^-3)$]'.format(name_of_elements[i]))
        ax[i].set_ylabel('Altitude [km]')



m = 1.67 * 10**(-27) # atomic mass unit (amu)
mN2 = 28*m # nitrogen mass, kilograms
mO2 = 32*m # oxygen2
mOx = 16*m # oxygen
N2_n_0 = 1 * 10**19 # density of nitrogen
O2_n_0 = 0.3 * N2_n_0 # density of oxygen
Ox_n_0 = 1 * 10**18 # density of oxygen


# T = np.linspace(200, 1000, nPts) # temperature array
# nPts = 100 # number of points
# alt = np.linspace(100, 500, nPts) # altitude array

# N2_dens = build_hydro(mN2, N2_n_0, T)
# O2_dens = build_hydro(mO2, O2_n_0, T)
# Ox_dens = build_hydro(mOx, Ox_n_0, T)

# dens_array = [N2_dens, O2_dens, Ox_dens]
# name_of_elements = ['Nitrogen', 'Oxygen 2', 'Oxygen']

# plot_of_density(dens_array, name_of_elements)