import numpy as np
import matplotlib.pyplot as plt
from tridiagonal import solve_tridiagonal
import density as ds

# ----------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------

if __name__ == "__main__":

    dx = 4.0

    # set x with 1 ghost cell on both sides:
    x = np.arange(100-dx, 500 + 2 * dx, dx)

    t_upper = 1000.0

    nPts = len(x)

    # set default coefficients for the solver:
    a = np.zeros(nPts) - 1
    b = np.zeros(nPts) + 2
    c = np.zeros(nPts) - 1
    d = np.zeros(nPts)
    
    # set Q background
    Q = np.zeros(nPts)
    Q[np.logical_and((x<400), (x>200))] = 0.4
    lambda_ = 80.0
    dz = x[1] - x[0]
    dzZ = dz * dz
    # set Q of EUV
    Q_euv = np.zeros(nPts)
    
    # introducing changing time, no longer fixzed local time
    n_days = 27
    # hours:
    dt = 1
    # hours:
    times = np.arange(0, n_days*24, dt)
    lon = 0.0 # longitude
    
    # an empty temp 2D array
    temp_2D = np.zeros((nPts, len(times)))
    times_2D = np.zeros((nPts, len(times)))
    alts_2D = np.zeros((nPts, len(times)))
    N2_2D = np.zeros((nPts, len(times)))
    O2_2D = np.zeros((nPts, len(times)))
    Ox_2D = np.zeros((nPts, len(times)))
    
    # set amplitude and phase
    Amp_Diurnal = 10.0
    Amp_Semi_Diurnal = 5.0
    Phase_Diurnal = np.pi/2
    Phase_Semi_Diurnal = 3 * np.pi/2
    
    # set F10.7
    F10_7 = 100 + 50/(24*365)*times + 25*np.sin(times/(27*24)*2*np.pi)
    
    # plot
    nrows = 2; ncols = 2
    fig, ax = plt.subplots(nrows, ncols, figsize = (12,12))
    
    for i, hour in enumerate(times):
        ut = hour % 24
        local_time = lon/15 + ut
        
        fac = -np.cos(local_time/24 * 2 * np.pi)
        if fac < 0:
            fac = 0
        
        # Add a source term:
        Sun_Heat = F10_7[i] * 0.4/100
        Q_euv[np.logical_and((x<400), (x>200))] = Sun_Heat * fac # Qeuv is a function of fac which is a function of time
        Q_back = Q
        d = (Q_back + Q_euv) * dzZ/lambda_
        
        t_lower = 200.0 + Amp_Diurnal * np.sin(local_time/24 * 2*np.pi + Phase_Diurnal) \
                        + Amp_Semi_Diurnal * np.sin(local_time/24 * 2*np.pi *2 + Phase_Semi_Diurnal)
        
        # boundary conditions (bottom - fixed):
        a[0] = 0
        b[0] = 1
        c[0] = 0
        d[0] = t_lower
    
        # top - floating boundary conditions:
        a[-1] = 1
        b[-1] = -1
        c[-1] = 0
        d[-1] = 0
    
        # solve for Temperature:
        t = solve_tridiagonal(a, b, c, d)
        
        densN2 = ds.build_hydro(ds.mN2, ds.N2_n_0, t, x)
        densO2 = ds.build_hydro(ds.mO2, ds.O2_n_0, t, x)
        densOx = ds.build_hydro(ds.mOx, ds.Ox_n_0, t, x)
        
        temp_2D[:,i] = t
        times_2D[:,i] = hour / 24
        alts_2D[:,i]= x
        N2_2D[:,i] = densN2
        O2_2D[:,i] = densO2
        Ox_2D[:,i] = densOx
    
    # plots
    z_array = [temp_2D, np.log10(N2_2D), np.log10(O2_2D), np.log10(Ox_2D)]
    cbar_names = ['Temperature', '$N_2$ Density', '$O_2$ Density', '$O$ Density']
    for a, z, cbar_name in zip(ax.flatten(), z_array, cbar_names):
        cmap = a.contourf(times_2D, alts_2D, z)
        cbar = fig.colorbar(cmap, ax=a)
        cbar.ax.set_ylabel(cbar_name, fontsize=13)
        a.set_ylabel('Altitude (km)', fontsize=15)
        a.set_xlabel('Local Time (Days)', fontsize=15)
        a.set_title('Plot of Time, Altitude and {}'.format(cbar_name), fontsize=15)
        
    plt.tight_layout()
    
    plotfile = 'conduction_v6.png'
    print('writing : ',plotfile)    
    fig.savefig(plotfile)
    # plt.close()
    