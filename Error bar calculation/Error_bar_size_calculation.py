#%%
import numpy as np
import matplotlib.pyplot as plt

import functions as f

R   = 8.314 # universal gas constant [J/mol/K]
k_B = 1.380649e-23 # Boltzmann constant [J/K]
N_A = 6.022e23 # Avogadro's number [#/mole]
G = 6.674e-11 # Gravitational constant [m^3/kg/s^2]
M_earth = 5.972e24 # Earth mass [kg]
r_earth = 6371 # Earth radius [km]

# https://kauai.ccmc.gsfc.nasa.gov/instantrun/msis
filename='50 - Year.txt'
filename='50 - Day of Year.txt'
filename='50 - Hour of Day.txt'
filename='50 - Latitude.txt'
filename='50 - Longitude.txt'
data = f.load_NRLMSIS_data(filename)

Type = np.array(('Year','Day of Year','Hour of Day','Latitude','Longitude'))

for i in range(5):
    print('========')
    print('Type: '+Type[i])
    filename = '100 - '+Type[i]+'.txt'
    data = f.load_NRLMSIS_data(filename)

    H = data[:,5]         # altitude [km]
    n_N2 = data[:,9]*1e6  # N2 number density [m^-3]
    n_O2 = data[:,10]*1e6 # O2 number density [m^-3]
    n_N = data[:,17]*1e6  # N  number density [m^-3]
    n_O = data[:,8]*1e6   # O  number density [m^-3]
    n_Ar = data[:,15]*1e6 # Ar number density [m^-3]
    n_He = data[:,14]*1e6 # He number density [m^-3]
    n_H = data[:,16]*1e6  # H  number density [m^-3]
    T = data[:,12]        # temperature [K]

    gas = np.array(('N2','O2','O','N','Ar','He','H'))
    n_gas = np.array((n_N2,n_O2,n_O,n_N,n_Ar,n_He,n_H))
    for j in range(3):
        print(' -------')
        print(' Gas = '+gas[j])
        print('  max: {:.2e} m^-3'.format(max(n_gas[j])))
        print('  ave: {:.2e} m^-3'.format(np.average(n_gas[j])))
        print('  min: {:.2e} m^-3'.format(min(n_gas[j])))