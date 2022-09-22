#%%
import numpy as np
import matplotlib.pyplot as plt

import functions

k_B = 1.380649e-23 # Boltzmann constant [J/K]

filename='msis20output_81312020-2.txt'
data = functions.load_NRLMSIS_data(filename)

h = data[:,5]     # altitude [km]
n_N2 = data[:,9]  # N2 number density [cm^-3]
n_O2 = data[:,10] # O2 number density [cm^-3]
n_N = data[:,17]  # N  number density [cm^-3]
n_O = data[:,8]   # O  number density [cm^-3]
n_Ar = data[:,15] # Ar number density [cm^-3]
n_He = data[:,14] # He number density [cm^-3]
n_H = data[:,16]  # H  number density [cm^-3]
T = data[:,12]    # temperature [K]

# total number density [cm^-3]
n_total = n_N2 + n_O2 + n_N + n_O + n_Ar + n_He + n_H

# pressure [Pa]
P = np.zeros(len(T))
for i in range(0,len(T)):
    P[i] = n_total[i]*1e6*k_B*T[i]

# Plot
fig, ax = plt.subplots(figsize=(4.75,6.33))
ax.semilogy(h,P)
ax.set_xlim(50,400)
ax.set_xlabel('Altitude, km')
ax.set_ylim(1e-7,1e2)
ax.set_yticks(np.logspace(-7,2,10))
ax.set_ylabel('Pressure, Pa')
ax2 = ax.twinx()
ax2.semilogy(h,P*760/101325)
ax2.set_ylim(1e-7*760/101325,1e2*760/101325)
ax2.set_yticks(np.logspace(-9,0,10))
ax2.set_ylabel('Pressure, Torr')
ax.grid()
fig.tight_layout()
fig.savefig('Pressure vs Altitude',dpi=300)
plt.show()
#%%
x = 6
print(h[x])
print(T[x])
print(P[x])
print(P[x]*760/101325)