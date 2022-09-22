#%%
import numpy as np
import matplotlib.pyplot as plt

import functions

R_a = 8.314 # universal gas constant [J/mol/K]
k_B = 1.380649e-23 # Boltzmann constant [J/K]
N_A = 6.022e23 # Avogadro's number [#/mole]
G = 6.674e-11 # Gravitational constant [m^3/kg/s^2]
M_earth = 5.972e24 # Earth mass [kg]
r_earth = 6371 # Earth radius [km]

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

# total mass density [kg/m^-3]
m_total = (n_N2*28 + n_O2*32 + n_N*14 + n_O*16 + n_Ar*40 + n_He*2 + n_H*1)/1000/N_A

# average "air" molar mass [kg/mole]
M_ave   = (n_N2*28 + n_O2*32 + n_N*14 + n_O*16 + n_Ar*40 + n_He*2 + n_H*1)/n_total/1000

# orbital speed [m/s]
v_sc = np.zeros(len(h))
for i in range(0,len(h)):
    v_sc[i] = np.sqrt(G*M_earth/((r_earth+h[i])*1000))

# mass flow rate [kg/s]
A_sc = 1 # cross-sectional area [m^2]
m_dot = np.zeros(len(h))
V_dot = np.zeros(len(h))
for i in range(0,len(h)):
    m_dot[i] = A_sc*v_sc[i]*n_total[i]*M_ave[i] # [kg/s]
    V_dot[i] = m_dot[i]*R_a*300/(M_ave[i]*N_A)/101325*60*1e6 # [sccm]
    

fig, ax = plt.subplots(figsize=(4.75,6.33))
ax.semilogy(h,m_dot/A_sc)
ax.set_xlim(50,400)
ax.set_xlabel('Altitude, km')
ax.set_ylim(1e-14,1e-5)
ax.set_yticks(np.logspace(-14,-5,10))
ax.set_ylabel('Mass flow rate/A, kg/s/$m^2$')
ax2 = ax.twinx()
ax2.semilogy(h,V_dot)
ax2.set_ylim(1e-6,1e3)
ax2.set_ylabel('Mass flow rate/A, sccm/$m^2$')
ax.grid()
fig.tight_layout()
fig.savefig('Mass flow rate vs Altitude',dpi=300)
plt.show()

#%%
x = 1
print(h[x])
print(v_sc[x])
print(T[x])
print(n_total[x]*1e6)
print(n_N2[x]/n_total[x])
print(n_O2[x]/n_total[x])
print(n_N[x]/n_total[x])
print(n_O[x]/n_total[x])
print(n_Ar[x]/n_total[x])
print(n_He[x]/n_total[x])
print(n_H[x]/n_total[x])
print(V_dot[x])
