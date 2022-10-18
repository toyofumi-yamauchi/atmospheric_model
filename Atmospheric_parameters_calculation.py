#%%
import numpy as np
import matplotlib.pyplot as plt

from functions import load_NRLMSIS_data, P_from_n_and_T, neutral_flux
import functions as f

R   = 8.314 # universal gas constant [J/mol/K]
k_B = 1.380649e-23 # Boltzmann constant [J/K]
N_A = 6.022e23 # Avogadro's number [#/mole]
G = 6.674e-11 # Gravitational constant [m^3/kg/s^2]
M_earth = 5.972e24 # Earth mass [kg]
r_earth = 6371 # Earth radius [km]

filename='msis20output_81312020-2.txt'
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

# total number density [m^-3]
n_total = n_N2 + n_O2 + n_N + n_O + n_Ar + n_He + n_H

# total mass density [kg/m^-3]
m_total = (n_N2*28 + n_O2*32 + n_N*14 + n_O*16 + n_Ar*40 + n_He*2 + n_H*1)/1000/N_A

# average "air" molar mass [kg/mol]
M_ave   = (n_N2*28 + n_O2*32 + n_N*14 + n_O*16 + n_Ar*40 + n_He*2 + n_H*1)/n_total/1000

# average "air" mass [kg]
m_air   = M_ave/N_A

# atmospheric parameter
A_sc = 1.0 # cross-sectional area [m^2]
P = np.zeros(len(H))     # pressure [Pa]
Γ = np.zeros(len(H))     # neutral flux [#/m^2/s]
v_sc = np.zeros(len(H))  # orbital speed [m/s]
m_dot = np.zeros(len(H)) # mass flow rate [sccm]
V_dot = np.zeros(len(H)) # volumetric flow rate [sccm]
for i in range(0,len(H)):
    P[i] = f.P_from_n_and_T(n_total[i],T[i])
    Γ[i], v_sc[i] = f.neutral_flux(H[i],n_total[i])
    m_dot[i] = m_air[i]*Γ[i]*A_sc
    V_dot[i] = m_dot[i]/m_air[i]/N_A*(R*300/101325)*60*1e6

# plots
figure_size = (4.75,6.33)

fig, ax = plt.subplots(figsize=figure_size)
l1 = ax.semilogy(H,m_dot/A_sc,label='$\dot{m}$')
ax.set_xlim(50,400)
ax.set_xlabel('Altitude, km')
ax.set_ylim(1e-8,1e1)
ax.set_yticks(np.logspace(-8,1,10))
ax.set_ylabel('Mass flow rate/A, kg/s/$m^2$')
ax2 = ax.twinx()
l2 = ax2.semilogy(H,V_dot,'--',label='$\dot{V}$')
ax2.set_ylim(1e0,1e9)
ax2.set_yticks(np.logspace(0,9,10))
ax2.set_ylabel('Volmetric flow rate/A, sccm/$m^2$')
ax.grid()
ax.set_title('Flow rate vs Altitude')
lns = l1+l2
labs = [l.get_label() for l in lns]
ax.legend(lns,labs,framealpha=1.0)
fig.tight_layout()
fig.savefig('Mass flow rate vs Altitude',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=figure_size)
ax.semilogy(H,n_total,'-',label='total')
ax.semilogy(H,n_N2,'--',label='N2')
ax.semilogy(H,n_O2,'--',label='O2')
ax.semilogy(H,n_N,'--',label='N')
ax.semilogy(H,n_O,'--',label='O')
ax.semilogy(H,n_Ar,'--',label='Ar')
ax.semilogy(H,n_He,'--',label='He')
ax.semilogy(H,n_H,'--',label='H')
ax.axhspan(1e14,1e15,color='k',alpha=0.2,label='WVC operation limit')
ax.set_xlim(50,400)
ax.set_xlabel('Altitude, km')
ax.set_ylim(1e14,1e23)
ax.set_yticks(np.logspace(14,23,10))
ax.set_ylabel('Number density, m$^-$$^3$')
ax.grid()
ax.legend(framealpha=1.0)
ax.set_title('Density vs Altitude')
fig.tight_layout()
fig.savefig('Density vs Altitude',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=figure_size)
ax.semilogy(H,P)
ax.set_xlim(50,400)
ax.set_xlabel('Altitude, km')
ax.set_ylim(1e-7,1e2)
ax.set_yticks(np.logspace(-7,2,10))
ax.set_ylabel('Pressure, Pa')
ax2 = ax.twinx()
ax2.semilogy(H,P*760/101325)
ax2.axhspan(1e-10,1e-8,color='k',alpha=0.2,label='WVC operation limit')
ax2.set_ylim(1e-7*760/101325,1e2*760/101325)
ax2.set_yticks(np.logspace(-9,0,10))
ax2.set_ylabel('Pressure, Torr')
ax.grid()
ax.set_title('Pressure vs Altitude')
ax2.legend(framealpha=1.0)
fig.tight_layout()
fig.savefig('Pressure vs Altitude',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=figure_size)
ax.plot(H,T)
ax.set_xlim(50,400)
ax.set_xlabel('Altitude, km')
ax.set_ylim(150,850)
ax.set_yticks(np.linspace(150,850,15))
ax.set_ylabel('Temperature, K')
ax.grid()
ax.set_title('Temperature vs Altitude')
fig.tight_layout()
fig.savefig('Temperature vs Altitude',dpi=300)
plt.show()
#%%
h = 105 # km
x = int(f.index_from_H(h,H))
print('h       = {} km'.format(H[x]))
print('v_sc    = {:.2f} m/s'.format(v_sc[x]))
print('T       = {:.1f} K'.format(T[x]))
print('P       = {:.1e} Pa'.format(P[x]))
print('n_total = {:.2e} m^-3'.format(n_total[x]))
print('n_N2    = {:.2e} m^-3 ({:.1f}%)'.format(n_N2[x],n_N2[x]/n_total[x]*100))
print('n_O2    = {:.2e} m^-3 ({:.1f}%)'.format(n_O2[x],n_O2[x]/n_total[x]*100))
print('n_N     = {:.2e} m^-3 ({:.1f}%)'.format(n_N[x],n_N[x]/n_total[x]*100))
print('n_O     = {:.2e} m^-3 ({:.1f}%)'.format(n_O[x],n_O[x]/n_total[x]*100))
print('n_Ar    = {:.2e} m^-3 ({:.1f}%)'.format(n_Ar[x],n_Ar[x]/n_total[x]*100))
print('n_He    = {:.2e} m^-3 ({:.1f}%)'.format(n_He[x],n_He[x]/n_total[x]*100))
print('n_H     = {:.2e} m^-3 ({:.1f}%)'.format(n_H[x],n_H[x]/n_total[x]*100))
print('Γ       = {:.2e} 1/m^2/s'.format(Γ[x]))
print('m_dot/A = {:.2e} kg/s/m^2'.format(m_dot[x]))
print('V_dot/A = {:.2e} sccm/m^2'.format(V_dot[x]))
