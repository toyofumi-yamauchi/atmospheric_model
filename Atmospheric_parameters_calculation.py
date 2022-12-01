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
ax.semilogy(H,n_total,'-',label='total',linewidth=3)
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
l1 = ax.semilogy(H,n_total,'-',label='total',linewidth=3)
l2 = ax.semilogy(H,n_N2,'--',label='N2')
l3 = ax.semilogy(H,n_O2,'--',label='O2')
#ax.semilogy(H,n_N,'--',label='N')
l4 = ax.semilogy(H,n_O,'--',label='O')
#ax.semilogy(H,n_Ar,'--',label='Ar')
#ax.semilogy(H,n_He,'--',label='He')
#ax.semilogy(H,n_H,'--',label='H')
ax.axhspan(1e14,1e15,color='k',alpha=0.2,label='WVC operation limit')
ax.set_xlim(50,250)
ax.set_xlabel('Altitude, km')
ax.set_ylim(1e14,1e23)
ax.set_yticks(np.logspace(14,23,10))
ax.set_ylabel('Number density, m$^-$$^3$')
ls = l1+l2+l3+l4
labs = [l.get_label() for l in ls]
ax.legend(ls,labs,framealpha=1.0)
ax.grid()
ax.set_title('Density vs Altitude')
fig.tight_layout()
fig.savefig('Density vs Altitude (Closeup)',dpi=300)
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

#%%
M = np.array([14.0067*2, 15.999*2, 15.999]) # molar mass for N2, O2, O [amu]
m = M/1000/N_A                              # mass for N2, O2, O [kg]
fig, ax = plt.subplots(figsize=(9.5,6.33))
l1 = ax.plot(H,T,label='temperature',linewidth=3.0)
ax.set_xlim(50,400)
ax.set_xlabel('Altitude, km')
ax.set_ylim(150,850)
ax.set_yticks(np.linspace(150,850,15))
ax.set_ylabel('Temperature, K')
ax2 = ax.twinx()
l2 = ax2.semilogy(H,v_sc,'k',label='orbital speed')
ax2.plot([0,0],[0,0])
l3 = ax2.semilogy(H,f.v_th(T,m[0]),label='thermal speed: N2')
l4 = ax2.semilogy(H,f.v_th(T,m[1]),label='thermal speed: O2')
l5 = ax2.semilogy(H,f.v_th(T,m[2]),label='thermal speed: O')
ax2.set_ylim(100,10000)
ax2.set_ylabel('Speed, m/s')
ax3 = ax.twinx()
ax3.plot([0,0],[0,0])
l6 = ax3.plot(H,np.arctan(f.v_th(T,m[0])/v_sc)*180/np.pi,'--',label='max $\\theta$: N2')
l7 = ax3.plot(H,np.arctan(f.v_th(T,m[1])/v_sc)*180/np.pi,'--',label='max $\\theta$: O2')
l8 = ax3.plot(H,np.arctan(f.v_th(T,m[2])/v_sc)*180/np.pi,'--',label='max $\\theta$: O')
ax3.set_ylim(0,45)
ax3.set_ylabel('Maximum angle in flow, $\degree$')
ax3.spines["right"].set_position(("axes",1.25))
ls = l1+l2+l3+l4+l5+l6+l7+l8
labs = [l.get_label() for l in ls]
ax3.legend(ls,labs,framealpha=1.0,loc='lower left', bbox_to_anchor=(1.45, 0.65))
ax.grid()
ax.set_title('Speed, Temperature vs Altitude')
fig.tight_layout()
fig.savefig('Speed, Temperature vs Altitude',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(6.0,4.25))
l1 = ax.semilogy(H,n_total*v_sc,label='total',linewidth=3)
l2 = ax.plot(H,n_N2*v_sc,'--',label='N2')
l3 = ax.plot(H,n_O2*v_sc,'--',label='O2')
l4 = ax.plot(H,n_O*v_sc,'--',label='O')
ax.set_xlim(50,400)
ax.set_xlabel('Altitude, km')
ax.set_ylim(1e18,1e26)
ax.set_yticks(np.logspace(16,26,11))
ax.set_ylabel('Flux, #/$m^2$s')
ls = l1+l2+l3+l4
labs = [l.get_label() for l in ls]
ax.legend(ls,labs,framealpha=1.0,loc='best')
ax.grid()
ax.set_title('Flux vs Altitude')
fig.tight_layout()
fig.savefig('Flux vs Altitude',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=figure_size)
l1 = ax.plot(H,(n_N2+n_O2+n_O)/n_total*100,label='%(N2+O2+O)')
l2 = ax.plot(H,n_N2/n_total*100,label='%N2')
l3 = ax.plot(H,n_O2/n_total*100,label='%O2')
l4 = ax.plot(H,n_O/n_total*100,label='%O')
ax.set_xlim(50,400)
ax.set_xlabel('Altitude, km')
ax.set_ylim(0,100)
ax.set_yticks(np.linspace(0,100,11))
ax.set_ylabel('Fraction, %')
ls = l1+l2+l3+l4
labs = [l.get_label() for l in ls]
ax.legend(ls,labs,framealpha=1.0,loc='best')
ax.grid()
ax.set_title('N2,O2,O vs Altitude')
fig.tight_layout()
fig.savefig('N2,O2,O vs Altitude',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(4.75,4.75))
l1 = ax.plot([0,0],[0,0])
l2 = ax.plot([0,0],[0,0])
l3 = ax.plot(H,n_O2/(n_O2+n_O)*100,label='%O2')
l4 = ax.plot(H,n_O/(n_O2+n_O)*100,label='%O')
ax.set_xlim(50,400)
ax.set_xlabel('Altitude, km')
ax.set_ylim(0,100)
ax.set_yticks(np.linspace(0,100,11))
ax.set_ylabel('Fraction, %')
ls = l1+l2+l3+l4
labs = [l.get_label() for l in ls]
ax.legend(ls,labs,framealpha=1.0,loc='best')
ax.grid()
ax.set_title('O2,O vs Altitude')
fig.tight_layout()
fig.savefig('O2,O vs Altitude',dpi=300)
plt.show()
#%% 
h = 400 # km
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
#%% 
n_x = 1.00e19
h_x = f.H_from_n(n_x,n_N2 + n_O2 + n_O,H)
print('h = {:.1f} km for n = {:.2e} m^-3'.format(h_x,n_x))
#%%
Γ_x = 0.625e21          # #/m^2/s
I_x = Γ_x*1.602e-19 # A/m^2
h_x = f.H_from_n(Γ_x,(n_N2 + n_O2 + n_O)*7800,H)
print('h = {:.1f} km for Γ = {:.2e} m^-3'.format(h_x,Γ_x))
print('I_x = {:.1f} mA/cm^2'.format(I_x*1000/100/100))