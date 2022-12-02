def load_NRLMSIS_data(filename):
    '''
    Converting the .txt file into the array
    .txt file is downloaded from https://kauai.ccmc.gsfc.nasa.gov/instantrun/msis
    '''

    # Import libraries
    import numpy as np
    import re

    # Importing the .txt file
    with open(filename,'r') as f:
        lines = f.readlines()
    #print(lines[0])
    #print(lines[1])

    # Converting the lines into rows
    num_rows = np.shape(lines)[0]-1
    num_cols = np.shape(re.split(' +',lines[1]))[0]
    rows = np.zeros((num_rows,num_cols))
    for i in range(0,np.shape(lines)[0]-1):
        rows[i,:] = re.split(' +',lines[i+1])
        #print(rows[i,:])

    data = rows

    return data

def P_from_n_and_T(n,T):
    '''
    Output: 
    P = pressure [Pa]
    Input: 
    n = total number density [m^-3]
    T = gas temperature [K]
    '''
    k_B = 1.380649e-23 # Boltzmann constant [J/K]
    P = n*k_B*T        # Pressure [Pa]
    return P

def neutral_flux(h,n):
    '''
    Output:
    v_sc = orbital speed [m/s]
    Γ    = neutral flux [#/m^2/s]
    Input: 
    h = altitude [kg]
    n = total number density [m^-3]
    '''
    import numpy as np

    G = 6.674e-11      # Gravitational constant [m^3/kg/s^2]
    M_earth = 5.972e24 # Earth mass [kg]
    r_earth = 6371     # Earth radius [km]

    v_sc = np.sqrt(G*M_earth/((r_earth+h)*1000))
    Γ = n*v_sc
    return Γ, v_sc

def index_from_H(h,H):
    '''
    Output: 
    x = index for a given h
    Input:
    h = the altitude [km]
    H = the list of altitude [km]
    '''
    ΔH = H[1]-H[0]
    x = (h-H[0])/ΔH
    return x

def H_from_n(n_x,N,H):
    '''
    Output: 
    h_x = altitude giving n_x [km]
    Input:
    n_x = the arbitrary density [m^-3]
    N = the list of density [m^-3]
    H = the list of altitude [km]
    '''
    import numpy as np
    
    n_total = N
    n_closest = min(n_total,key=lambda x:abs(x-n_x))
    index_n = np.where(n_total == n_closest)[0][0]
    h_closest = H[index_n]
    if n_closest > n_x:
        n_1 = n_closest
        h_1 = h_closest
        n_2 = n_total[index_n+1]
        h_2 = H[index_n+1]
    elif n_closest < n_x:
        n_1 = n_total[index_n-1]
        h_1 = H[index_n-1]
        n_2 = n_closest
        h_2 = h_closest

    h_x = h_1 + (n_x-n_1)/(n_2-n_1)*(h_2-h_1)

    return h_x

def H_from_flux(Γ_x,Γ,H):
    '''
    Output: 
    h_x = altitude giving Γ_x [km]
    Input:
    Γ_x = the arbitrary flux [#/m^2/s]
    Γ = the list of flux [#/m^2/s]
    H = the list of altitude [km]
    '''
    import numpy as np
    
    Γ_total = Γ
    Γ_closest = min(Γ_total,key=lambda x:abs(x-Γ_x))
    index_Γ = np.where(Γ_total == Γ_closest)[0][0]
    h_closest = H[index_Γ]
    if Γ_closest > Γ_x:
        Γ_1 = Γ_closest
        h_1 = h_closest
        Γ_2 = Γ_total[index_Γ+1]
        h_2 = H[index_Γ+1]
    elif Γ_closest < Γ_x:
        Γ_1 = Γ_total[index_Γ-1]
        h_1 = H[index_Γ-1]
        Γ_2 = Γ_closest
        h_2 = h_closest

    h_x = h_1 + (Γ_x-Γ_1)/(Γ_2-Γ_1)*(h_2-h_1)
    return Γ_x

def v_th(T,m):
    '''
    Output: 
    v = thermal speed [m/s]
    Input:
    T = temperature [K]
    m = mass [kg]
    '''
    import numpy as np
    k_B = 1.380649e-23 # Boltzmann constant [J/K]
    v = np.sqrt(k_B*T/m)
    return v