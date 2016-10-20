import numpy as np
from scipy.integrate import quad

# Define some universal constants
ry2ev = 13.605698066
#kb = 8.6173324E-5      #eV/K
kb = 8.6173324E-5/ry2ev #Ry/K

def interpolate_y(x,y,x0):
    """Linear interpolation to find arbitrary y0 for a given x0 in x vs. y data."""
    #
    # First we must between which sws entries sws0 lies:
    index0 = 0
    for i in range(len(x[:])-1):
        xMin = x[i]
        xMax = x[i+1]
        yMin = y[i]
        yMax = y[i+1]
        if x0 > xMin and x0 < xMax:
            break
    # Construct line connecting (xMin,yMin) and (xMax,yMax)
    A = (yMin-yMax)/(xMin-xMax)
    B = (xMin*yMax-xMax*yMin)/(xMin-xMax)
    y0 = A*x0 + B
    return y0

def el_ent_parser(filename):
    """Extract electronic entropies from the KGRN input file.
    
    Returns the entropy and the temperature at which it was calculated.
    Units are eV and Kelvin.
    """
    file = open(filename,'r')
    lines = file.readlines()
    nxtsws_tag = " NXTSWS:   IQ IT ITA MMT Type  NZ ION  ELN   QTR   SPLIT  FIX  CONC"
    alat_tag = "Alat ="
    entrop_tag = " ENTROP:"
    #nxtsws_tag = "NXTSWS:"
    for i in range(len(lines)):
        if nxtsws_tag in lines[i]:
            indMin = i+2
        elif alat_tag in lines[i]:
            indMax = i-2
            break
    #
    concs = np.zeros(indMax + 1 - indMin)
    entrops = np.zeros(indMax + 1 - indMin)
    its = np.zeros(indMax + 1 - indMin)
    #print('num_sites = ',num_sites)
    ind_tmp = 0
    for i in range(indMin,indMax+1):
        concs[ind_tmp] = float(lines[i].split()[-1])
        its[ind_tmp] = int(lines[i].split()[1])
        #print(concs[ind_tmp])
        #it,ita = int(lines[i].split()[1]),int(lines[i].split()[2])
        ind_tmp += 1
    #
    num_sites = np.max(its)
    #
    ind_tmp = 0
    for i in range(indMax,len(lines)):
        if entrop_tag in lines[i]:
            entrops[ind_tmp] = float(lines[i].split()[-4])
            #print(entrops[ind_tmp])
            ind_tmp += 1
            if ind_tmp == len(entrops):
                ind_tfermi = i+2
                break
    #
    tfermi = float(lines[ind_tfermi].split()[-2])
    #
    ent_tot = 0.0
    for i in range(len(concs)):
        ent_tot += concs[i]*entrops[i]/tfermi*ry2ev
    ent_tot /= num_sites
    return ent_tot,tfermi

def dos_parser(filename):
    """Extract DOS at Fermi level from the KGRN input file.
    
    Returns the entropy and the temperature at which it was calculated.
    Units are eV and Kelvin.
    """
    file = open(filename,'r')
    lines = file.readlines()
    nxtsws_tag = " NXTSWS:   IQ IT ITA MMT Type  NZ ION  ELN   QTR   SPLIT  FIX  CONC"
    alat_tag = "Alat ="
    dos_tag = " Dos(Ef)    ="
    mag_tag = " Magn. mom. ="
    hop_tag = " Hopfield   ="
    for i in range(len(lines)):
        if nxtsws_tag in lines[i]:
            indMin = i+2
        elif alat_tag in lines[i]:
            indMax = i-2
            break
    #
    concs = np.zeros(indMax + 1 - indMin)
    doses = np.zeros(indMax + 1 - indMin)
    its = np.zeros(indMax + 1 - indMin)
    #
    ind_tmp = 0
    for i in range(indMin,indMax+1):
        concs[ind_tmp] = float(lines[i].split()[-1])
        its[ind_tmp] = int(lines[i].split()[1])
        ind_tmp += 1
    #
    print("concs = ",concs)
    print("its = ",its)
    #
    num_sites = np.max(its)
    ind_tmp = len(doses) - 1
    # Because KGRN output file is different for non- and magnetic calculations,
    # we have to do some additional checks to make sure we are reading the right
    # values.
    for i in range(len(lines)-1,indMax,-1):
        if dos_tag in lines[i]:
            #print(lines[i])
            if mag_tag in lines[i+1] or hop_tag in lines[i+1]:
                #print(lines[i+1])
                doses[ind_tmp] = float(lines[i].split()[-1])
                ind_tmp -= 1
            if ind_tmp == -1:
                break
    #
    for i in range(len(doses)):
        print(doses[i])
    #
    dos_tot = 0.0
    for i in range(len(concs)):
        dos_tot += concs[i]*doses[i]
    dos_tot /= num_sites
    print('dos_tot = ',dos_tot)
    dos_tot /= ry2ev
    return dos_tot  

def S_conf(concs):
    """
    
    """
    summa = 0.0
    for i in range(len(concs)):
        summa += concs[i]*np.log(concs[i])
    summa *= -kb
    return summa

def F_conf(S_conf,T):
    """
    S_conf = configurational entropy in units of Ry/Kelvin,
    T      = temperature in Kelvin.
    """
    return -T*S_conf

def S_elec(dos,T):
    """Electronic entropy. Linear approximation is used, where S(T=0K) = 0
    and S(T=T0) = S0.
    
    dos should be given in units of 1/Ry.
    """
    return np.pi**2/3*kb**2*T*dos

def F_elec(dos,T):
    #print('T, Selec = ',T,S_elec(dos,T)/kb)
    return -T*S_elec(dos,T)

def S_vib(r,r0,B0,M,grun,T):
    """Vibrational entropy."""
    x = debye_temp(r,r0,B0,M,grun) / T
    return 3*kb*(4.0/3*debye_func(x) - np.log(1-np.exp(-x)))
    
def debye_integrand(x):
    return x**3/(np.exp(x) - 1)

def debye_func(x):
    return 3.0/x**3*quad(debye_integrand,0,x)[0]

def debye_temp(r,r0,B0,M,grun):
    """
    r    = WS-radius should be given in bohr,
    r0   = Eq. WS-radius should be given in bohr,
    B0     should be given in GPa, 1GPa = 10 kbar,
    M      should be given in a.m.u's,
    grun = Gruneisen parameter, unitless.
    """
    B0 = 10*B0 # Convert input B0 in GPa to kbar
    C = 41.63 # In units of K*s
    V0 = 4.0*np.pi*r0**3/3
    V  = 4.0*np.pi*r**3/3
    DT0 = C*np.sqrt(r0*B0/M)
    return DT0 * (V0/V)**grun

def E_D(r,r0,B0,M,grun,T):
    x = debye_temp(r,r0,B0,M,grun) / T
    return 9.0/8*kb*debye_temp(r,r0,B0,M,grun) + 3*kb*T*debye_func(x)

def F_vib(r,r0,B0,M,grun,T):
    """
    r    = WS-radius should be given in bohr,
    r0   = Eq. WS-radius should be given in bohr,
    B0     should be given in GPa, 1GPa = 10 kbar,
    M      should be given in a.m.u's,
    grun = Gruneisen parameter,
    T    = temperature in Kelvin.
    """
    return - T*S_vib(r,r0,B0,M,grun,T) + E_D(r,r0,B0,M,grun,T)

def S_mag(concs,moms):
    """
    r     = WS-radius should be given in bohr,
    concs = list of concentrations,
    moms  = list of mag. moms.
    """
    summa = 0.0
    for i in range(lenconcs):
        summa += concs[i]*np.log(np.abs(moms[i])+1)
    summa /= n_sites
    summa *= kb
    return summa
            
#def F_mag(T,concs,moms):
#    """
#    r    = WS-radius should be given in bohr,
#    concs = list of concentrations,
#    moms  = list of mag. moms.
#    T    = temperature in Kelvin.
#    """
#    return -T*S_mag(r,concs,moms)

def F_mag(S_mag,T):
    """
    S_mag = magnetic entropy in units of Ry/Kelvin,
    T     = temperature in Kelvin.
    """
    return -T*S_mag

def T_Curie(E_DLM,E_FM):
    """
    """
    return 2.0/3*(E_DLM-E_FM)/kb

#def F_gibbs(EDFT,):
#    """
#    """
#    return EDFT + F_elec(dos,T) + F_vib(r,r0,B0,M,grun,T) + F_mag(T,concs,moms)
