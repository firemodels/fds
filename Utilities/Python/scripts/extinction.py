"""
extinction.py
Converted from MATLAB script extinction.m to Python.
-----------------------
C Weinschenk
Verification of Extinction Model
Fuel: Methane
6/2012
-----------------------
"""
import os
import numpy as np
import pandas as pd
import fdsplotlib

#-----------------------
# Initialize species
#-----------------------

#-------------
# vector key
# (1) = nitrogen
# (2) = methane
# (3) = oxygen
# (4) = carbon dioxide
# (5) = water vapor
#-------------
T_crit = 1780          # critical flame temperature for methane, Beyler T(OI) [K]
Y_O2_0 = 0.156         # Beyler OI converted to a mass fraction
mass = 0.1             # mass [kg]
volume = 0.001         # volume [m^3]
R = 8.3145             # gas constant [J/mol K]

# Parameters
nu = [0, -1, -2, 1, 2] # stoichiometric coefficients
y_MW = [28.0134, 16.042460, 31.9988, 44.0095, 18.01528] # [g/mol]
y_hf = [0.0, -74873, 0.0, -393522, -241826] # [J/mol]
pres_0 = 1.013253421185575e+05 # initial presure [Pa]

# -----------------------------
# Paths and identifiers
# -----------------------------
firemodels_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),'..','..','..','..'))
fds_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),'..','..','..'))
#outdir = os.path.join(firemodels_dir, 'out','Convection','')
results_dir = os.path.join(fds_dir, 'Verification','Extinction','')
pltdir = os.path.join(fds_dir, 'Manuals','FDS_User_Guide','SCRIPT_FIGURES','')
os.makedirs(pltdir, exist_ok=True)

# Initialize variable arrays
num_samp = 100         # number of input sets
ignite = np.zeros((num_samp,2))
extinct_o2 = np.zeros((num_samp,2))
fds_ext_o2 = np.zeros((num_samp,2))
extinct_fuel = np.zeros((num_samp,2))
fuel_seed = np.zeros((num_samp,))
temp_seed = np.zeros((num_samp,))
phi_fuel = np.zeros((num_samp,))
phi = np.zeros((num_samp,5))
T0 = np.zeros((num_samp,))
phi_st = np.zeros((3,))
phi_r = np.zeros((3,))
RR = np.zeros((3,))
d_phi = np.zeros((num_samp,5))
phi_new = np.zeros((num_samp,5))
coeff_ref = np.zeros((7,5))
coeff0 = np.zeros((7,5))
coeff = np.zeros((7,5))
del_h_init = np.zeros((5,))
del_h = np.zeros((5,))
Q = np.zeros((num_samp,))
hrr_ext = np.zeros((26, num_samp))
temp_ext = np.zeros((26, num_samp))
o2_ext = np.zeros((26, num_samp))
fu_ext = np.zeros((26, num_samp))
fds_ignite = np.zeros((num_samp,2))

mult = 1.
seed = np.round(np.linspace(0.0029,0.0546, 10), decimals=4)  # array of fuel mass fractions
seed = seed/mult
for ii in range(0, 10):
    for jj in range(0, 10):
        fuel_seed[jj+10*(ii-1)] = seed[ii]

seed2 = np.linspace(300,1875,10)  # array of temperatures (K)
for ii in range(0, 10):
    for jj in range(0, 10):
        temp_seed[jj+10*(ii-1)] = seed2[jj]

o2_mult = 4.*mult

for i in range(0, num_samp):
    phi_fuel[i] = fuel_seed[i] # initial mass fraction of fuel
    phi[i,:] = [(1-(o2_mult+1)*fuel_seed[i]), fuel_seed[i], o2_mult*fuel_seed[i], 0.0, 0.0] # initial mass fractions of all species
    T0[i] = temp_seed[i] # initial temperature [K]


#-----------------------
# Determine max change in species
#-----------------------
for jj in range(0,num_samp):
    nu_mw_sum = 0
    phi_sum = 0

    for i in range(0, len(nu)):
        if nu[i] < 0:
            nu_mw_sum = nu_mw_sum + abs(nu[i])*y_MW[i]
            phi_sum = phi_sum + phi[jj,i]

    for j in range(0, len(nu)):
        if nu[i] < 0:
            phi_st[i] = abs(nu[i])*y_MW[i]/nu_mw_sum
            phi_r[i] = phi[jj,i]/phi_sum
    for i in range(0, len(RR)):
        if phi_st[i] > 0:
            RR[i] == phi_r[i]/phi_st[i]
            
    ILR = 2
    for i in range(2, len(nu)):
        if nu[i] < 0:
            if RR[i] < RR[ILR]:
                ILR = i

    for i in range(0, len(nu)):
        d_phi[jj,i] = nu[i]*y_MW[i]/(abs(nu[ILR])*y_MW[ILR])*phi[jj,ILR]
        
    phi_new[jj,:] = phi[jj,:] + d_phi[jj,:]

    #-----------------------
    # Determine Extinction
    #-----------------------
    #---------
    #NOTE: coefficients were found from NIST Webbook
    #---------
    #reference temperature coeffs
    coeff_ref[:,0] = [28.98641,1.853978,-9.647459,16.63537,0.000117,-8.671914,0.0]
    coeff_ref[:,1] = [-0.703029,108.4773,-42.52157,5.862788,0.678565,-76.84376,-74.87310]
    coeff_ref[:,2] = [31.32234,-20.23531,57.86644,-36.50624,-0.007374,-8.903471,0.0]
    coeff_ref[:,3] = [24.99735,55.18696,-33.69137,7.948387,-0.136638,-403.6075,-393.5224]
    coeff_ref[:,4] = [30.09200,6.832514,6.793435,-2.534480,0.082139,-250.8810,-241.8264]

    #initial temperature coeffs
    #nitrogen cp coeffs [J/mol K]
    if T0[jj] <=500:
        coeff0[:,0] = [28.98641,1.853978,-9.647459,16.63537,0.000117,-8.671914,0.0]
    elif 500 < T0[jj] <= 2000:
        coeff0[:,0] = [19.50583,19.88705,-8.598535,1.369784,0.527601,-4.935202,0.0]
    else:
        coeff0[:,0] = [35.51872,1.128728,-0.196103,0.014662,-4.553760,-18.97091,0.0]
    
    #methane cp coeffs [J/mol K]
    if T0[jj] <= 1300:
        coeff0[:,1] = [-0.703029,108.4773,-42.52157,5.862788,0.678565,-76.84376,-74.87310]
    else:
        coeff0[:,1] = [85.81217,11.26467,-2.114146,0.138190,-26.42221,-153.5327,-74.87310]
    
    #oxygen cp coeffs [J/mol K]
    if T0[jj] <=700:
        coeff0[:,2] = [31.32234,-20.23531,57.86644,-36.50624,-0.007374,-8.903471,0.0]
    elif (T0[jj] <= 2000):
        coeff0[:,2] = [30.03235,8.772972,-3.988133,0.788313,-0.741599,-11.32468,0.0]
    else:
        coeff0[:,2] = [20.91111,10.72071,-2.020498,0.146449,9.245722,5.337651,0.0]
    
    #carbon dooxide cp coeffs [J/mol K]
    if T0[jj] <= 1200:
        coeff0[:,3] = [24.99735,55.18696,-33.69137,7.948387,-0.136638,-403.6075,-393.5224]
    else:
        coeff0[:,3] = [58.16639,2.720074,-0.492289,0.038844,-6.447293,-425.9186,-393.5224]
    
    #water vapor cp coeffs [J/mol K]
    if T0[jj] <= 1700:
        coeff0[:,4] = [30.09200,6.832514,6.793435,-2.534480,0.082139,-250.8810,-241.8264]
    else:
        coeff0[:,4] = [41.96246,8.622053,-1.499780,0.098199,-11.15764,-272.1797,-241.8264]
    

    #critical temperature coeffs
    #nitrogen cp coeffs [J/mol K]
    if T_crit <=500:
        coeff[:,0] = [28.98641,1.853978,-9.647459,16.63537,0.000117,-8.671914,0.0]
    elif T_crit <= 2000:
        coeff[:,0] = [19.50583,19.88705,-8.598535,1.369784,0.527601,-4.935202,0.0]
    else:
        coeff[:,0] = [35.51872,1.128728,-0.196103,0.014662,-4.553760,-18.97091,0.0]
    
    #methane cp coeffs [J/mol K]
    if T_crit <= 1300:
        coeff[:,1] = [-0.703029,108.4773,-42.52157,5.862788,0.678565,-76.84376,-74.87310]
    else:
        coeff[:,1] = [85.81217,11.26467,-2.114146,0.138190,-26.42221,-153.5327,-74.87310]
    
    #oxygen cp coeffs [J/mol K]
    if T_crit <=700:
        coeff[:,2] = [31.32234,-20.23531,57.86644,-36.50624,-0.007374,-8.903471,0.0]
    elif T_crit <= 2000:
        coeff[:,2] = [30.03235,8.772972,-3.988133,0.788313,-0.741599,-11.32468,0.0]
    else:
        coeff[:,2] = [20.91111,10.72071,-2.020498,0.146449,9.245722,5.337651,0.0]
    
    #carbon dooxide cp coeffs [J/mol K]
    if T_crit <= 1200:
        coeff[:,3] = [24.99735,55.18696,-33.69137,7.948387,-0.136638,-403.6075,-393.5224]
    else:
        coeff[:,3] = [58.16639,2.720074,-0.492289,0.038844,-6.447293,-425.9186,-393.5224]
    
    #water vapor cp coeffs [J/mol K]
    if T_crit <= 1700:
        coeff[:,4] = [30.09200,6.832514,6.793435,-2.534480,0.082139,-250.8810,-241.8264]
    else:
        coeff[:,4] = [41.96246,8.622053,-1.499780,0.098199,-11.15764,-272.1797,-241.8264]
    
    t=T_crit/1000
    t0=T0[jj]/1000
    for i in range(0, 5):
        del_h_init[i] = 1000*(coeff0[0,i]*t0 + (1/2)*coeff0[1,i]*(t0**2) + (1/3)*coeff0[2,i]*(t0**3) + (1/4)*coeff0[3,i]*(t0**4) - coeff0[4,i]/t0 + coeff0[5,i] - coeff0[6,i])
        del_h[i] = 1000*(coeff[0,i]*t + (1/2)*coeff[1,i]*(t**2) + (1/3)*coeff[2,i]*(t**3) + (1/4)*coeff[3,i]*(t**4) - coeff[4,i]/t + coeff[5,i] - coeff[6,i])
    

    Q[jj] = 0
    h = 0
    h0 = 0
    for i in range(0, len(y_hf)):
        Q[jj] = Q[jj] - y_hf[i]*d_phi[jj,i]/y_MW[i]
        h0 = h0 +phi[jj,i]*del_h_init[i]/y_MW[i]
        h = h + phi[jj,i]*del_h[i]/y_MW[i]

    if h0 + Q[jj]  > h:
        ignite[jj,:] = [T0[jj],phi[jj,2]]
    else:
        extinct_o2[jj,:] = [T0[jj],phi[jj,2]]

#-----------------------
# Import FDS results
#-----------------------
epsilon = 1e-10

fds_file = [os.path.join(results_dir, 'extinction_1_devc.csv'),
            os.path.join(results_dir + 'extinction_2_devc.csv')]

X_ignite=np.zeros((num_samp,2))
X_fds_ignite=np.zeros((num_samp,2))
X_extinct_o2=np.zeros((num_samp,2))
X_fds_ext_o2=np.zeros((num_samp,2))

inds = [1, 0]
for ifile in inds:
    
    if os.path.exists(fds_file[ifile]) is False:
        print('Error: File %s does not exist. Skipping case.'%(fds_file[ifile]))
    
    extinct = pd.read_csv(fds_file[ifile], header=1)
    extinct_1=extinct.values # data
    
    for i in range(0, 100):
       hrr_ext[:,i] = extinct_1[:,4*(i+1)-3]
       temp_ext[:,i] = extinct_1[:,4*(i+1)-2]+273.15
       o2_ext[:,i] = extinct_1[:,4*(i+1)-1]
       fu_ext[:,i] = extinct_1[:,4*(i+1)]
    
    for i in range(0, 100):
        if np.sum(hrr_ext[:,i]) > epsilon:
            fds_ignite[i,:] = [temp_ext[0,i],o2_ext[0,i]]
            fds_ext_o2[i,:] = [0,0]
        else:
            fds_ext_o2[i,:] = [temp_ext[0,i],o2_ext[0,i]]
            fds_ignite[i,:] = [0,0]
            
    # Simple Extinction Model
    
    if ifile==0:
       simple_o2 = np.array([Y_O2_0,0.0939,0])
       simple_temp = np.array([293.15,873,873])
    else:
       simple_o2 = np.array([Y_O2_0,0])
       simple_temp = np.array([293.15,T_crit])

    
    # Make the plot
    tmpm = 273
    X_simple_o2 = (simple_o2/32)/(simple_o2/32 + (1-simple_o2)/28)
    X_ignite[:,1] = (ignite[:,1]/32)/(ignite[:,1]/32 + (1-ignite[:,1])/28)
    X_fds_ignite[:,1] = (fds_ignite[:,1]/32)/(fds_ignite[:,1]/32 + (1-fds_ignite[:,1])/28)
    X_extinct_o2[:,1] = (extinct_o2[:,1]/32)/(extinct_o2[:,1]/32 + (1-extinct_o2[:,1])/28)
    X_fds_ext_o2[:,1] = (fds_ext_o2[:,1]/32)/(fds_ext_o2[:,1]/32 + (1-fds_ext_o2[:,1])/28)
    
    
    if ifile==0: # Modify expected behavior
       for jj in range(0, num_samp):
          if ((extinct_o2[jj,0]-tmpm>600) and (X_extinct_o2[jj,1]>0)):
             X_ignite[jj,1]=X_extinct_o2[jj,1]
             X_extinct_o2[jj,1]=0
             ignite[jj,0]=extinct_o2[jj,0]
             extinct_o2[jj,0]=0
    
    git_file = os.path.join(results_dir, 'extinction_1_git.txt')
    version_string = fdsplotlib.get_version_string(git_file)
    fig = fdsplotlib.plot_to_fig(x_data=simple_temp-tmpm, y_data=X_simple_o2, marker_style='k-', data_label='Simple Model', linewidth=1,
                                 x_min=0.0, x_max=1700, y_min=0, y_max=0.21, xticks=[0, 500, 1000, 1500], yticks=np.linspace(0, 0.2, 11),
                                 revision_label=version_string, figure_handle=None, 
                                 x_label='Temperature (°C)',
                                 y_label='Oxygen Volume Fraction',
                                 figure_size=(6,6), plot_size=(4.5,4.5), plot_origin=(1.125,1.0))
    fdsplotlib.plot_to_fig(x_data=ignite[:,0]-tmpm, y_data=X_ignite[:,1], marker_style='rs', data_label='Expected Burning', 
                           marker_fill_color=(1,1,1,0.0), markersize=6, markeredgewidth=0.5, figure_handle=fig)
    fdsplotlib.plot_to_fig(x_data=fds_ignite[:,0]-tmpm, y_data=X_fds_ignite[:,1], marker_style='r+', data_label='FDS Burning', 
                           marker_fill_color=(1,1,1,0.0), markersize=4, markeredgewidth=0.5, figure_handle=fig)
    fdsplotlib.plot_to_fig(x_data=extinct_o2[:,0]-tmpm, y_data=X_extinct_o2[:,1], marker_style='bo', data_label='Expected Extinction', 
                           marker_fill_color=(1,1,1,0.0), markersize=7, markeredgewidth=0.5, figure_handle=fig)
    fdsplotlib.plot_to_fig(x_data=fds_ext_o2[:,0]-tmpm, y_data=X_fds_ext_o2[:,1], marker_style='b*', data_label='FDS Extinction', 
                           markersize=4, markeredgewidth=0.5, figure_handle=fig)
    fig.axes[0].xaxis.set_tick_params(pad=10)
    fig.axes[0].yaxis.set_tick_params(pad=10)
    local_pdf = os.path.join(pltdir, 'extinction_%d.pdf'%(ifile+1))
    fig.savefig(local_pdf)
    
    
