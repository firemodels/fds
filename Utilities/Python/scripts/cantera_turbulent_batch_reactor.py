"""
McDermott
21 May 2024

Turbulent Batch Reactor (TBR) model
Constant-pressure, adiabatic kinetics simulation.
"""

import sys
import os
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Helper function to correct mass fractions
def save_csv_file(fullFilePath, csvdata):
    csvdata = csvdata.apply(lambda x: x.str.strip() if x.dtype == 'object' else x)
    csvdata.replace('nan','',inplace=True)
    csvdata.replace('NAN','',inplace=True)
    csvdata = csvdata.map(lambda x: f'{x:10.5f}' if isinstance(x, (float)) else x)
    csvdata.to_csv(fullFilePath,index=False)


# EDC turbulent batch reactor
def EDC_batch_reactor(mechanismFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName, saveStateOpt, writeInterval):
    print("Running case:"+caseName)
    current_working_directory = os.getcwd()
    Cantera_DIR = os.path.join(current_working_directory, "../../Input_Libraries/Chemical_Mechanisms/Cantera/")
    Chemistry_DIR = os.path.join(current_working_directory, "../../../Verification/Chemistry/")
    
    gas = ct.Solution(Cantera_DIR+mechanismFile)
    mixingstates = ct.SolutionArray(gas, extra=['t', 'internal_energy', 'enthalpy','C_MF', 'H_MF', 'N_MF', 'O_MF'])
    cellstates = ct.SolutionArray(gas, extra=['t', 'internal_energy', 'enthalpy', 'C_MF', 'H_MF', 'N_MF', 'O_MF'])
    
    elements = ['C', 'H', 'N', 'O']
    element_MFs = {el: 0.0 for el in elements}
    
    # Initial conditions
    equiv_ratio = 0.6
    init_temp = T0
    init_pr=P0
        
    # initial state of the cell
    C = ct.Quantity(gas)
    C.TP = init_temp, init_pr
    C.set_equivalence_ratio(equiv_ratio, {fuel:1}, {'O2':1,'N2':3.76})
    C.mass = 1.0
    internal_energy_per_unit_mass = C.int_energy_mass
    volume = C.v
    

    
    # combustion substep time integration
    t_cfd = 0.
    while t_cfd < t_end:
        
        U = ct.Quantity(gas)
        M = ct.Quantity(gas,constant="HP")
        J = ct.Quantity(gas)
        U = zeta_0*C
        M = (1 - zeta_0)*C
        zeta = zeta_0
    
        t_sub = 0.
        J.TPX = U.TPX
        while t_sub < dt_cfd:
            t_sub += dt_sub_chem
            zeta_last = zeta
            zeta = zeta_0*np.exp(-t_sub/tau_mix)
    
            J.mass = (zeta_last-zeta)*C.mass      # adjust mass of stream before adding
            M += J                                # mix U into M
            U.mass -= J.mass                      # reduce unmixed mass
            
            # define reactor using only the mixed portion of the cell
            gas.TPX = M.T, M.P, M.X
            r = ct.IdealGasConstPressureReactor(gas)
            sim = ct.ReactorNet([r])
            sim.advance(dt_sub_chem)
            M.TPX = gas.TPX # now set the mixed zone based on reactor solution
            
            element_MFs = {el: 0.0 for el in elements}
            for species in gas.species_names:
                species_index = gas.species_index(species)
                species_mass = M.Y[species_index] # Total mass of this species
                for el in elements:
                    element_MFs[el] += species_mass * gas.n_atoms(species, el) * gas.atomic_weight(el)/ gas.molecular_weights[species_index]
               
            mixingstates.append(M.state, 
                      t=(t_cfd+t_sub), 
                      internal_energy=M.int_energy_mass,
                      enthalpy=M.enthalpy_mass,
                      C_MF=element_MFs['C'], 
                      H_MF=element_MFs['H'], 
                      N_MF=element_MFs['N'], 
                      O_MF=element_MFs['O'])
    
        t_cfd += dt_cfd   
        
        # Solve for temperature while keeping internal energy fixed
        C = ct.Quantity(gas)
        newY=(U.mass*U.Y+M.mass*M.Y)/(U.mass+M.mass)   
        C.UVY=internal_energy_per_unit_mass, volume, newY
        
        element_MFs = {el: 0.0 for el in elements}
        for species in gas.species_names:
            species_index = gas.species_index(species)
            species_mass = C.Y[species_index] # Total mass of this species
            for el in elements:
                element_MFs[el] += species_mass * gas.n_atoms(species, el) * gas.atomic_weight(el)/ gas.molecular_weights[species_index]
    
        cellstates.append(C.state, 
                      t=t_cfd, 
                      internal_energy=C.int_energy_mass,
                      enthalpy=C.enthalpy_mass,
                      C_MF=element_MFs['C'], 
                      H_MF=element_MFs['H'], 
                      N_MF=element_MFs['N'], 
                      O_MF=element_MFs['O'])
    
    
        if saveStateOpt == 1:
            statesToSave = mixingstates
        else:
            statesToSave = cellstates
    
    columnNames = ['Time']+['TMP']+['Pressure']+['Enthalpy']+['OH']+['O2']+[fuel]+['C_MF']+['H_MF']+['N_MF']+['O_MF']
    statesForCSV = np.array([
        statesToSave.t, 
        statesToSave.T-273.15,
        statesToSave.P,
        statesToSave.enthalpy,
        statesToSave.Y[:, gas.species_index('OH')],
        statesToSave.Y[:, gas.species_index('O2')],
        statesToSave.Y[:, gas.species_index(fuel)],
        statesToSave.C_MF,
        statesToSave.H_MF,
        statesToSave.N_MF,
        statesToSave.O_MF,
    ]).T
    
    # Reduce the number of rows of statesForCSV
    statesForCSVReduced = []
    tCount = 0
    for j in range(len(statesForCSV)):
        if (tCount%writeInterval ==0):
            statesForCSVReduced.append(statesForCSV[j])
        tCount = tCount+1
    csvdata = pd.DataFrame(statesForCSVReduced, columns=columnNames)
    outfilename = caseName+'_soln.csv'
    save_csv_file(Chemistry_DIR+outfilename,csvdata)

#------------------------------------------------------------------
#------------------------------------------------------------------
## Vary Zeta0: Methane GRI cases with mixing substeps (only one cfd timestep).
#------------------------------------------------------------------
#------------------------------------------------------------------


# Helper function to get default CFD step parameters
def get_OneCFDStep_vary_zeta0_default_parameters():
    equiv_ratio = 0.6
    T0 = 1200  # Initial temperature [K]
    P0 = ct.one_atm  # Initial pressure [Pa]
    tau_mix = 0.01
    t_end = 1e-1
    dt_cfd = 1e-1
    dt_sub_chem = 1e-4
    writeInterval=10
    return equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval

# Define mechanism and fuel
mechFile = 'Methane_grimech30.yaml'
fuel = 'CH4'

caseName='EDC_OneCFDStep_Methane_grimech30_Zeta1p0'
equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval = get_OneCFDStep_vary_zeta0_default_parameters()
zeta_0 = 1.0
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,1,writeInterval)

caseName='EDC_OneCFDStep_Methane_grimech30_Zeta0p75'
equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval = get_OneCFDStep_vary_zeta0_default_parameters()
zeta_0 = 0.75
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,1,writeInterval)

caseName='EDC_OneCFDStep_Methane_grimech30_Zeta0p5'
equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval = get_OneCFDStep_vary_zeta0_default_parameters()
zeta_0 = 0.5
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,1,writeInterval)

caseName='EDC_OneCFDStep_Methane_grimech30_Zeta0p25'
equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval = get_OneCFDStep_vary_zeta0_default_parameters()
zeta_0 = 0.25
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,1,writeInterval)

caseName='EDC_OneCFDStep_Methane_grimech30_Zeta0p0'
equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval = get_OneCFDStep_vary_zeta0_default_parameters()
zeta_0 = 0.0
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,1,writeInterval)


# #------------------------------------------------------------------
# #------------------------------------------------------------------
# ## Vary Zeta0: Ignition delay Methane GRI cases (multiple cfd timestep).
# #------------------------------------------------------------------
# #------------------------------------------------------------------

# Helper function to get default CFD step parameters
def get_MultiCFDStep_vary_zeta0_default_parameters():
    equiv_ratio = 0.6
    T0 = 900  # Initial temperature [K]
    P0 = ct.one_atm  # Initial pressure [Pa]
    tau_mix = 0.01
    t_end = 20
    dt_cfd = 1e-2
    dt_sub_chem = 1e-4
    writeInterval=10
    return equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval

mechFile='Methane_grimech30.yaml'
fuel = 'CH4'

caseName='EDC_MultiCFDStep_Methane_grimech30_Zeta1p0'
equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval = get_MultiCFDStep_vary_zeta0_default_parameters()
zeta_0 = 1.0
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,2,writeInterval)

caseName='EDC_MultiCFDStep_Methane_grimech30_Zeta0p75'
equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval = get_MultiCFDStep_vary_zeta0_default_parameters()
zeta_0 = 0.75
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,2,writeInterval)

caseName='EDC_MultiCFDStep_Methane_grimech30_Zeta0p5'
equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval = get_MultiCFDStep_vary_zeta0_default_parameters()
zeta_0 = 0.5
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,2,writeInterval)

caseName='EDC_MultiCFDStep_Methane_grimech30_Zeta0p25'
equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval = get_MultiCFDStep_vary_zeta0_default_parameters()
zeta_0 = 0.25
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,2,writeInterval)

caseName='EDC_MultiCFDStep_Methane_grimech30_Zeta0p0'
equiv_ratio, T0, P0, tau_mix, t_end, dt_cfd, dt_sub_chem,writeInterval = get_MultiCFDStep_vary_zeta0_default_parameters()
zeta_0 = 0.0
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,2,writeInterval)


#------------------------------------------------------------------
#------------------------------------------------------------------
## Vary tau_mix: Methane GRI cases with mixing substeps (only one cfd timestep).
#------------------------------------------------------------------
#------------------------------------------------------------------


# Helper function to get default CFD step parameters
def get_OneCFDStep_vary_taumix_default_parameters():
    equiv_ratio = 0.6
    T0 = 1200  # Initial temperature [K]
    P0 = ct.one_atm  # Initial pressure [Pa]
    zeta_0 = 1.0
    t_end = 1e-1
    dt_cfd = 1e-1
    dt_sub_chem = 1e-4
    writeInterval=10
    return equiv_ratio, T0, P0, zeta_0, t_end, dt_cfd, dt_sub_chem,writeInterval

# Define mechanism and fuel
mechFile = 'Methane_grimech30.yaml'
fuel = 'CH4'

caseName='EDC_OneCFDStep_Methane_grimech30_taumix0p1'
equiv_ratio, T0, P0, zeta_0, t_end, dt_cfd, dt_sub_chem,writeInterval = get_OneCFDStep_vary_taumix_default_parameters()
tau_mix = 0.1
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,1,writeInterval)

caseName='EDC_OneCFDStep_Methane_grimech30_taumix0p01'
equiv_ratio, T0, P0, zeta_0, t_end, dt_cfd, dt_sub_chem,writeInterval = get_OneCFDStep_vary_taumix_default_parameters()
tau_mix = 0.01
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,1,writeInterval)

caseName='EDC_OneCFDStep_Methane_grimech30_taumix0p001'
equiv_ratio, T0, P0, zeta_0, t_end, dt_cfd, dt_sub_chem,writeInterval = get_OneCFDStep_vary_taumix_default_parameters()
tau_mix = 0.001
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,1,writeInterval)

caseName='EDC_OneCFDStep_Methane_grimech30_taumix0p0001'
equiv_ratio, T0, P0, zeta_0, t_end, dt_cfd, dt_sub_chem,writeInterval = get_OneCFDStep_vary_taumix_default_parameters()
tau_mix = 0.0001
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,1,writeInterval)

caseName='EDC_OneCFDStep_Methane_grimech30_taumix0p00001'
equiv_ratio, T0, P0, zeta_0, t_end, dt_cfd, dt_sub_chem,writeInterval = get_OneCFDStep_vary_taumix_default_parameters()
tau_mix = 0.00001
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName,1,writeInterval)


