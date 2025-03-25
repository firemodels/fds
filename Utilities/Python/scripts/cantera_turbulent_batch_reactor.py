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
def EDC_batch_reactor(mechanismFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName):
    print("Running case:"+caseName)
    current_working_directory = os.getcwd()
    Cantera_DIR = os.path.join(current_working_directory, "../Input_Libraries/Chemical_Mechanisms/Cantera/")
    Chemistry_DIR = os.path.join(current_working_directory, "../../Verification/Chemistry/")
    
    gas = ct.Solution(Cantera_DIR+mechanismFile)
    mixingstates = ct.SolutionArray(gas, extra=['t', 'internal_energy', 'enthalpy','C_mass', 'H_mass', 'N_mass', 'O_mass'])
    cellstates = ct.SolutionArray(gas, extra=['t', 'internal_energy', 'enthalpy', 'C_mass', 'H_mass', 'N_mass', 'O_mass'])
    
    elements = ['C', 'H', 'N', 'O']
    element_masses = {el: 0.0 for el in elements}
    
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
        zeta = zeta_0
        U = ct.Quantity(gas)
        M = ct.Quantity(gas,constant="HP")
        J = ct.Quantity(gas)
        U = zeta_0*C
        M = (1 - zeta_0)*C
    
        t_sub = 0.
        J.TPX = U.TPX
        while t_sub < dt_cfd:
            t_sub += dt_sub_chem
            zeta_last = zeta
            zeta = zeta_0*np.exp(-t_sub/tau_mix)
    
            J.mass = (zeta_last-zeta)*C.mass      # adjust mass of stream before adding
            M += J                                # mix U into M
            U.mass -= J.mass                      # reduce unmixed mass
            #input("Press Enter to continue...")
            # define reactor using only the mixed portion of the cell
            gas.TPX = M.T, M.P, M.X
            r = ct.IdealGasConstPressureReactor(gas)
            sim = ct.ReactorNet([r])
            sim.advance(dt_sub_chem)
            M.TPX = gas.TPX # now set the mixed zone based on reactor solution
            
            element_masses = {el: 0.0 for el in elements}
            for species in gas.species_names:
                species_index = gas.species_index(species)
                species_mass = M.Y[species_index] * C.mass  # Total mass of this species
                for el in elements:
                    element_masses[el] += species_mass * gas.n_atoms(species, el) * gas.atomic_weight(el)/ gas.molecular_weights[species_index]
               
            mixingstates.append(M.state, 
                      t=(t_cfd+t_sub)*1e3, 
                      internal_energy=M.int_energy_mass,
                      enthalpy=M.enthalpy_mass,
                      C_mass=element_masses['C'], 
                      H_mass=element_masses['H'], 
                      N_mass=element_masses['N'], 
                      O_mass=element_masses['O'])
    
        t_cfd += dt_cfd   
        
        # Solve for temperature while keeping internal energy fixed
        C = ct.Quantity(gas)
        newY=(U.mass*U.Y+M.mass*M.Y)/(U.mass+M.mass)   
        C.UVY=internal_energy_per_unit_mass, volume, newY
        
        element_masses = {el: 0.0 for el in elements}
        for species in gas.species_names:
            species_index = gas.species_index(species)
            species_mass = C.Y[species_index] * C.mass  # Total mass of this species
            for el in elements:
                element_masses[el] += species_mass * gas.n_atoms(species, el) * gas.atomic_weight(el)/ gas.molecular_weights[species_index]
    
        cellstates.append(C.state, 
                      t=t_cfd*1e3, 
                      internal_energy=M.int_energy_mass,
                      enthalpy=M.enthalpy_mass,
                      C_mass=element_masses['C'], 
                      H_mass=element_masses['H'], 
                      N_mass=element_masses['N'], 
                      O_mass=element_masses['O'])
    
    
    statesToSave=mixingstates
    
    columnNames = ['Time'] + ['OH'] + ['H']+['O2'] + [fuel]+ ['C2H2'] + ['C_mass']+['H_mass']+['N_mass']+['O_mass']+['Temperature'] + ['Pressure'] + ['Internal_energy']  + ['Enthalpy'] 
    statesForCSV = np.array([
        statesToSave.t, 
        statesToSave.Y[:, gas.species_index('OH')],
        statesToSave.Y[:, gas.species_index('H')],
        statesToSave.Y[:, gas.species_index('O2')],
        statesToSave.Y[:, gas.species_index(fuel)],
        statesToSave.Y[:, gas.species_index('C2H2')],
        statesToSave.C_mass,
        statesToSave.H_mass,
        statesToSave.N_mass,
        statesToSave.O_mass,
        statesToSave.T,
        statesToSave.P,
        statesToSave.internal_energy,
        statesToSave.enthalpy
    ]).T
    csvdata = pd.DataFrame(statesForCSV, columns=columnNames)
    outfilename = caseName+'_soln.csv'
    save_csv_file(Chemistry_DIR+outfilename,csvdata)

#------------------------------------------------------------------
#------------------------------------------------------------------
## Methane GRI cases with mixing substeps (only one cfd timestep).
#------------------------------------------------------------------
#------------------------------------------------------------------

mechFile='Methane_grimech30.yaml'
fuel = 'CH4'
equiv_ratio = 0.6
T0 = 1200  # Initial temperature [K]
P0 = ct.one_atm  # Initial pressure [Pa]
tau_mix = 0.01
t_end = 1e-1
dt_cfd = 1e-1
dt_sub_chem = 1e-4


# EDC_Methane_grimech30_T1200K_Phi0p6_Zeta1p0
caseName='EDC_Methane_grimech30_T1200K_Phi0p6_Zeta1p0'
zeta_0 = 1.0
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName)

# EDC_Methane_grimech30_T1200K_Phi0p6_Zeta1p0
caseName='EDC_Methane_grimech30_T1200K_Phi0p6_Zeta0p75'
zeta_0 = 0.75
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName)

# EDC_Methane_grimech30_T1200K_Phi0p6_Zeta1p0
caseName='EDC_Methane_grimech30_T1200K_Phi0p6_Zeta0p5'
zeta_0 = 0.5
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName)

# EDC_Methane_grimech30_T1200K_Phi0p6_Zeta1p0
caseName='EDC_Methane_grimech30_T1200K_Phi0p6_Zeta0p25'
zeta_0 = 0.25
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName)

# EDC_Methane_grimech30_T1200K_Phi0p6_Zeta1p0
caseName='EDC_Methane_grimech30_T1200K_Phi0p6_Zeta0p0'
zeta_0 = 0.0
EDC_batch_reactor(mechFile, T0, P0, fuel, equiv_ratio, zeta_0, tau_mix, t_end, dt_cfd, dt_sub_chem, caseName)


#------------------------------------------------------------------
#------------------------------------------------------------------
## Ignition delay Methane GRI cases (many cfd timestep).
#------------------------------------------------------------------
#------------------------------------------------------------------
