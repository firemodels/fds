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



current_working_directory = os.getcwd()
Cantera_DIR = os.path.join(current_working_directory, "../../Input_Libraries/Chemical_Mechanisms/Cantera/")
Chemistry_DIR = os.path.join(current_working_directory, "../../../Verification/Chemistry/")

gas = ct.Solution(Cantera_DIR+"Methane_grimech30.yaml")
mixingstates = ct.SolutionArray(gas, extra=['t'])

cellgas = ct.Solution(Cantera_DIR+"Methane_grimech30.yaml")
cellstates = ct.SolutionArray(cellgas, extra=['t', 'C_mass', 'H_mass', 'N_mass', 'O_mass'])
#cellstates = ct.SolutionArray(cellgas, extra=['t'])

elements = ['C', 'H', 'N', 'O']
element_masses = {el: 0.0 for el in elements}



# Initial conditions
equiv_ratio = 0.6
init_temp = 900
init_pr=ct.one_atm
zeta_0 = 0.25 # initial "unmixed fraction"


# initial state of the cell
C = ct.Quantity(gas)
C.TP = init_temp, init_pr
C.set_equivalence_ratio(equiv_ratio, {'CH4':1}, {'O2':1,'N2':3.76})
C.mass = 1.0
internal_energy_per_unit_mass = C.int_energy_mass
volume = C.v
print(C.Y)
for species_index, species in enumerate(gas.species_names):
    species_mass = C.Y[species_index] * C.mass  # Mass of the species in the mixture
    for el in elements:
        element_masses[el] += species_mass * gas.n_atoms(species, el) * gas.atomic_weight(el) / gas.molecular_weights[species_index]

print("\nElemental Balance (Total Mass in kg):")
for el, mass in element_masses.items():
    print(f" {el}: {mass:.6f} kg")
print(f"Sum of element mass = {sum(element_masses.values()):.6f} kg")


# print(C.report())
# sys.exit()

# combustion substep parameters
tau_mix = 1e-2
t_end = 10.0
dt_cfd = 1e-1
dt_sub_chem = 1e-4

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
        
        #print('Time: {m}'.format(m=t_cfd+t_sub))
        print('C mass: {m}'.format(m=M.mass+U.mass))
        print('U mass: {m}'.format(m=U.mass))
        print('M mass: {m}'.format(m=M.mass))
        print('C TMP: {m}'.format(m=C.T))
        print('U TMP: {m}'.format(m=U.T))
        print('M TMP: {m}'.format(m=M.T))
        print('C Pr: {m}'.format(m=C.P))
        print('U Pr: {m}'.format(m=U.P))
        print('M Pr: {m}'.format(m=M.P))
        
        zeta_last = zeta
        zeta = zeta_0*np.exp(-t_sub/tau_mix)
        print(f'b t_sub: {t_cfd+t_sub:.3f},U.T: {U.T:.3f}, U.P: {U.P:.3f}')
        
        print(f'b t_sub: {t_cfd+t_sub:.3f},M.v: {M.v:.3f}, M.volume: {M.volume:.3f}, M.T: {M.T:.3f}, M.P: {M.P:.3f}')
        J.mass = (zeta_last-zeta)*C.mass      # adjust mass of stream before adding
        M += J                                # mix U into M
        U.mass -= J.mass                      # reduce unmixed mass
        print(f'a t_sub: {t_cfd+t_sub:.3f},M.v: {M.v:.3f}, M.volume: {M.volume:.3f}, M.T: {M.T:.3f}, M.P: {M.P:.3f}')
        print(f'zeta_last: {zeta_last:.3f},zeta: {zeta:.3f}')
        #input("Press Enter to continue...")
        # define reactor using only the mixed portion of the cell
        gas.TPX = M.T, M.P, M.X
        r = ct.IdealGasConstPressureReactor(gas)
        sim = ct.ReactorNet([r])
        #sim.verbose = True
        sim.advance(dt_sub_chem)
        M.TPX = gas.TPX # now set the mixed zone based on reactor solution
       

    t_cfd += dt_cfd   
    print(f't_cfd: {t_cfd:.3f}, zeta_0: {zeta_0:.3f}, zeta: {zeta:.3f}')
    
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

    mixingstates.append(r.thermo.state, t=t_cfd*1e3)
    cellstates.append(C.state, 
                  t=t_cfd*1e3, 
                  C_mass=element_masses['C'], 
                  H_mass=element_masses['H'], 
                  N_mass=element_masses['N'], 
                  O_mass=element_masses['O'])
    print('Time: {m}'.format(m=t_cfd+t_sub))
    print('C.Volume: {m}'.format(m=C.v))
    print('C.T: {m}'.format(m=C.T))
    print('C.P: {m}'.format(m=C.P))
    #input("Press Enter to continue...")
    #print('U CH4: {m}'.format(m=U.Y[gas.species_index('CH4')]))
    #print('M CH4: {m}'.format(m=M.Y[gas.species_index('CH4')]))
    #print('C CH4: {m}'.format(m=C.Y[gas.species_index('CH4')]))
    #print('M TMP: {m}'.format(m=M.T))
    #print('C TMP: {m}'.format(m=C.T))

statesToSave=cellstates

columnNames = ['Time'] + ['OH'] + ['H']+['O2'] + ['CH4']+ ['C2H2'] + ['C_mass']+['H_mass']+['N_mass']+['O_mass']+['Temperature'] + ['Pressure'] 
statesForCSV = np.array([
    statesToSave.t, 
    statesToSave.Y[:, gas.species_index('OH')],
    statesToSave.Y[:, gas.species_index('H')],
    statesToSave.Y[:, gas.species_index('O2')],
    statesToSave.Y[:, gas.species_index('CH4')],
    statesToSave.Y[:, gas.species_index('C2H2')],
    statesToSave.C_mass,
    statesToSave.H_mass,
    statesToSave.N_mass,
    statesToSave.O_mass,
    statesToSave.T,
    statesToSave.P
]).T
csvdata = pd.DataFrame(statesForCSV, columns=columnNames)
save_csv_file(Chemistry_DIR+"ign_delay_Methane_grimech30_EDC_T900K_Phi0p6_soln.csv",csvdata)

