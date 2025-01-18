import numpy as np
import os
import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd


# Helper function to correct mass fractions
def correct_mass_fractions(gas):
    mass_fractions = gas.Y
    mass_fractions[mass_fractions < 0] = 0.0  # Set negative values to 0
    total = np.sum(mass_fractions)
    if total > 0:
        mass_fractions /= total  # Normalize the mass fractions
    gas.Y = mass_fractions


# Helper function to correct mass fractions
def save_csv_file(fullFilePath, csvdata):
    csvdata = csvdata.apply(lambda x: x.str.strip() if x.dtype == 'object' else x)
    csvdata.replace('nan','',inplace=True)
    csvdata.replace('NAN','',inplace=True)
    csvdata = csvdata.map(lambda x: f'{x:10.5f}' if isinstance(x, (float)) else x)
    csvdata.to_csv(fullFilePath,index=False)


# Perform the reaction calculation
def calc_reaction_rate(mechanismFile, T0, P0, Y0, tstart, tend, dt, caseName, conditionalStop=False, condtionalSpeciesName=''):
    print("Running case:"+caseName)
    current_working_directory = os.getcwd()
    Cantera_DIR = os.path.join(current_working_directory, "../Input_Libraries/Chemical_Mechanisms/Cantera/")
    Species_DIR = os.path.join(current_working_directory, "../../Verification/Species/")
    
    gas = ct.Solution(os.path.join(Cantera_DIR, mechanismFile))
    gas.TPY = T0, P0, Y0
    reactor = ct.IdealGasReactor(gas)
    sim = ct.ReactorNet([reactor]) # Create a reactor network

    # Arrays to store time and temperature
    csvdata = pd.DataFrame()
    stateArr = []

    # Integrate the reactor
    time = tstart  # Initial time [s]
    while time < tend:
        time += dt
        if (conditionalStop):
            if (gas[condtionalSpeciesName].Y[0] >= 1e-20):
               sim.advance(time)
        else:
            sim.advance(time)
        stateArr.append([
            time,
            *gas.Y,
            gas.T - 273.15,
            gas.P - ct.one_atm
        ])

    tot_write_point=25
    writeInterval = max(1,int(tend/dt/tot_write_point))
    stateArrReduced = stateArr[::writeInterval]
    columnNames = ['Time'] + gas.species_names  + ['Temperature'] + ['Pressure'] 
    csvdata = pd.DataFrame(stateArrReduced, columns=columnNames)
    save_csv_file(Species_DIR+caseName+"_soln.csv",csvdata)



print("Running reactionrate_arrhenius cases...")

#------------------------------------------------------------------
# reactionrate_arrhenius_0order_1step
#------------------------------------------------------------------
caseName='reactionrate_arrhenius_0order_1step'
mechFile='one_step_CO2_formation.yaml'
T0 = 293.15  # Initial temperature [K]
P0 = ct.one_atm  # Initial pressure [Pa]
Y0 = {"CO": 0.289655, "O2": 0.165517,"N2": 0.544828}
tstart, tend, dt = 0.0, 10, 1e-3
conditionalStop=True
condtionalSpeciesName="O2"
calc_reaction_rate(mechFile, T0, P0, Y0, tstart, tend, dt, caseName, conditionalStop, condtionalSpeciesName)


#------------------------------------------------------------------
# reactionrate_arrhenius_2order_1step
#------------------------------------------------------------------
caseName='reactionrate_arrhenius_2order_1step'
mechFile='one_step_propane_oxidation.yaml'
T0 = 293.15  # Initial temperature [K]
P0 = ct.one_atm  # Initial pressure [Pa]
Y0 = {"C3H8": 0.06, "O2": 0.219,"N2": 0.721}
tstart, tend, dt = 0.0, 5, 0.05
calc_reaction_rate(mechFile, T0, P0, Y0, tstart, tend, dt, caseName)


#------------------------------------------------------------------
# reactionrate_arrhenius_1p75order_2step
#------------------------------------------------------------------
caseName='reactionrate_arrhenius_1p75order_2step'
mechFile='two_step_propane_oxidation.yaml'
T0 = 293.15  # Initial temperature [K]
P0 = ct.one_atm  # Initial pressure [Pa]
Y0 = {"C3H8": 0.06, "O2": 0.219,"N2": 0.721}
tstart, tend, dt = 0.0, 10, 0.05
calc_reaction_rate(mechFile, T0, P0, Y0, tstart, tend, dt, caseName)

#------------------------------------------------------------------
# reactionrate_arrhenius_1p75order_2stepr
#------------------------------------------------------------------
caseName='reactionrate_arrhenius_1p75order_2stepr'
mechFile='two_step_propane_oxidation_reverse.yaml'
T0 = 293.15  # Initial temperature [K]
P0 = ct.one_atm  # Initial pressure [Pa]
Y0 = {"C3H8": 0.06, "O2": 0.219,"N2": 0.721}
tstart, tend, dt = 0.0, 10, 0.05
calc_reaction_rate(mechFile, T0, P0, Y0, tstart, tend, dt, caseName)

#------------------------------------------------------------------
# reactionrate_arrhenius_equilibrium
#------------------------------------------------------------------
caseName='reactionrate_arrhenius_equilibrium'
mechFile='Westbrook_Dryer_Propane.yaml'
T0 = 350+273.15  # Initial temperature [K]
P0 = ct.one_atm  # Initial pressure [Pa]
Y0 = {"C3H8": 0.060321,"O2": 0.218851,"N2": 0.720828}
tstart, tend, dt = 0.0, 4, 0.00025
calc_reaction_rate(mechFile, T0, P0, Y0, tstart, tend, dt, caseName)


#------------------------------------------------------------------
# reactionrate_arrhenius_jones_lindstedt
#------------------------------------------------------------------
caseName='reactionrate_arrhenius_jones_lindstedt'
mechFile='Jones_Lindstedt.yaml'
T0 = 450+273.15  # Initial temperature [K]
P0 = ct.one_atm  # Initial pressure [Pa]
Y0 = {"C3H8": 0.060321,"O2": 0.218851,"H2O": 0.00722, "CO2": 0.000591, "N2": 0.713017}
tstart, tend, dt = 0.0, 5, 0.001
calc_reaction_rate(mechFile, T0, P0, Y0, tstart, tend, dt, caseName)
