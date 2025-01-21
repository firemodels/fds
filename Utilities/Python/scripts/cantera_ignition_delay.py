import numpy as np
import time
import os
import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd


def calc_ignition_delay(gas, equivRatios, T, fuel,estimated_ignition_delay_times, csvdata, caseCount):
    reactor_pressure = 101325  # Pascals
    reference_species = "OH"
    writeInterval = 20
    
    for phi in equivRatios:
        # Now create a SolutionArray out of these
        ignition_delays = ct.SolutionArray(gas, shape=T.shape, extra={"tau": estimated_ignition_delay_times})
        
        ignition_delays.set_equivalence_ratio(phi, fuel=fuel, oxidizer={"O2": 1.0, "N2": 3.76})
        #print("Y=",gas.Y) # To know the mass fractions
        ignition_delays.TP = T, reactor_pressure
    
        for i, state in enumerate(ignition_delays):
            caseCount = caseCount +1
            stateArr = []
            stateArrReduced = []
            # Setup the gas and reactor
            gas.TPX = state.TPX
            r = ct.IdealGasReactor(contents=gas, name="Batch Reactor")
            reactor_network = ct.ReactorNet([r])
        
            reference_species_history = []
            time_history = []
        
            t0 = time.time()
        
            t = 0
            while t < estimated_ignition_delay_times[i]:
                t = reactor_network.step()
                time_history.append(t)
                reference_species_history.append(gas[reference_species].X[0])
                tCust = t
                if (tCust > estimated_ignition_delay_times[i]):
                    tCust = estimated_ignition_delay_times[i]
                stateArr.append([tCust, gas[reference_species].Y[0],gas.T-273.15 ]) 
            i_ign = np.array(reference_species_history).argmax()
            tau = time_history[i_ign]
            t1 = time.time()
            
            # Reduce number of timesteps to reduce file size
            # Before ignition delay+0.1 s write after every writeInterval variable 
            # After that write all the times.
            tCount = 0
            for j in range(len(stateArr)):
                t= stateArr[j][0]
                if (t < tau + 0.1):
                    if (tCount%writeInterval ==0):
                        stateArrReduced.append(stateArr[j])
                else:
                    stateArrReduced.append(stateArr[j])
                tCount = tCount+1
            stateArrDF = pd.DataFrame(stateArrReduced)
            
            caseIndx = str(caseCount)
            if caseCount ==1:
                csvdata = stateArrDF
            else:    
                #Pad the smaller DataFrame with NaNs to match the number of rows
                if csvdata.shape[0] > stateArrDF.shape[0]:
                    pad_rows = csvdata.shape[0] - stateArrDF.shape[0]
                    stateArrDF = pd.concat([stateArrDF, pd.DataFrame(np.full((pad_rows, stateArrDF.shape[1]), np.nan))], axis=0, ignore_index=True)
                elif stateArrDF.shape[0] > csvdata.shape[0]:
                    pad_rows = stateArrDF.shape[0] - csvdata.shape[0]
                    csvdata = pd.concat([csvdata, pd.DataFrame(np.full((pad_rows, csvdata.shape[1]), np.nan))], axis=0, ignore_index=True)
    
                csvdata = pd.concat([csvdata, stateArrDF], axis=1, ignore_index=True)
                
            # print(f"Computed Ignition Delay: {tau:.3e} seconds for phi={phi} T={ignition_delays[i].T}K. Took {t1 - t0:3.2f}s to compute")
        
            ignition_delays[i].tau = tau
            
    return csvdata, caseCount        

# End of calc_ignition_delay function


current_working_directory = os.getcwd()
Cantera_DIR = current_working_directory+"/../Input_Libraries/Chemical_Mechanisms/Cantera/"
Chemistry_DIR = current_working_directory+"/../../Verification/Chemistry/"
csvdata = pd.DataFrame()
caseCount = 0

# Metahne_grimech 12 
print("Solving for Metahne_grimech:")
gas = ct.Solution(Cantera_DIR+"Methane_grimech30.yaml")
equivRatios = np.array([0.6, 1.0,  1.4])
T = np.array([900, 1000, 1100, 1200])
estimated_ignition_delay_times = np.ones_like(T, dtype=float)
estimated_ignition_delay_times[:] = 10
csvdata, caseCount=calc_ignition_delay(gas, equivRatios, T, "CH4", estimated_ignition_delay_times, csvdata, caseCount)
print("Done Metahne_grimech, caseCount=", str(caseCount))

# Methane_TianfengLu
print("Solving for Methane_TianfengLu:")
gas = ct.Solution(Cantera_DIR+"Methane_TianfengLu.yaml")
equivRatios = np.array([1.0])
T = np.array([1100])
estimated_ignition_delay_times = np.ones_like(T, dtype=float)
estimated_ignition_delay_times[:] = 1
csvdata, caseCount=calc_ignition_delay(gas, equivRatios, T, "CH4", estimated_ignition_delay_times, csvdata, caseCount)
print("Done Methane_TianfengLu, caseCount=", str(caseCount))

# Methane_Smooke
print("Solving for Methane_Smooke:")
gas = ct.Solution(Cantera_DIR+"Methane_Smooke.yaml")
equivRatios = np.array([1.0])
T = np.array([1100])
estimated_ignition_delay_times = np.ones_like(T, dtype=float)
estimated_ignition_delay_times[:] = 1
csvdata, caseCount=calc_ignition_delay(gas, equivRatios, T, "CH4", estimated_ignition_delay_times, csvdata, caseCount)
print("Done Methane_Smooke, caseCount=", str(caseCount))

# Ethylene_TianfengLu
print("Solving for Ethylene_TianfengLu:")
gas = ct.Solution(Cantera_DIR+"Ethylene_TianfengLu.yaml")
equivRatios = np.array([1.0])
T = np.array([1100])
estimated_ignition_delay_times = np.ones_like(T, dtype=float)
estimated_ignition_delay_times[:] = 1
csvdata, caseCount=calc_ignition_delay(gas, equivRatios, T, "C2H4", estimated_ignition_delay_times, csvdata, caseCount)
print("Done Ethylene_TianfengLu, caseCount=", str(caseCount))

# Propane_USC
print("Solving for Propane_USC:")
gas = ct.Solution(Cantera_DIR+"Propane_USC.yaml")
equivRatios = np.array([1.0])
T = np.array([1100])
estimated_ignition_delay_times = np.ones_like(T, dtype=float)
estimated_ignition_delay_times[:] = 1
csvdata, caseCount=calc_ignition_delay(gas, equivRatios, T, "C3H8", estimated_ignition_delay_times, csvdata, caseCount)
print("Done Propane_USC, caseCount=", str(caseCount))

# Propane_Z66
print("Solving for Propane_Z66:")
gas = ct.Solution(Cantera_DIR+"Propane_Z66.yaml")
equivRatios = np.array([1.0])
T = np.array([1100])
estimated_ignition_delay_times = np.ones_like(T, dtype=float)
estimated_ignition_delay_times[:] = 1
csvdata, caseCount=calc_ignition_delay(gas, equivRatios, T, "C3H8", estimated_ignition_delay_times, csvdata, caseCount)
print("Done Propane_Z66, caseCount=", str(caseCount))

# nHeptane_Chalmers
print("Solving for nHeptane_Chalmers:")
gas = ct.Solution(Cantera_DIR+"nHeptane_Chalmers.yaml")
equivRatios = np.array([1.0])
T = np.array([1100])
estimated_ignition_delay_times = np.ones_like(T, dtype=float)
estimated_ignition_delay_times[:] = 1
csvdata, caseCount=calc_ignition_delay(gas, equivRatios, T, "C7H16", estimated_ignition_delay_times, csvdata, caseCount)
print("Done nHeptane_Chalmers, caseCount=", str(caseCount))

# Set Dataframe column headers and writing format
timeFormat = '{:10.5f}'
OHFormat = '{:.3E}'
TMPFormat = '{:10.2f}'
colPerCase = 3 # Time, OH, TMP
for i in range(caseCount):
    caseIndx = str(i+1)
    csvdata = csvdata.rename(columns={i*colPerCase: 'Time'+caseIndx, i*colPerCase+1: 'OH'+caseIndx, i*colPerCase+2: 'TMP'+caseIndx})
    csvdata['Time'+caseIndx] = csvdata['Time'+caseIndx].map(lambda x: timeFormat.format(x))
    csvdata['OH'+caseIndx] = csvdata['OH'+caseIndx].map(lambda x: OHFormat.format(x))
    csvdata['TMP'+caseIndx] = csvdata['TMP'+caseIndx].map(lambda x: TMPFormat.format(x))

csvdata = csvdata.apply(lambda x: x.str.strip() if x.dtype == 'object' else x)
csvdata.replace('nan','',inplace=True)
csvdata.replace('NAN','',inplace=True)
csvdata.to_csv(Chemistry_DIR+"cantera_ignition_delay.csv",index=False)

