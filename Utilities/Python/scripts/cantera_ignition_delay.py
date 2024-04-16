import numpy as np
import time
import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd

print(f"Runnning Cantera version: {ct.__version__}")

gas = ct.Solution("./mechanisms/GRIMECH/grimech30.yaml")
reactor_pressure = 101325  # Pascals
reference_species = "OH"

# Make a list of all the equivalence ratios and temperatures we would like to run simulations at
equivRatios = np.array([0.6, 1.0,  1.4])
T = np.array([900, 1000, 1100, 1200])


csvdata = pd.DataFrame()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylabel("Ignition Delay (s)")
ax.set_xlabel(r"$\frac{10000}{T (K)}$", fontsize=18)
    
caseCount = 0

for phi in equivRatios:
    estimated_ignition_delay_times = np.ones_like(T, dtype=float)
    estimated_ignition_delay_times[:] = 10

    # Now create a SolutionArray out of these
    ignition_delays = ct.SolutionArray(
        gas, shape=T.shape, extra={"tau": estimated_ignition_delay_times}
    )
    
    ignition_delays.set_equivalence_ratio(
        phi, fuel="CH4", oxidizer={"O2": 1.0, "N2": 3.76}
    )
    #print("Y=",gas.Y) # To know the mass fractions
    ignition_delays.TP = T, reactor_pressure
    

    for i, state in enumerate(ignition_delays):
        caseCount = caseCount +1
        stateArr = []
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
            stateArr.append([t, gas[reference_species].X[0],gas.T ]) 
    
        i_ign = np.array(reference_species_history).argmax()
        tau = time_history[i_ign]
        t1 = time.time()
        stateArrDF = pd.DataFrame(stateArr)
        

        
        caseIndx = str(caseCount)
        if caseCount ==1:
            csvdata = stateArrDF
        else:    
            # Pad the smaller DataFrame with NaNs to match the number of rows
            if csvdata.shape[0] > stateArrDF.shape[0]:
                pad_rows = csvdata.shape[0] - stateArrDF.shape[0]
                stateArrDF = pd.concat([stateArrDF, pd.DataFrame(np.full((pad_rows, stateArrDF.shape[1]), np.nan))], axis=0, ignore_index=True)
            elif stateArrDF.shape[0] > csvdata.shape[0]:
                pad_rows = stateArrDF.shape[0] - csvdata.shape[0]
                csvdata = pd.concat([csvdata, pd.DataFrame(np.full((pad_rows, csvdata.shape[1]), np.nan))], axis=0, ignore_index=True)
            csvdata = pd.concat([csvdata, stateArrDF], axis=1, ignore_index=True)
            
        #csvdata = csvdata.rename(columns={0: 'Time'+caseIndx, 1: 'OH'+caseIndx, 2: 'TMP'+caseIndx})
        print(
            f"Computed Ignition Delay: {tau:.3e} seconds for phi={phi} T={ignition_delays[i].T}K. Took {t1 - t0:3.2f}s to compute"
        )
    
        ignition_delays[i].tau = tau
        

    ax.semilogy(10000/ ignition_delays.T, ignition_delays.tau, "o-",label = "Phi = {}".format(phi))
   
# Set Dataframe column headers
colPerCase = 3 # Time, OH, TMP
for i in range(caseCount):
    caseIndx = str(i+1)
    csvdata = csvdata.rename(columns={i*colPerCase: 'Time'+caseIndx, i*colPerCase+1: 'OH'+caseIndx, i*colPerCase+2: 'TMP'+caseIndx})

csvdata.to_csv('../../../Verification/Chemistry/cantera_ignition_delay.csv',index=False)

# Show plot legend
plt.legend() 

# Add a second axis on top to plot the temperature for better readability
ax2 = ax.twiny()
ticks = np.hstack((np.arange(800, 1300, 100))) #ax.get_xticks()
#ticks = ax.get_xticks()
ax2.set_xticks((10000/ ticks))
ax2.set_xticklabels(ticks)
ax2.set_xlim(ax.get_xlim())
ax2.set_xlabel("Temperature: $T(K)$");


# plt.savefig('ignDelay_temp_sweep.png')
plt.show();

