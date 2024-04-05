#!$FIREMODELS/fds/.github/fds_python_env/bin/python3
# Chandan Paul
# 2 April 2024

import os
import numpy as np
import pandas as pd
import time
import cantera as ct
import matplotlib.pyplot as plt

# print(f"Runnning Cantera version: {ct.__version__}")

# get the current working directory
current_working_directory = os.getcwd()

# print output to the console
# print(current_working_directory)
Cantera_DIR = current_working_directory+"/../Input_Libraries/Chemical_Mechanisms/Cantera/"
Chemistry_DIR = current_working_directory+"/../../Verification/Chemistry/"

gas = ct.Solution(Cantera_DIR+"grimech30.yaml")
reactor_temperature = 1273.15  # Kelvin
reactor_pressure = 101325  # Pascals
reference_species = "OH"

gas.TP = reactor_temperature, reactor_pressure

# Define the fuel, oxidizer and set the stoichiometry
gas.set_equivalence_ratio(phi=1.0, fuel="CH4", oxidizer={"O2": 1.0, "N2": 3.76})
# print("Init Y=",gas.Y)

# Create a batch reactor object and add it to a reactor network
# In this example, the batch reactor will be the only reactor
# in the network
r = ct.IdealGasReactor(contents=gas, name="Batch Reactor")
reactor_network = ct.ReactorNet([r])

# use the above list to create a DataFrame
time_history = ct.SolutionArray(gas, extra="t")

def ignition_delay(states, species):
    """
    This function computes the ignition delay from the occurence of the
    peak in species' concentration.
    """
    i_ign = states(species).Y.argmax()
    return states.t[i_ign]

# print("State=",r.thermo.state)
# print("T=",gas.TP)
# print("Y=",gas.Y)
#print("P=",r.thermo.state.P)
#print("O2=",r.thermo.statei("O2").Y)

# Tic
t0 = time.time()

# This is a starting estimate. If you do not get an ignition within this time, increase it
estimated_ignition_delay_time = 0.1
t = 0

counter = 1
while t < estimated_ignition_delay_time:
    t = reactor_network.step()
    if not counter % 10:
        # We will save only every 10th value. Otherwise, this takes too long
        # Note that the species concentrations are mass fractions
        time_history.append(r.thermo.state, t=t)
    counter += 1

# We will use the 'OH' species to compute the ignition delay
tau = ignition_delay(time_history, reference_species)

# Toc
t1 = time.time()

# print(f"Computed Ignition Delay: {tau:.3e} seconds. Took {t1-t0:3.2f}s to compute")

fig, (ax1, ax2) = plt.subplots(2) 

ax1.plot(time_history.t, time_history(reference_species).Y, "-o")
ax1.set(xlabel='Time (s)', ylabel='$Y_{OH}$')
ax1.set_xlim([0, 0.05])
ax1.arrow(
    0,
    0.008,
    tau,
    0,
    width=0.0001,
    head_width=0.0005,
    head_length=0.001,
    length_includes_head=True,
    color="r",
    shape="full",
)
ax1.annotate(
    r"$Ignition Delay: \tau_{ign}$", xy=(0, 0), xytext=(0.01, 0.0082), fontsize=16
);

ax2.plot(time_history.t, time_history.T, "-o")
ax2.set(xlabel='Time (s)', ylabel='$T (K)$')
ax2.set_xlim([0, 0.05])

# plt.savefig('ignitionDelay_fixedTemp.pdf')

# Write data to csv file for post-processing

df = pd.DataFrame({'Time (s)': time_history.t, 'T (K)': time_history.T})

df.to_csv(Chemistry_DIR+"ignition_delay_Cantera_results.csv", sep=',', header=True, index=False, index_label=False)



