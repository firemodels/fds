#!$FIREMODELS/fds/.github/fds_python_env/bin/python

import subprocess
import fdsplotlib
import importlib
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)

# Scripts to run prior to dataplot

#print("catchpole_spread_rates...");   subprocess.run(["python","./scripts/catchpole_spread_rates.py"])

# Dataplot and scatplot options

# Statistics output options

# # Run dataplot and scatplot scripts

# fdsplotlib.dataplot(config_filename='../Matlab/FDS_validation_dataplot_inputs.csv',
#                     expdir='../../../exp/',
#                     cmpdir='../../../out/',
#                     pltdir='../../Manuals/',
#                     close_figs=True,
#                     verbose=True,
#                     plot_list=['all'])

# Special cases

#print("Sandia_Pools...");       subprocess.run(["python","./scripts/Sandia_Pools.py"])
print("TUS_Facade...");         subprocess.run(["python","./scripts/TUS_Facade_contours.py"])
print("Wu_Bakar_Tunnels...");   subprocess.run(["python","./scripts/Wu_Bakar_Tunnels.py"])

print("Python validation scripts completed successfully!")
