#!$FIREMODELS/fds/.github/fds_python_env/bin/python

import subprocess
import fdsplotlib
import importlib
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)

# Scripts to run prior to dataplot

print("catchpole_spread_rates...");   subprocess.run(["python","./scripts/catchpole_spread_rates.py"])

# Dataplot and scatplot options

# Statistics output options

# Run dataplot and scatplot scripts

# fdsplotlib.dataplot(config_filename='../Matlab/FDS_verification_dataplot_inputs.csv',
#                     expdir='../../Verification/',
#                     cmpdir='../../Verification/',
#                     pltdir='../../Manuals/',
#                     close_figs=True,
#                     verbose=True,
#                     plot_range=[2])

# Special cases

print("Python validation scripts completed successfully!")
