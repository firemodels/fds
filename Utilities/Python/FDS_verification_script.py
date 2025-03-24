#!$FIREMODELS/fds/.github/fds_python_env/bin/python
# McDermott
# 2 April 2024

import subprocess
import fdsplotlib
import importlib
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)

# Scripts to run prior to dataplot

# print("ignition_delay...");   subprocess.run(["python","./scripts/cantera_ignition_delay.py"])
# print("reaction_rates...");   subprocess.run(["python","./scripts/cantera_reaction_rates.py"])

# Dataplot and scatplot options

# Statistics output options

# Run dataplot and scatplot scripts

fdsplotlib.dataplot(config_filename='../Matlab/FDS_verification_dataplot_inputs.csv',
                    expdir='../../Verification/',
                    cmpdir='../../Verification/',
                    pltdir='../../Manuals/',
                    close_figs=True,
                    verbose=True,
                    plot_range=[2])

# Special cases

print("pyrolysis...");   subprocess.run(["python","./scripts/pyrolysis.py"])

print("verification scripts completed successfully!")
