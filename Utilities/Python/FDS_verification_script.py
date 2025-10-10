#!$FIREMODELS/fds/.github/fds_python_env/bin/python
# McDermott
# 2 April 2024

import subprocess
import fdsplotlib
import importlib
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)
print("Using:", fdsplotlib.__file__)

# Scripts to run prior to dataplot

# print("ignition_delay...");   subprocess.run(["python","./scripts/cantera_ignition_delay.py"])
# print("reaction_rates...");   subprocess.run(["python","./scripts/cantera_reaction_rates.py"])
# print("turbulent_batch_reactor...");   subprocess.run(["python","./scripts/cantera_turbulent_batch_reactor.py"])

# Dataplot and scatplot options

# Statistics output options

# Run dataplot and scatplot scripts

fdsplotlib.dataplot(config_filename='../Matlab/FDS_verification_dataplot_inputs.csv',
                    expdir='../../Verification/',
                    cmpdir='../../Verification/',
                    pltdir='../../Manuals/',
                    close_figs=True,
                    verbose=True,
                    plot_range=[2,2]) # plot_range[start, end], optionally instead use plot_list['Dataname']

# Special cases

print("pyrolysis...");      subprocess.run(["python","./scripts/pyrolysis.py"])
print("turb_model...");     subprocess.run(["python","./scripts/turb_model.py"])

print("verification scripts completed successfully!")
