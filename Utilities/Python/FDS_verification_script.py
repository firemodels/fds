#!$FIREMODELS/fds/.github/fds_python_env/bin/python
# McDermott
# 2 April 2024

import subprocess
import fdsplotlib
import importlib
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)

# Scripts to run prior to dataplot

print("ignition_delay...");   subprocess.run(["python","./scripts/cantera_ignition_delay.py"])

# Dataplot and scatplot options

# Statistics output options

# Run dataplot and scatplot scripts

fdsplotlib.dataplot(config_filename='FDS_verification_dataplot_config.csv',
                    revision='test',
                    expdir='../../Verification/',
                    cmpdir='../../Verification/',
                    pltdir='../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/',
                    close_figs=True,
                    verbose=True,
                    plot_list=['all'])

# Special cases

print("verification scripts completed successfully!")
