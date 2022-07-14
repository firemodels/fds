# Theobald Hose Stream fds Input Files

This directory contains experimental parameters, an fds file template, and a python script to build fds input files.

## Mass Producing Input Files from a Template

The python script `build_spray_input_files.py` reads `theobald_effect_1981_fds.csv`, wich contains the paramaters for the 90 spray and converts them into paramaters to be inserted into `theobald_Template.fds`, making the file `paramfile.csv`. `build_spray_input_files.py` then calls `swaps.py` from the fds Utilties directory, wich uses `theobald_Template.fds` and `paramfile.csv` to generate the fds input files. 

Each fds input file is numbered based on it's row index in `theobald_effect_1981_fds.csv`