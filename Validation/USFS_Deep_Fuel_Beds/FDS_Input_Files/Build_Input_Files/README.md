# USFS Deep Fuel Beds fds Input Files

This directory contains an fds file template and a python script to build fds input files using experimental parameters from the exp repository.

## Mass Producing Input Files from a Template

The python script `build_input_files.py` reads `exp_params.csv` in the exp repository, wich contains the paramaters for the 108 experamental burns and converts them into paramaters to be inserted into `Template.fds`, making the file `paramfile.csv`. `build_input_files.py` then calls `swaps.py` from the fds Utilties directory, wich uses `Template.fds` and `paramfile.csv` to generate the fds input files and then moves these files up one level to FDS_Input_Files.

## Usage

    $ python build_input_files.py

