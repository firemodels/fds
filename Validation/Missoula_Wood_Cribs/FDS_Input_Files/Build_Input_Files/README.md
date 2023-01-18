# Missoula Wood Cribs FDS Input Files

This directory contains an fds file template and a python script to build fds input files.

The data for the input files is taken from Table 1 of

Sara McAllister and Mark Finney. Burning Rates of Wood Cribs with Implications for Wildland Fires. Fire Technology, 52, 1755-1777, 2016.


## Mass Producing Input Files from a Template

The python script `build_input_files.py` reads `paramfile.csv` and runs `swaps.py` to build the input files.

## Usage

    $ python build_input_files.py

