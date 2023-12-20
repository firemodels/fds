# High Energy Arc Fault Simulations

The subdirectories listed above contain FDS input files for simulations of high energy arc faults within electrical enclosures. The results of these simulations are reported in [NUREG-2262](https://www.nrc.gov/docs/ML2310/ML23108A113.pdf). 

Each subdirectory contains a single FDS input file that serves as a template for generating the entire set. To generate all the input files, `cd` into the appropriate subdirectory and issue the command
```
python ../../../Input_File_Tools/swaps.py
```
where `swaps.py` is a Python script found in `Utilities/Input_File_Tools`. The script opens the file called `paramfile.csv` and reads the name of the file in the first row and column. This is the template for all of the input files. The parameter values in `paramfile.csv` are "swapped" for the dummy arguments in the template.

Note that these simulations are fragile; that is, they are vulnerable to numerical instability because of the rapid release of energy from a concentrated source of energy.
