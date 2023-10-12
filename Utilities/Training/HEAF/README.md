# High Energy Arc Fault Simulations

The subdirectories listed above contain FDS input files for simulations of high energy arc faults within electrical enclosures. The results of these simulations are reported in [NUREG-2262](https://www.nrc.gov/docs/ML2310/ML23108A113.pdf). 

Each subdirectory contains a single FDS input file that serves as a template for generating the entire set. To generate all the input files, use the command
```
python swaps.py
```
where `swaps.py` is a Python script found in `Utilities/Input_File_Tools`. 
