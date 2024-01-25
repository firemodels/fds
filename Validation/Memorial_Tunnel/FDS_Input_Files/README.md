## Notes for Preparing Memorial Tunnel Input Files

There are 98 Memorial Tunnel experiments. The first 81 are referred to as the "tranvserse ventilation" tests; the last 17 the "longitudinal ventilation" tests. The former employ a flat ceiling approximately 4.4 m high. 
The latter replace the flat ceiling with jet fans under the semi-circular vaulted ceiling. 

The input files for the transverse cases are constructed using two python scripts, which are run in sequence:
```
python ~/firemodels/fds/Utilities/Input_File_Tools/swaps.py
python write_ramps.py
```
The first of these scripts reads the skeleton file called `Template.fds` and replaces certain character strings with parameters listed in the file called `paramfile.csv`. The second script reads the measured heat release and 
ventilation rates from the experimental files contained in the folder `firemodels/exp/Memorial_Tunnel` and adds `RAMP`s to the input file. Many of the input parameters for the various sets of experiments are contained in the `.txt` files that 
are catenated into the input files using the `&CATF` feature. Thus, when the cases are run, the input and output files will be appended with the string `cat` to indicate that these are the complete, catenated files.
