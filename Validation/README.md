# FDS Validation

The folders above contain the FDS input files for the FDS validation suite. These cases are run and processed with each minor release of FDS. The output from the validation simulations are contained in a separate repository [firemodels/out](https://github.com/firemodels/out). The experimental data is found in another repository [firemodels/exp](https://github.com/firemodels/exp). Publically available test reports may be found [here](https://drive.google.com/drive/folders/0B-EZ4HlrI6VDT2R5SjNFOGtIdTg).

Each folder contains a bash script called `Run_All.sh` that creates a new folder called `Current_Results`, copies the FDS input files from the folder called `FDS_Input_Files`, and runs them. Another script called `Process_Output.sh` can be used to copy the necessary output files to the out repository.

Detailed instructions for running and processing the validation cases is found in the wiki called [FDS Validation Process](https://github.com/firemodels/fds/wiki/FDS-Validation-Process).
