# FDS Validation

The sub-folders listed in this directory contain the input files for the FDS validation suite. These cases are run and processed with each minor release of FDS. The output from the validation simulations are contained in a separate repository [firemodels/out](https://github.com/firemodels/out). The experimental data is found in another repository [firemodels/exp](https://github.com/firemodels/exp). Publically available test reports may be found [here](https://drive.google.com/drive/folders/0B-EZ4HlrI6VDT2R5SjNFOGtIdTg).

Each sub-folder contains a bash script called `Run_All.sh` that creates a new folder called `Current_Results`, copies the FDS input files from the folder called `FDS_Input_Files`, and runs them. Another script called `Process_Output.sh` can be used to copy the necessary output files to the out repository.

There are two scripts in the current folder called `Run_Serial.sh` and `Run_Parallel.sh`. Each test series should be listed in one or the other. The first one lists all the series that contain only single process jobs. Typically, this script is run first because a large amount of the jobs are completed within a day, and these results act as a bellweather to identify problems that might have popped up since the last time the cases were run. `Run_Parallel.sh` lists the longer, multi-process jobs. This script is run after the serial cases are well along.

The script `Process_All_Output.sh` can be run at any time while the cases are running. This script scans the standard out and error files to determine if all the cases belonging to a given test series are completed; and if they are, copies the output to the `out` repository. All of the test series should be listed in this script.

The Python script called `Compile_Timings.py` scans the standard out files and generates a comma-delimited spreadsheet file that lists each completed case, its wall clock time, its number of processes, and the product of these two in units of CPU-hours.

Detailed instructions for running and processing the validation cases is found in the wiki called [FDS Validation Process](https://github.com/firemodels/fds/wiki/FDS-Validation-Process).
