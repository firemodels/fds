# FDS Verification

The folders above contain the FDS verification input files that are run each night by the firebot script contained in the repository [firemodels/bot](https://github.com/firemodels/bot). Detailed instructions for running and processing the verification cases are found in the wiki [FDS Verification Process](https://github.com/firemodels/fds/wiki/FDS-Verification-Process).

For those making substantial changes to the FDS source code, it is a good idea to run the verification cases in debug mode for a few time steps to quickly identify easy-to-make and easy-to-fix bugs, like array overshoots and divide by zeros. To do this, follow these steps:
1. Compile the debug versions of FDS, with and without OpenMP.
2. cd to the `Verification` folder and do a `git clean -dxf` which will erase all uncommitted files within `Verification` and its subfolders, including old output files. Do not do this if there are files that you have not yet committed. Save these somewhere else.
3. cd to the `scripts` subfolder and issue the command `./Run_FDS_Cases.sh -m 2 -d -q firebot`. The `-m` option instructs FDS to only run the cases for 2 time steps. The `-d` option tells FDS to use the debug version of FDS. The `-q` option is the queue to use. At NIST, we use the `firebot` queue for this kind of quick testing.
4. cd back up to the `Verification` folder and issue the command `grep forrtl */*err` which will search all the diagnostic output files for Fortran run-time errors. Deal with these and then repeat.
5. When all the debug cases have run successfully, you can compile the release versions of FDS (with and without OpenMP), clean the `Verication` folder again, cd to the `scripts` folder, and issue `./Run_FDS_Cases.sh -q firebot` which will run the cases in release mode. This takes a few hours. 
6. When the cases are done, run a few cases to test the restart feature with the command `./Run_FDS_Cases.sh -r`. These few cases take less than a minute to run.
7. Run the Matlab script in `Utilities/Matlab` called `FDS_verification_script.m`. This takes 10 to 15 minutes
8. When the Matlab script completes, cd to `Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots` and open `verification_scatterplot_output.csv` using a spreadsheet editor and look for cases that have failed. Correct these and then repeat the process.
