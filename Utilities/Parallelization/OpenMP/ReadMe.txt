ReadMe for Parallelization Test "Suite" (PTS)
(C) C. Rogsch

The PTS is based on a Linux GUI scripting interface (zenity, runs without problems on Ubuntu). Aim of this PTS is to test the correct parallelization with OpenMP. Tests are based on the Intel ThreadChecker. Some FAQ are listed below, if anything is not clear, do not hesitate to contact...

######## -------- ######## -------- ######## -------- ######## -------- ######## -------- ########

How to start PTS?
PTS should be started via commannd line (shell), type ./Start_Tests.sh >& Screen.txt &
It is recommended to write all screen outputs to a file to ensure/verify the correct compilation process. Actually, no automatic check of correct compilation is implemented, thus you have to check this visually by looking at the Screen.txt file.

######## -------- ######## -------- ######## -------- ######## -------- ######## -------- ########

What features are implemented, how to proceed?
After starting PTS the user can choose which kind of FDS-Version should be compiled. There are different choices about OpenMP-builds, thus the user can choose which files should be compiled with OpenMP/Intel ThreadChecker-flags. This is important, because an FDS OpenMP/ThreadChecker run is very time consuming. If only subroutines are edited it is easier to check only this subroutines instead of the complete code.

Different FDS-Versions can be compiled and tested in immediate succession. Different testfiles can be chosen, too. This means that all selected FDS-Versions will run all selected testfiles. At the end of the complete testrun, a message will appear (Calculation finished). During calculation in each Test_Cases subdirectory output files are generated (txt-files) based on testcase name and FDS-Version. After each run an Outputfile with Intel ThreadChecker results is generated, too. Open this files and check if any error-messages are included!!!

######## -------- ######## -------- ######## -------- ######## -------- ######## -------- ########

Which files should I choose for testing?
Compare update FDS-Source files with the last edition of FDS-Source files, where OpenMP edits are done. Subroutines, which are changed should be tested. Have a look at the "DEBUG_OpenMP_Overview" file. This will be added soon...


######## -------- ######## -------- ######## -------- ######## -------- ######## -------- ########

Where can I find compiled FDS-Versions?
All compiled FDS-Versions are copied to FDS_OpenMP-Versions directory.

######## -------- ######## -------- ######## -------- ######## -------- ######## -------- ########

Why is PTS so storage consuming?
That is because you do not delete old files or FDS output files (SliceFiles, PL3D-Files...) Also, each compiled FDS_OpenMP-Version is copied to the selected Test_Cases directories. Select "Delete files after run" in the last check-box of PTS, than all files are deleted after the calculation, except the .txt-files for screen output and ThreadChecker output, but they are not very storage consuming. If you do not need them any more, delete them "by hand". The option "not to delete files" is added to check the calculation in SmokeView, e. g. to analyse results "visually".

######## -------- ######## -------- ######## -------- ######## -------- ######## -------- ########

How to add new test cases?
Create in folder "Test_Cases" a new subdirectory, e. g. Test000
Copy the new file in this directory, e.g. Test000.fds
Name file and folder in the same way (like directory = Test000, file = Test000.fds), do not use any spaces in the name

IMPORTANT: DO NOT COPY MORE THAN ONE TESTFILE IN EACH SUBDIRECTORY OF TEST_CASES!!! CHECK THAT DIRECTORY AND FILENAME (without .fds) HAVE THE SAME NAME!!!

Edit file Start_Test.sh
 - go to function "select_calculation_filename"
   - add a new line to this function, see Start_Test.sh file
 - go to function "run_calculation"
   - add an entry like for each new file
        Test000)       -> Edit this to new file-/foldername without .fds
          run_fds      -> Do not edit
          ;;           -> do not edit

Finally, test it if it works...! You do not need do do an OpenMP/ThreadChecker calculation, just a calculation to see if the file is selected correctly. Look at your Screen.txt file! If it runs successfully, edit "DEBUG_OpenMP_Overview" file.
