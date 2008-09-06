#!/bin/csh -f

# specify location of repository root
setenv SVNROOT ~/FDS-SMV

# VVVVVVVVVVV Do not change these line VVVVVVVVVVVVVV
setenv BASEDIR `pwd`
setenv FDSDATA ../FDS_Output_Files
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

cp ../MCFRS/MCFRS_Flashover/MCFRS_Flashover_00_devc.csv $FDSDATA/.
cp ../MCFRS/MCFRS_Flashover/MCFRS_Flashover_00_hrr.csv $FDSDATA/.
cp ../MCFRS/MCFRS_Flashover/MCFRS_Flashover_01_devc.csv $FDSDATA/.
cp ../MCFRS/MCFRS_Flashover/MCFRS_Flashover_01_hrr.csv $FDSDATA/.
cp ../MCFRS/MCFRS_Flashover/MCFRS_Flashover_02_devc.csv $FDSDATA/.
cp ../MCFRS/MCFRS_Flashover/MCFRS_Flashover_02_hrr.csv $FDSDATA/.
cp ../MCFRS/MCFRS_Flashover/MCFRS_Flashover_03_devc.csv $FDSDATA/.
cp ../MCFRS/MCFRS_Flashover/MCFRS_Flashover_03_hrr.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_00_hrr.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_00_devc.csv $FDSDATA/.
cp ../Demonstrations/2Room_Ranch/ranch_01_hrr.csv $FDSDATA/.
cp ../Demonstrations/2Room_Ranch/ranch_02_hrr.csv $FDSDATA/.
cp ../Demonstrations/2Room_Ranch/ranch_03_hrr.csv $FDSDATA/.
cp ../Demonstrations/2Room_Ranch/ranch_04_hrr.csv $FDSDATA/.
