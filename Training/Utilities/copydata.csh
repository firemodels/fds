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
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_01_hrr.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_02_hrr.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_03_hrr.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_04_hrr.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_05_hrr.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_06_hrr.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_07_hrr.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_01_devc.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_02_devc.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_03_devc.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_04_devc.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_05_devc.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_06_devc.csv $FDSDATA/.
cp ../MFRI/MFRI_Training_Tower/MFRI_Training_Tower_07_devc.csv $FDSDATA/.

cp ../Demonstrations/2Room_Ranch/ranch_01_hrr.csv $FDSDATA/.
cp ../Demonstrations/2Room_Ranch/ranch_02_hrr.csv $FDSDATA/.
cp ../Demonstrations/2Room_Ranch/ranch_03_hrr.csv $FDSDATA/.
cp ../Demonstrations/2Room_Ranch/ranch_04_hrr.csv $FDSDATA/.
cp ../Demonstrations/2Room_Ranch/ranch_01_devc.csv $FDSDATA/.
cp ../Demonstrations/2Room_Ranch/ranch_02_devc.csv $FDSDATA/.
cp ../Demonstrations/2Room_Ranch/ranch_03_devc.csv $FDSDATA/.
cp ../Demonstrations/2Room_Ranch/ranch_04_devc.csv $FDSDATA/.

cp ../MCFRS/MCFRS_Ranch/MCFRS_Ranch_00_hrr.csv $FDSDATA/.
cp ../MCFRS/MCFRS_Ranch/MCFRS_Ranch_00_devc.csv $FDSDATA/.
