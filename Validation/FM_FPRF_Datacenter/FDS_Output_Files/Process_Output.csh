#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/FM_FPRF_Datacenter/FDS_Output_Files
set DDIR=$FDSSMV/Validation/FM_FPRF_Datacenter/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR


