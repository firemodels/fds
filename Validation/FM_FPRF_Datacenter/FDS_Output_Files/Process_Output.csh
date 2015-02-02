#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/FM_FPRF_Datacenter/FDS_Output_Files
set DDIR=$SVNROOT/Validation/FM_FPRF_Datacenter/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR


