#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/Ulster_SBI/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Ulster_SBI/Current_Results
cd $DDIR
$PDIR/Ulster_SBI_fds2ascii < $WDIR/Ulster_SBI_30_kW_fds2ascii.input
$PDIR/Ulster_SBI_fds2ascii < $WDIR/Ulster_SBI_45_kW_fds2ascii.input
$PDIR/Ulster_SBI_fds2ascii < $WDIR/Ulster_SBI_60_kW_fds2ascii.input
cp Ulster*fds2ascii.csv $WDIR
cd $WDIR

