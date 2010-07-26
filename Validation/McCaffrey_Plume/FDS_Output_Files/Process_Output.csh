#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/McCaffrey_Plume/FDS_Output_Files
set DDIR=$SVNROOT/Validation/McCaffrey_Plume/Current_Results
cd $DDIR
$PDIR/fds2ascii < $WDIR/McCaffrey_14_kW_fds2ascii.input
$PDIR/fds2ascii < $WDIR/McCaffrey_22_kW_fds2ascii.input
$PDIR/fds2ascii < $WDIR/McCaffrey_33_kW_fds2ascii.input
$PDIR/fds2ascii < $WDIR/McCaffrey_45_kW_fds2ascii.input
$PDIR/fds2ascii < $WDIR/McCaffrey_57_kW_fds2ascii.input
cp McC*fds2ascii.csv $WDIR

