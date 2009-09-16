#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/Bryant_Doorway/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Bryant_Doorway/Current_Results
cd $DDIR
$PDIR/fds2ascii $WDIR/Bryant_034_kW_fds2ascii.input
$PDIR/fds2ascii $WDIR/Bryant_065_kW_fds2ascii.input
$PDIR/fds2ascii $WDIR/Bryant_096_kW_fds2ascii.input
$PDIR/fds2ascii $WDIR/Bryant_128_kW_fds2ascii.input
$PDIR/fds2ascii $WDIR/Bryant_160_kW_fds2ascii.input
$PDIR/fds2ascii $WDIR/Bryant_320_kW_fds2ascii.input
$PDIR/fds2ascii $WDIR/Bryant_511_kW_fds2ascii.input
cp Br*fds2ascii.csv $WDIR
cd $WDIR


