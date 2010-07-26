#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/Restivo_Experiment/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Restivo_Experiment/Current_Results
cd $DDIR
$PDIR/fds2ascii < $WDIR/Restivo_fds2ascii_hori.input
$PDIR/fds2ascii < $WDIR/Restivo_fds2ascii_vert.input
cp Res*fds2ascii_hori.csv $WDIR
cp Res*fds2ascii_vert.csv $WDIR

