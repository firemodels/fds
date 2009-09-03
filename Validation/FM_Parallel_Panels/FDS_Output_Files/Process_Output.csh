#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/FM_Parallel_Panels/FDS_Output_Files
set DDIR=$SVNROOT/Validation/FM_Parallel_Panels/Current_Results
cd $DDIR
$PDIR/fds2ascii $WDIR/FM_Parallel_Panel_1_fds2ascii.input
$PDIR/fds2ascii $WDIR/FM_Parallel_Panel_2_fds2ascii.input
$PDIR/fds2ascii $WDIR/FM_Parallel_Panel_3_fds2ascii.input
$PDIR/fds2ascii $WDIR/FM_Parallel_Panel_4_fds2ascii.input
$PDIR/fds2ascii $WDIR/FM_Parallel_Panel_5_fds2ascii.input
$PDIR/fds2ascii $WDIR/FM_Parallel_Panel_6_fds2ascii.input
cp FM*fds2ascii.csv $WDIR
cd $WDIR


