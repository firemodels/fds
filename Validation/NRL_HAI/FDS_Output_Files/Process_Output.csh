#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/NRL_HAI/FDS_Output_Files
set DDIR=$SVNROOT/Validation/NRL_HAI/Current_Results
cd $DDIR
$PDIR/fds2ascii < $WDIR/NRL_HAI_1_fds2ascii.input
$PDIR/fds2ascii < $WDIR/NRL_HAI_2_fds2ascii.input
$PDIR/fds2ascii < $WDIR/NRL_HAI_3_fds2ascii.input
$PDIR/fds2ascii < $WDIR/NRL_HAI_4_fds2ascii.input
$PDIR/fds2ascii < $WDIR/NRL_HAI_5_fds2ascii.input
$PDIR/fds2ascii < $WDIR/NRL_HAI_6_fds2ascii.input
$PDIR/fds2ascii < $WDIR/NRL_HAI_7_fds2ascii.input
$PDIR/fds2ascii < $WDIR/NRL_HAI_8_fds2ascii.input
$PDIR/fds2ascii < $WDIR/NRL_HAI_9_fds2ascii.input
cp NRL*fds2ascii.csv $WDIR
cd $WDIR


