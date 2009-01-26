#!/bin/csh -f
set PDIR=~/FDS_Repository/Utilities/Data_Processing
set WDIR=~/FDS_Repository/Validation/Ulster_SBI/FDS_Output_Files
set DDIR=~/VALIDATION/Ulster_SBI/FDS_5.3
cd $DDIR
$PDIR/Ulster_SBI_fds2ascii < $WDIR/Ulster_SBI_30_kW_fds2ascii.input
$PDIR/Ulster_SBI_fds2ascii < $WDIR/Ulster_SBI_45_kW_fds2ascii.input
$PDIR/Ulster_SBI_fds2ascii < $WDIR/Ulster_SBI_60_kW_fds2ascii.input
cp Ulster*fds2ascii.csv $WDIR
cd $WDIR

