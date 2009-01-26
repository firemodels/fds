#!/bin/csh -f
set PDIR=~/FDS_Repository/Utilities/Data_Processing
set WDIR=~/FDS_Repository/Validation/McCaffrey_Plume/FDS_Output_Files
set DDIR=~/VALIDATION/MCCAFFREY_PLUME/FDS_5.3
cd $DDIR
$PDIR/fds2ascii $WDIR/McCaffrey_14_kW_fds2ascii.input
$PDIR/fds2ascii $WDIR/McCaffrey_22_kW_fds2ascii.input
$PDIR/fds2ascii $WDIR/McCaffrey_33_kW_fds2ascii.input
$PDIR/fds2ascii $WDIR/McCaffrey_45_kW_fds2ascii.input
$PDIR/fds2ascii $WDIR/McCaffrey_57_kW_fds2ascii.input
cp McC*fds2ascii.csv $WDIR

