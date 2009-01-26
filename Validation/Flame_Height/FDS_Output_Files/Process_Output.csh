#!/bin/csh -f
set PDIR=~/FDS_Repository/Utilities/Data_Processing
set WDIR=~/FDS_Repository/Validation/Flame_Height/FDS_Output_Files
set DDIR=~/VALIDATION/FLAME_HEIGHT/FDS_5.3
cd $DDIR
$PDIR/fds2ascii $WDIR/Qs=10000_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=1000_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=100_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=10_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=1_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=2000_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=200_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=2_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=5000_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=500_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=50_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=5_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=p1_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=p2_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=p5_fds2ascii.input
cp Qs*fds2ascii.csv $WDIR
cd $WDIR
$PDIR/flame_height > $WDIR/FDS_Flame_Height.csv
cd ./

