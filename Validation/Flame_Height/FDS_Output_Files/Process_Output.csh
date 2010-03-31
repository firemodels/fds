#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/Flame_Height/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Flame_Height/Current_Results
cd $DDIR
$PDIR/fds2ascii $WDIR/Qs=10000_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=1000_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=100_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=10_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=1_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=2000_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=200_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=20_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=2_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=5000_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=500_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=50_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=5_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=p1_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=p2_RI=05_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=p5_RI=05_fds2ascii.input

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

$PDIR/fds2ascii $WDIR/Qs=10000_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=1000_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=100_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=10_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=1_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=2000_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=200_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=20_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=2_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=5000_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=500_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=50_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=5_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=p1_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=p2_RI=20_fds2ascii.input
$PDIR/fds2ascii $WDIR/Qs=p5_RI=20_fds2ascii.input
cp Qs*fds2ascii.csv $WDIR
cd $WDIR
$PDIR/flame_height > $WDIR/FDS_Flame_Height.csv
cd ./

