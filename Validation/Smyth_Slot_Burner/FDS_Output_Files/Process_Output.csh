#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/Smyth_Slot_Burner/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Smyth_Slot_Burner/Current_Results
cd $DDIR
$PDIR/fds2ascii < $WDIR/Smyth_Slot_Burner_7mm_fds2ascii.inp
$PDIR/fds2ascii < $WDIR/Smyth_Slot_Burner_9mm_fds2ascii.inp
$PDIR/fds2ascii < $WDIR/Smyth_Slot_Burner_11mm_fds2ascii.inp
cp Smy*fds2ascii.csv $WDIR

