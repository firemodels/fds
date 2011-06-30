#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/Smyth_Slot_Burner/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Smyth_Slot_Burner/Current_Results
cp $DDIR/Smy*line.csv $WDIR
cp $DDIR/*svn.txt $WDIR

