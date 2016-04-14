#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Smyth_Slot_Burner/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Smyth_Slot_Burner/Current_Results
cp $DDIR/Smy*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

