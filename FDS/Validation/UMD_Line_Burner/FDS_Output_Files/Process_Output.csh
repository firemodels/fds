#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/UMD_Line_Burner/Current_Results
set WDIR=$FDSSMV/Validation/UMD_Line_Burner/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

