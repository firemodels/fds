#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/LEMTA_Spray/Current_Results
set WDIR=$FDSSMV/Validation/LEMTA_Spray/FDS_Output_Files
cp $DDIR/LEMT*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

