#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/UL_NIST_Vents/FDS_Output_Files
set DDIR=$SVNROOT/Validation/UL_NIST_Vents/Current_Results
cd $DDIR
#$PDIR/layer_height < $WDIR/ATF_Corridors_050_kW_HGL_1.input
#cp UL*HGL.csv $WDIR
cp UL*devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR


