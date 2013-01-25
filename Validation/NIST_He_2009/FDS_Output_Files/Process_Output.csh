#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/NIST_He_2009/FDS_Output_Files
set DDIR=$SVNROOT/Validation/NIST_He_2009/Test
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

