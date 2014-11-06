#!/bin/bash

CURDIR=`pwd`
HOST=blaze.el.nist.gov
FDSREPO=FDS-SMV
DIRS="FDS_Technical_Reference_Guide FDS_User_Guide FDS_Verification_Guide"

for dir in $DIRS
do
  scp $HOST\:~firebot/$FDSREPO/Manuals/$dir/SCRIPT_FIGURES/* ~/$FDSREPO/Manuals/$dir/SCRIPT_FIGURES/.
done

dir=FDS_Verification_Guide/
scp $HOST\:~firebot/$FDSREPO/Manuals/$dir/SCRIPT_FIGURES/Scatterplots/* ~/$FDSREPO/Manuals/$dir/SCRIPT_FIGURES/Scatterplots/.

cd ../FDS_Validation_Guide/SCRIPT_FIGURES
for dir in *
do
  scp $HOST\:~firebot/$FDSREPO/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/$dir/* ~/$FDSREPO/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/$dir/.
done
cd $CURDIR
