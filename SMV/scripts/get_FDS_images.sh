#!/bin/bash
echo copying FDS User Guide images from firebot
GUIDE=FDS-SMV/Manuals/FDS_User_Guide/SCRIPT_FIGURES
FROMDIR=~firebot/$GUIDE
TODIR=~/$GUIDE
cp $FROMDIR/*.* $TODIR/.

echo copying FDS Verification Guide images from firebot
GUIDE=FDS-SMV/Manuals/FDS_Verification_Guide/SCRIPT_FIGURES
FROMDIR=~firebot/$GUIDE
TODIR=~/$GUIDE
cp $FROMDIR/*.* $TODIR/.
cp $FROMDIR/Scatterplots/* $TODIR/Scatterplots/.

GUIDE=FDS-SMV/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES
FROMDIR=~firebot/$GUIDE
TODIR=~/$GUIDE

echo copying FDS Validation Guide images from firebot
for D in $FROMDIR/*; do
   if [ -d "${D}" ]; then
      SUBDIR=$(basename $D)
      cp $FROMDIR/$SUBDIR/* $TODIR/$SUBDIR/.
   fi
done
