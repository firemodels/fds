#!/bin/bash

FDSREPO=FDS-SMV

while getopts 'r:' OPTION
do
case $OPTION  in
  r)
   FDSREPO="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

echo copying FDS User Guide images from firebot
FIGUREDIR=$FDSREPO/Manuals/FDS_User_Guide/SCRIPT_FIGURES
cp ~firebot/$FIGUREDIR/*.* ~/$FIGUREDIR/.

echo copying FDS Verification Guide images from firebot
FIGUREDIR=$FDSREPO/Manuals/FDS_Verification_Guide/SCRIPT_FIGURES
cp ~firebot/$FIGUREDIR/*.* ~/$FIGUREDIR/.
cp ~firebot/$FIGUREDIR//Scatterplots/* ~/$FIGUREDIR/Scatterplots/.

FIGUREDIR=$FDSREPO/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES
echo copying FDS Validation Guide images from firebot
for D in $FROMDIR/*; do
   if [ -d "${D}" ]; then
      SUBDIR=$(basename $D)
      cp ~firebot/$FIGUREDIR/$SUBDIR/* ~/$FIGUREDIR/$SUBDIR/.
   fi
done
