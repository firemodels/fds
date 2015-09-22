#!/bin/bash

export FIREBOTROOT=/home/firebot/SMOKEVIS2/FDS-SMVgitclean
export FIREBOTMANS=$FIREBOTROOT/Manuals
export FBUG=$FIREBOTMANS/FDS_User_Guide
export FBVG=$FIREBOTMANS/FDS_Verification_Guide
export FBVAL=$FIREBOTMANS/FDS_Validation_Guide
export BASEDIR=`pwd`

# Copy User's Guide Figures
cp $FBUG/SCRIPT_FIGURES/* $BASEDIR/FDS_User_Guide/SCRIPT_FIGURES/
echo Users Guide Figures Copied

# Copy Verification Guide Figures
cp $FBVG/SCRIPT_FIGURES/*.pdf $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/.
cp $FBVG/SCRIPT_FIGURES/*.png $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/.
cp $FBVG/SCRIPT_FIGURES/Scatterplots/*.tex $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/.
echo Verification Figures Copied

# Copy Validation Guide Figures
#cp -R $FBVAL/SCRIPT_FIGURES/* $BASEDIR/FDS_Validation_Guide/SCRIPT_FIGURES/ &> /dev/null
rsync -r --exclude=*.git $FBVAL/SCRIPT_FIGURES/* $BASEDIR/FDS_Validation_Guide/SCRIPT_FIGURES/
echo Validation Guide Figures Copied

