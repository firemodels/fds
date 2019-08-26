#!/bin/bash

export FIREBOTROOT=/home2/smokevis2/firebot/FireModels_clone/fds/ # blaze firebot
#FIREBOTROOT=/home4/firebot/FireModels_clone/fds/ # burn firebot
export FIREBOTMANS=$FIREBOTROOT/Manuals/
export FIREBOTVER=$FIREBOTROOT/Verification/
export FIREBOTVAL=$FIREBOTROOT/Validation/
export FBTG=$FIREBOTMANS/FDS_Technical_Reference_Guide/
export FBUG=$FIREBOTMANS/FDS_User_Guide/
export FBVG=$FIREBOTMANS/FDS_Verification_Guide/
export FBVAL=$FIREBOTMANS/FDS_Validation_Guide/
export BASEDIR=`pwd`

# Copy Tech Guide Figures
cp $FBTG/SCRIPT_FIGURES/* $BASEDIR/FDS_Technical_Reference_Guide/SCRIPT_FIGURES/
echo Tech Guide Figures Copied

# Copy User's Guide Figures
cp $FBUG/SCRIPT_FIGURES/* $BASEDIR/FDS_User_Guide/SCRIPT_FIGURES/
echo Users Guide Figures Copied

# Copy Verification Guide Figures
cp $FBVG/SCRIPT_FIGURES/*.pdf $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/.
cp $FBVG/SCRIPT_FIGURES/*.png $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/.
cp $FBVG/SCRIPT_FIGURES/*.tex $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/.
cp $FBVG/SCRIPT_FIGURES/Scatterplots/*.tex $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/.
echo Verification Figures Copied

# Copy Validation Guide Figures
#cp -R $FBVAL/SCRIPT_FIGURES/* $BASEDIR/FDS_Validation_Guide/SCRIPT_FIGURES/ &> /dev/null
rsync -r --exclude=*.git $FBVAL/SCRIPT_FIGURES/* $BASEDIR/FDS_Validation_Guide/SCRIPT_FIGURES/
echo Validation Guide Figures Copied

# Copy Verification Results
#rsync -v -r --include '*/' --include '*_git.txt' --include '*.csv' --include '*.prt5' --exclude '*' $FIREBOTVER/* $BASEDIR/../Verification/
#cp $FIREBOTVER/Miscellaneous/mesh_transformation.smv $BASEDIR/../Verification/Miscellaneous/.
#echo Verification Results Copied

# Copy Validation Results
#rsync -v -r --include '*/' --include '*_git.txt' --include '*.csv' --exclude '*' $FIREBOTVAL/* $BASEDIR/../Validation/
#echo Validation Results Copied

