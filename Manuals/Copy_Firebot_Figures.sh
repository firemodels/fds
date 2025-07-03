#!/bin/bash

USE_SSH=

# uncomment following USE_SSH line to use ssh to copy figures
# use when you can't cross mount directories containing firebot
#USE_SSH=1

# set firebot host to spark
# if using USE_SSH=1 over VPN add user name to HOST, e.g., HOST=usrname@spark.nist.gov
HOST=spark.nist.gov
HOSTDIR=/home/firebot/FireModels_clone/fds

# shouldn't have to change lines below

if [ "$USE_SSH" == "" ]; then
  CP=cp
  export FIREBOTROOT=$HOSTDIR
fi

if [ "$USE_SSH" == "1" ]; then
  CP="scp -q"
  export FIREBOTROOT=$HOST:$HOSTDIR
fi

export FIREBOTMANS=$FIREBOTROOT/Manuals/
export FIREBOTVER=$FIREBOTROOT/Verification/
export FIREBOTVAL=$FIREBOTROOT/Validation/
export FBTG=$FIREBOTMANS/FDS_Technical_Reference_Guide/
export FBUG=$FIREBOTMANS/FDS_User_Guide/
export FBVG=$FIREBOTMANS/FDS_Verification_Guide/
export FBVAL=$FIREBOTMANS/FDS_Validation_Guide/
export BASEDIR=`pwd`

# Copy Tech Guide Figures
$CP $FBTG/SCRIPT_FIGURES/* $BASEDIR/FDS_Technical_Reference_Guide/SCRIPT_FIGURES/
echo Tech Guide Figures Copied

# Copy User's Guide Figures
$CP $FBUG/SCRIPT_FIGURES/* $BASEDIR/FDS_User_Guide/SCRIPT_FIGURES/
echo Users Guide Figures Copied

# Copy Verification Guide Figures
$CP $FBVG/SCRIPT_FIGURES/*.pdf $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/.
$CP $FBVG/SCRIPT_FIGURES/*.png $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/.
$CP $FBVG/SCRIPT_FIGURES/*.tex $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/.
$CP $FBVG/SCRIPT_FIGURES/Scatterplots/*.tex $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/.
echo Verification Figures Copied

# Copy Validation Guide Figures
#cp -R $FBVAL/SCRIPT_FIGURES/* $BASEDIR/FDS_Validation_Guide/SCRIPT_FIGURES/ &> /dev/null
rsync -r --exclude=*.git $FBVAL/SCRIPT_FIGURES/* $BASEDIR/FDS_Validation_Guide/SCRIPT_FIGURES/
echo Validation Guide Figures Copied

# Copy Verification Results
#rsync -v -r --include '*/' --include '*_git.txt' --include '*.csv' --include '*.err' --exclude '*' $FIREBOTVER/* $BASEDIR/../Verification/
#$CP $FIREBOTVER/Miscellaneous/mesh_transformation.smv $BASEDIR/../Verification/Miscellaneous/.
#$CP $FIREBOTVER/Sprinklers_and_Sprays/fluid_part_mom_x_1.prt5 $BASEDIR/../Verification/Sprinklers_and_Sprays/.
#$CP $FIREBOTVER/Sprinklers_and_Sprays/fluid_part_mom_y_1.prt5 $BASEDIR/../Verification/Sprinklers_and_Sprays/.
#$CP $FIREBOTVER/Sprinklers_and_Sprays/fluid_part_mom_z_1.prt5 $BASEDIR/../Verification/Sprinklers_and_Sprays/.
#$CP $FIREBOTVER/Sprinklers_and_Sprays/terminal_velocity_dt_0_0001_1.prt5 $BASEDIR/../Verification/Sprinklers_and_Sprays/.
#$CP $FIREBOTVER/Sprinklers_and_Sprays/terminal_velocity_dt_0_001_1.prt5 $BASEDIR/../Verification/Sprinklers_and_Sprays/.
#$CP $FIREBOTVER/Sprinklers_and_Sprays/terminal_velocity_dt_0_01_1.prt5 $BASEDIR/../Verification/Sprinklers_and_Sprays/.
#$CP $FIREBOTVER/Sprinklers_and_Sprays/terminal_velocity_dt_0_1_1.prt5 $BASEDIR/../Verification/Sprinklers_and_Sprays/.
#$CP $FIREBOTVER/Sprinklers_and_Sprays/terminal_velocity_dt_1_0_1.prt5 $BASEDIR/../Verification/Sprinklers_and_Sprays/.
#echo Verification Results Copied

# Copy Validation Results
#rsync -v -r --include '*/' --include '*_git.txt' --include '*.csv' --exclude '*' $FIREBOTVAL/* $BASEDIR/../Validation/
#echo Validation Results Copied

