#!/bin/bash

# Common_Run_All.sh
# Kristopher Overholt
# 12/4/2013

# This is a common script that is sourced by all of the individual
# Run_All.sh scripts for each validation case. To avoid code duplication,
# this script contains options and functions that are global to all of
# the Run_All. scripts, such as option flags, directory creation, and
# FDS file copying.

CURDIR=`pwd`

cd $SVNROOT/Utilities/Scripts/
SCRIPTDIR=`pwd`
export QFDS=$SCRIPTDIR/qfds.sh
cd $CURDIR

export BASEDIR=`pwd`
export INDIR=Current_Results

function usage {
echo "Run_All.sh [ -b -h -o output_dir -q queue_name -s -x ]"
echo "Runs FDS validation set"
echo ""
echo "Options"
echo "-b - use debug version of FDS"
echo "-h - display this message"
echo "-o output_dir - specify output directory"
echo "     default: Current_Results"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "     other options: fire60s, fire70s, vis"
echo "-s - stop FDS runs"
echo "-x - do not copy FDS input files"
echo "-y - overwrite existing files"
exit
}

DEBUG=$OPENMP
while getopts 'bho:q:sxy' OPTION
do
case $OPTION in
  b)
   DEBUG="-b $OPENMP"
   ;;
  h)
  usage;
   ;;
  o)
   INDIR="$OPTARG"
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
  s)
   export STOPFDS=1
   ;;
  x)
   export DONOTCOPY=1
   ;;   
  y)
   export OVERWRITE=1
   ;;   
esac
done

if [ "$QUEUE" != "" ]; then
   QUEUE="-q $QUEUE"
fi

##############################################################

# Skip if STOPFDS (-s option) is specified
if [ ! $STOPFDS ] ; then
  # Check for existence of $INDIR (Current_Results) directory
  if [ -d "$INDIR" ]; then
      # Check for files in $INDIR (Current_Results) directory
      if [[ "$(ls -A $INDIR)" && ! $OVERWRITE ]]; then
          echo "Directory $INDIR already exists with files."
          echo "Use the -y option to overwrite files."
          echo "Exiting."
          exit
      elif [[ "$(ls -A $INDIR)" && $OVERWRITE ]]; then
        # Continue along
        :
      fi
  # Create $INDIR (Current_Results) directory if it doesn't exist
  else
     mkdir $INDIR
  fi
fi

if [ ! $DONOTCOPY ] ; then
  # Copy FDS input files to $INDIR (Current_Results) directory
  cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR
fi
