#!/bin/bash

# Common_Run_All.sh
# Kristopher Overholt
# 12/4/2013

# This is a common script that is sourced by all of the individual
# Run_All.sh scripts for each validation case. To avoid code duplication,
# this script contains options and functions that are global to all of
# the Run_All. scripts, such as option flags, directory creation, and
# FDS file copying.

export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

function usage {
echo "Run_All.sh [ -d -h -o output_dir -q queue_name -s -x ]"
echo "Runs FDS validation set"
echo ""
echo "Options"
echo "-d - use debug version of FDS"
echo "-h - display this message"
echo "-o output_dir - specify output directory"
echo "     default: Current_Results"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "     other options: fire60s, fire70s, vis"
echo "-s - stop FDS runs"
echo "-x - do not copy FDS input files"
exit
}

while getopts 'dho:q:sx' OPTION
do
case $OPTION in
  d)
   DEBUG="-b"
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
      if [ "$(ls -A $INDIR)" ]; then
          while true; do
              echo "Directory $INDIR already exists with files."
              read -p "Do you wish to overwrite the files in this directory? [y]es or [n]o `echo $'\n> '`" yn
              case $yn in
                  [Yy]* ) break;; # Continue along
                  [Nn]* ) echo "Exiting."; exit;;
                  * ) echo "Please answer [y]es or [n]o. `echo $'\n '`";;
              esac
          done
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
