#!/bin/bash

# This is a common script that is sourced by all of the individual
# Run_All.sh scripts for each validation case. To avoid code duplication,
# this script contains options and functions that are global to all of
# the Run_All. scripts, such as option flags, directory creation, and
# FDS file copying.

CURDIR=`pwd`

cd $SVNROOT/..

cd $SVNROOT/Utilities/Scripts/
SCRIPTDIR=`pwd`
cd $CURDIR

export BASEDIR=`pwd`
export INDIR=Current_Results
JOB_PREFIX=
export STOPFDSMAXITER=
DV=
TCP=
EXE=
resource_manager=
walltime=
showcommandline=
showscript=
CHECK=

INTEL="-I"
# the mac doesn't have Intel MPI
if [ "`uname`" == "Darwin" ] ; then
  INTEL=
fi

function get_full_path {
  filepath=$1

  if [[ $filepath == /* ]]; then
    full_filepath=$filepath
  else
    dir_filepath=$(dirname  "${filepath}")
    filename_filepath=$(basename  "${filepath}")
    curdir=`pwd`
    cd $dir_filepath
    full_filepath=`pwd`/$filename_filepath
    cd $curdir
  fi
}

function usage {
echo "Run_All.sh [ -b -h -o output_dir -q queue_name -s -x ]"
echo "Runs FDS validation set"
echo ""
echo "Options"
echo "-b - use debug version of FDS"
echo "-C - check that case has run"
echo "-e exe - run using exe (full path to fds)."
echo "      Note: environment must be defined to use this executable"
echo "-h - display this message"
echo "-I - run with Intel MPI version of fds"
echo "-j job_prefix - specify job prefix"
echo "-m n - run cases only n time steps"
echo "-o output_dir - specify output directory"
echo "     default: Current_Results"
echo "-O - run with Open MPI version of fds"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "-r resource_manager - default: SLURM, other options: PBS/Torque"
echo "-s - stop FDS runs"
echo "-u - use development version of FDS"
echo "-v - show script run by qfds.sh for each validation case"
echo "-V - show qfds.sh command line for each validation case"
echo "-w walltime - default: empty, PBS: hh:mm:ss, SLURM: dd-hh:mm:ss"
echo "-x - do not copy FDS input files"
echo "-y - remove existing files"
echo "-z - overwrite existing files"
exit
}

DEBUG=
while getopts 'bCe:EhIj:m:o:Oq:r:suvVw:xyz' OPTION
do
case $OPTION in
  b)
   DEBUG="-T db "
   ;;
  C)
   CHECK=1
   ;;
  e)
   EXE="$OPTARG"
   INTEL=
   DV=
   DEBUG=
   ;;
  E)
   TCP="-E "
   ;;
  h)
  usage;
   ;;
  j)
   JOBPREFIX="-j $OPTARG"
   ;;
  I)
   INTEL="-I"
   EXE=
   ;;
  m)
   export STOPFDSMAXITER="$OPTARG"
   ;;
  o)
   INDIR="$OPTARG"
   ;;
  O)
   INTEL="-L"
   EXE=
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
  r)
   resource_manager="$OPTARG"
   ;;
  s)
   export STOPFDS=1
   ;;
  u)
  DV="-T dv"
   ;;
  V)
  showcommandline="-V"
   ;;
  v)
  showscript="-v"
   ;;
  w)
   walltime="-w $OPTARG"
   ;;
  x)
   export DONOTCOPY=1
   ;;   
  y)
   export DELETEFILES=1
   export OVERWRITE=
   ;;   
  z)
   export OVERWRITE=1
   export DELETEFILES=
   ;;   
esac
done

if [ "$EXE" != "" ]; then
  get_full_path $EXE
  EXE="-e $full_filepath"
fi

export QFDS="$SCRIPTDIR/qfds.sh $walltime $showcommandline $showscript $DV $INTEL $EXE"
if [ "$CHECK" != "" ]; then
  export QFDS="$SVNROOT/fds/Verification/scripts/Check_FDS_Cases.sh"
fi

if [ "$QUEUE" != "" ]; then
   QUEUE="-q $QUEUE"
fi
DEBUG="$DEBUG $JOBPREFIX"
DEBUG="$DEBUG $TCP"

if [ "$resource_manager" == "SLURM" ]; then
   export RESOURCE_MANAGER="SLURM"
else
   export RESOURCE_MANAGER="PBS/Torque"
fi
##############################################################

# Skip if STOPFDS (-s option) is specified
if [ ! $STOPFDS ] ; then
  # Check for existence of $INDIR (Current_Results) directory
  if [ -d "$INDIR" ]; then
      # Check for files in $INDIR (Current_Results) directory
      if [[ "$(ls -A $INDIR)" && "$OVERWRITE" == "" && "$DELETEFILES" == "" ]]; then
          echo "The directory $INDIR already exists with files."
          echo "Use the -y option to REMOVE existing files or -z to OVERWRITE existing files"
          echo "Exiting."
          exit
      elif [[ "$(ls -A $INDIR)" && "$DELETEFILES" == "1" ]]; then
          rm $INDIR/*
      fi
  # Create $INDIR (Current_Results) directory if it doesn't exist
  else
     mkdir $INDIR
  fi
fi

if [ ! $DONOTCOPY ] ; then
  # Copy FDS input files to $INDIR (Current_Results) directory
  cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR
  HAVE_TXT=`git ls-files $BASEDIR/FDS_Input_Files/*.txt`
  if [ "$HAVE_TXT" != "" ]; then
    cp $BASEDIR/FDS_Input_Files/*.txt $BASEDIR/$INDIR
  fi
fi
