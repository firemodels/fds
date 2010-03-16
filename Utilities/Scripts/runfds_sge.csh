#!/bin/csh -f
set scratchdir=$SVNROOT/Utilities/Scripts/tmp
set dir=$1
set infile=$2
set tasks=$3

set fulldir=$BASEDIR/$dir

if(! -e $FDS) then
  echo "The file $FDS does not exit. Run aborted"
  exit
endif
if(! -d $fulldir) then
  echo "The directory $fulldir does not exit. Run aborted."
  exit
endif

cd $fulldir
qsub -t $tasks sge-fds-array.sh
