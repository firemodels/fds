#!/bin/csh -f
set scratchdir=$SVNROOT/Utilities/Scripts/tmp
set dir=$1
set infile=$2
set host=bluesky.cfr.nist.gov

set fulldir=$BASEDIR/$dir
set in=$infile.fds
set out=$infile.err
set stopfile=$infile.stop

if(! -e $FDS) then
  echo "The file $FDS does not exit. Run aborted"
  exit
endif
if(! -d $fulldir) then
  echo "The directory $fulldir does not exit. Run aborted."
  exit
endif
if(! -e $fulldir/$in) then
  echo "The fds input  files, $fulldir/$in, does not exit. Run aborted."
  exit
endif
if($?STOPFDS) then
 echo "stopping case: $infile"
 touch $fulldir/$stopfile
 exit
endif
if(-e $fulldir/$stopfile) then
 rm $fulldir/$stopfile
endif
cd $fulldir
echo $in
$FDS $in >& $out
