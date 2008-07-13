#!/bin/csh -f
set bindir=~/bin
set dir=$1
set infile=$2

set fulldir=$JOBDIR/$dir
set in=$infile
set insmv=$infile.smv

if(! -e $SMV5) then
  echo "The file $SMV5 does not exit. Run aborted"
  exit
endif
if(! -d $fulldir) then
  echo "The directory $fulldir does not exit. Run aborted."
  exit
endif
if(! -e $fulldir/$insmv) then
  echo "The Smokeview file, $fulldir/$insmv does not exit. Run aborted."
  exit
endif
cd $fulldir
$SMV5 -runscript $in
if(! -d $REPORTFIGDIR) then
  echo "*** Warning, the directory $REPORTFIGDIR does not exist - "
  echo "    Smokeview generated figures not copied."
endif
mv *.png $REPORTFIGDIR/.
