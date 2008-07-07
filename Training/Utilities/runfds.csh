#!/bin/csh -f
set bindir=~/bin
set dir=$1
set infile=$2
set host=$3

set fulldir=$JOBDIR/$dir
set in=$infile.fds
set out=$infile.err
set stopfile=$infile.stop

set scriptfile=$bindir/script.$$
if(! -e $FDS5) then
  echo "The file $FDS5 does not exit. Run aborted"
  exit
endif
if(! -d $fulldir) then
  echo "The directory $fulldir does not exit"
  exit
endif
if(! -e $fulldir/$in) then
  echo "The fds input  file, $fulldir/$in does not exit"
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
cat << EOF > $scriptfile
#!/bin/csh -f
cd $fulldir
# note: The environment variable, FDS5, is defined 
#       in the calling script, runall.csh
$FDS5 $in >& $out
EOF
chmod +x $scriptfile
echo Running $in on $host
ssh -n $host $scriptfile
rm $scriptfile
