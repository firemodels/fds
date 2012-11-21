#!/bin/csh -f
set bindir=~/bin
set dir=$1
set infile=$2
set host=$3

set fulldir=$BASEDIR/$dir
set in=$infile.fds
set out=$infile.err
set stopfile=$infile.stop

set scriptfile=$bindir/script.$$
if(! -e $FDSMPI) then
  echo "The file $FDSMPI does not exist. Run aborted"
  exit
endif
if($?LAMNODES) then
else
  echo "The environment variable  LAMNODES is not defined. Run aborted."
  exit
endif
if(! -d $fulldir) then
  echo "The directory $fulldir does not exist. Run aborted."
  exit
endif
if(! -e $fulldir/$in) then
  echo "The fds input  file, $fulldir/$in does not exist. Run aborted."
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
mpirun $LAMNODES $FDSMPI $in >& $out
EOF
chmod +x $scriptfile
echo Running $in on $host
ssh -n $host $scriptfile
rm $scriptfile
