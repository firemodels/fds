#!/bin/csh -f
set bindir=~/bin
set dir=$1
set infile=$2
set host=$3
set user=`whoami`

set fulldir=$BASEDIR/$dir
set in=$infile.fds
set out=$infile.err
set stopfile=$infile.stop

set scriptfile=$bindir/script.$$
if(! -e $FDS) then
  echo "The file $FDS does not exist. Run aborted"
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
$FDS $in >& $out
EOF
set startdate=`date`
set pid=$$
chmod +x $scriptfile
cat << EOF  | Mail -s "FDS run status $pid STARTED $in" $user
      host: $host
 directory: $fulldir
start time: $startdate
 stop time:
     input: $in
EOF
echo Running $in on $host
ssh -n $host $scriptfile
rm $scriptfile
set stopdate=`date`
cat << EOF  | Mail -s "FDS run status $pid FINISHED $in" $user
      host: $host
 directory: $fulldir
start time: $startdate
 stop time: $stopdate
     input: $in
EOF
rm $scriptfile
