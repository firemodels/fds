#!/bin/bash -f
scratchdir=$SVNROOT/Utilities/Scripts/tmp
dir=$1
infile=$2

fulldir=$BASEDIR/$dir
in=$infile.fds
out=$infile.err
stopfile=$infile.stop

scriptfile=$scratchdir/script.$$
if ! [ -e $FDS ];  then
  echo "The file $FDS does not exit. Run aborted"
  exit
fi
if ! [ -d $fulldir ]; then
  echo "The directory $fulldir does not exit. Run aborted."
  exit
fi
if ! [ -e $fulldir/$in ]; then
  echo "The fds input file, $fulldir/$in, does not exit. Run aborted."
  exit
fi
if [ $STOPFDS ]; then
 echo "stopping case: $infile"
 touch $fulldir/$stopfile
 exit
fi
if [ -e $fulldir/$stopfile ]; then
 rm $fulldir/$stopfile
fi
cat << EOF > $scriptfile
#!/bin/bash -f
cd $fulldir
#$ -N Ver_'$infile'
#$ -wd $fulldir
#$ -e $out -o /dev/null 
$FDS $in 
EOF
chmod +x $scriptfile
echo Running $in 
qsub $scriptfile
rm $scriptfile
