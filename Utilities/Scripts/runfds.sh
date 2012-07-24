#!/bin/bash -f

# if a queue is not specified with the -q option then
# use the system's default queue

queue=
while getopts 'q:' OPTION
do
case $OPTION in
  q)
   queue="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

[ "$queue" != "" ] ; queue="-q $queue"

scratchdir=$SVNROOT/Utilities/Scripts/tmp
dir=$1
infile=$2

fulldir=$BASEDIR/$dir
in=$infile.fds
outerr=$fulldir/$infile.err
outlog=$fulldir/$infile.log
stopfile=$infile.stop

scriptfile=$scratchdir/script.$$
if ! [ -e $FDS ];  then
  echo "The file $FDS does not exist. Run aborted"
  exit
fi
if ! [ -d $fulldir ]; then
  echo "The directory $fulldir does not exist. Run aborted."
  exit
fi
if ! [ -e $fulldir/$in ]; then
  echo "The fds input file, $fulldir/$in, does not exist. Run aborted."
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
if [ -e $outlog ]; then
 rm $outlog
fi
cat << EOF > $scriptfile
#!/bin/bash -f
#\$ -S /bin/bash
#\$ -N VV_$infile -e $outerr -o $outlog
#PBS -N VV_$infile -e $outerr -o $outlog
cd $fulldir

echo Time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`

$FDS $in 
EOF
chmod +x $scriptfile
echo Running $in 
qsub $queue $scriptfile
rm $scriptfile
