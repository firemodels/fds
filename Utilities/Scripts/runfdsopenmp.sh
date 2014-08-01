#!/bin/bash

# defaults

queue=
background=no
QSUB=qsub
nnodes=1
nprocs=8
nthreads=1
benchmark=no

while getopts 'bn:o:q:' OPTION
do
case $OPTION in
  b)
   benchmark=yes
   ;;
  n)
   nprocs="$OPTARG"
   nthreads="$OPTARG"
   ;;
  o)
   nthreads="$OPTARG"
   ;;
  q)
   queue="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

# in benchmark mode run a case "alone" on a node
if [ "$benchmark" == "yes" ]; then
  nprocs=8
fi

if [ "$JOBPREFIX" == "" ]; then
  JOBPREFIX=VV_
fi

# If queue is "none" then use "background" to submit jobs
# instead of qsub (ie a queing system).

if [ "$queue" == "none" ]; then
  queue=
  QSUB="$BACKGROUND -u 75 -d 10 "
  background=yes;
  notfound=`$BACKGROUND -help 2>&1 | tail -1 | grep "not found" | wc -l`
  if [ "$notfound" == "1" ];  then
    echo "The program $BACKGROUND is not available. Run aborted"
    exit
  fi
fi
if [ "$queue" != "" ]; then
   queue="-q $queue"
fi

# setup parameters

scratchdir=$SVNROOT/Utilities/Scripts/tmp
dir=$1
infile=$2

fulldir=$BASEDIR/$dir
in=$infile.fds
outerr=$fulldir/$infile.err
outlog=$fulldir/$infile.log
stopfile=$infile.stop

scriptfile=$scratchdir/script.$$

# ensure that various files and directories exist

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

if [ -e $fulldir/$stopfile ]; then
  rm $fulldir/$stopfile
fi
if [ $STOPFDS ]; then
  echo "stopping case: $infile"
  touch $fulldir/$stopfile
  exit
fi
if [ $STOPFDSMAXITER ]; then
  echo "creating delayed stop file: $infile"
  echo $STOPFDSMAXITER > $fulldir/$stopfile
fi

if [ -e $outlog ]; then
  rm $outlog
fi

# create run script

cat << EOF > $scriptfile
#!/bin/bash
#\$ -S /bin/bash
#\$ -N $JOBPREFIX$infile -e $outerr -o $outlog
#PBS -l nodes=$nnodes:ppn=$nprocs
#PBS -N $JOBPREFIX$infile -e $outerr -o $outlog
cd $fulldir

echo Time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`

# Set number of OpenMP threads on target machine
export OMP_NUM_THREADS=$nthreads

$FDS $in 
EOF

echo Running `basename $FDS` $in 
if [ "$background" != "yes" ]; then
  chmod +x $scriptfile
  $QSUB $queue $scriptfile
else
  cd $fulldir
  $QSUB $FDS $in
fi

rm $scriptfile
