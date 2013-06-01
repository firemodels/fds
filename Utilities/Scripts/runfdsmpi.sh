#!/bin/bash

queue=
background=no
QSUB=qsub

function usage {
  echo "Usage: runfdsmpi.sh nthreads dir casename"
  echo ""
  echo "Runs a parallel FDS case on a Linux cluster using the "
  echo "PBS/SGE qsub batch queuing command"
  echo ""
  echo "nthreads - number of threads (usually number of &mesh lines)"
  echo "     dir - directory containing FDS case"
  echo "casename - fds case (without .fds extension)"
  echo
  exit
}
IB=
if [ "$FDSNETWORK" == "infiniband" ] ; then
IB=_IB
fi
if [ "$JOBPREFIX" == "" ]; then
  JOBPREFIX=VV_
fi

while getopts 'hq:' OPTION
do
case $OPTION in
  h)
  usage;
  ;;
  q)
   queue="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

# If queue is "none" then use "background" to submit jobs
# instead of qsub (ie a queing system).

if [ "$queue" == "none" ]; then
  queue=
  QSUB="$BACKGROUND -u 75 -d 10 "
  background=yes;
  if ! [ -e $BACKGROUND ];  then
    echo "The file $BACKGROUND does not exist. Run aborted"
    exit
  fi
fi
if [ "$queue" != "" ]; then
   queue="-q $queue"
fi

scratchdir=$SVNROOT/Utilities/Scripts/tmp
nthreads=$1
dir=$2
infile=$3

if test $nthreads -le 0
then
echo "Number of threads specified is $nthreads . Must be bigger than 0."
echo "Run aborted."
exit
fi
nnodes=$(echo "($nthreads-1)/8+1" | bc)
nprocs=$(echo "($nthreads-1)/$nnodes+1" | bc)

fulldir=$BASEDIR/$dir
in=$infile.fds
outerr=$fulldir/$infile.err
outlog=$fulldir/$infile.log
stopfile=$infile.stop

scriptfile=$scratchdir/script.$$
if ! [ -e $FDSMPI ];  then
  echo "The file $FDSMPI does not exist. Run aborted"
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
#!/bin/bash
#PBS -N $JOBPREFIX$infile(MPI$IB)
#PBS -l nodes=$nnodes:ppn=$nprocs
#PBS -S /bin/bash
#PBS -e $outerr
#PBS -o $outlog
#\$ -N $JOBPREFIX$infile(MPI$IB)
#\$ -pe mpi $nthreads
#\$ -S /bin/bash
#\$ -e $outerr
#\$ -o $outlog

cd $fulldir

echo Time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`

export LD_LIBRARY_PATH=$MPIDIST/lib:$FORTLIB:$LD_LIBRARY_PATH
$MPIDIST/bin/mpirun -np $nthreads $FDSMPI $in 
EOF
chmod +x $scriptfile
if [ "$background" != "yes" ]; then
  echo Running `basename $FDSMPI` $in 
  chmod +x $scriptfile
  $QSUB $queue $scriptfile
else
  echo Running `basename $FDS` $in 
  cd $fulldir
  $QSUB $FDS $in
fi

rm $scriptfile
