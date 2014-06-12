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

scratchdir=$SVNROOT/Utilities/Scripts/tmp
nthreads=$1
dir=$2
infile=$3

# Choose the submit and run commands

if [ "$RESOURCE_MANAGER" == "SLURM" ]; then
  QSUB="sbatch"
  QOPT="-p"
  RUNCMD="srun"
else
  QSUB="qsub"
  QOPT="-q "
  RUNCMD="$MPIDIST/bin/mpirun -np $nthreads"
fi

# If queue is "none" then use "background" to submit jobs
# instead of a resource manager (qsub or sbatch, ie a queing system).

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

if test $nthreads -le 0
then
echo "Number of threads specified is $nthreads . Must be bigger than 0."
echo "Run aborted."
exit
fi
# original method used by runmpifds.sh
#nnodes=$(echo "($nthreads-1)/16+1" | bc)
#nprocs=$(echo "($nthreads-1)/$nnodes+1" | bc)
# try running firebot using qfds method for 
# defining nodes/processes
nprocs=1
nnodes=$(echo "($nthreads-1)/$nprocs+1" | bc)
if test $nnodes -le 0
then
  nnodes=1
elif test $nnodes -gt 32
then
  nnodes=32
fi

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
if [[ $STOPFDS -gt 1 ]]; then
  echo "creating delayed stop file: $infile"
  echo $STOPFDS > $fulldir/$stopfile
  exit
elif [ $STOPFDS ]; then
  echo "stopping case: $infile"
  touch $fulldir/$stopfile
  exit
else
  if [ -e $fulldir/$stopfile ]; then
    rm $fulldir/$stopfile
  fi
fi
if [ -e $outlog ]; then
  rm $outlog
fi

# create run script

cat << EOF > $scriptfile
#!/bin/bash
#PBS -N $JOBPREFIX$infile(MPI$IB)
#PBS -l nodes=$nnodes:ppn=$nprocs
#PBS -l walltime=04:00:00
#PBS -S /bin/bash
#PBS -e $outerr
#PBS -l pvmem=1GB
#PBS -o $outlog
#\$ -N $JOBPREFIX$infile(MPI$IB)
#\$ -pe mpi $nthreads
#\$ -S /bin/bash
#\$ -e $outerr
#\$ -o $outlog
#SBATCH -J $JOBPREFIX$infile(MPI$IB)
#SBATCH --mem-per-cpu=1000
#SBATCH -n $nthreads
#SBATCH -t 04:00:00
#SBATCH -e $outerr
#SBATCH -o $outlog
#SBATCH -p $queue

cd $fulldir

echo Time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`

export LD_LIBRARY_PATH=$MPIDIST/lib:$FORTLIB:$LD_LIBRARY_PATH
$RUNCMD $FDSMPI $in 

# Run by $QSUB $QOPT $queue $scriptfile
EOF
chmod +x $scriptfile
if [ "$background" != "yes" ]; then
  echo Running `basename $FDSMPI` $in 
  chmod +x $scriptfile
  $QSUB $QOPT $queue $scriptfile
else
  echo Running `basename $FDS` $in 
  cd $fulldir
  $QSUB $FDS $in
fi

rm $scriptfile
