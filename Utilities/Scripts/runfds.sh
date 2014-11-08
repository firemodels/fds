#!/bin/bash

# defaults

queue=
background=no
QSUB=qsub
USEFDS=yes
nthreads=1

if [ "$JOBPREFIX" == "" ]; then
  JOBPREFIX=VV_
fi

while getopts 'n:o:q:w' OPTION
do
case $OPTION in
  n)
  nthreads="$OPTARG"
  ;;
  o)
  ignored="$OPTARG"
  ;;
  q)
   queue="$OPTARG"
   ;;
  w)
   FDS=$WFDS
   nthreads=1
   USEFDS=no
   ;;
esac
done
shift $(($OPTIND-1))

scratchdir=$SVNROOT/Utilities/Scripts/tmp
dir=$1
infile=$2
debug_flag=$3

# If queue is "none" then use "background" to submit jobs
# instead of qsub (ie a queing system).

# Choose the submit and run commands

if [ "$RESOURCE_MANAGER" == "SLURM" ]; then
  QSUB="sbatch"
  QOPT="-p"
  RUNCMD="srun"
else
  QSUB="qsub"
  QOPT="-q"
  RUNCMD=
fi

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

# setup parameters

if [ "$debug_flag" == "1" ]; then
echo debug_flag is set to 1
if [ "$FDS_DEBUG" == "1" ]; then
exit
fi
fi

fulldir=$BASEDIR/$dir
in=$infile.fds
outerr=$fulldir/$infile.err
outlog=$fulldir/$infile.log
stopfile=$infile.stop

scriptfile=$scratchdir/script.$$

# ensure that various files and directories exist

if ! [ -e $FDS ];  then
  if [ "$USEFDS" == "yes" ]; then
    echo "The file $FDS does not exist. Run aborted"
  else
    echo "Warning: The file $FDS does not exist."
  fi
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
#PBS -N $JOBPREFIX$infile -e $outerr -o $outlog
#PBS -l nodes=1:ppn=$nthreads
#PBS -l walltime=04:00:00
#PBS -l pvmem=1GB
#SBATCH -J $JOBPREFIX$infile
#SBATCH --mem-per-cpu=1000
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=$nthreads
#SBATCH -t 04:00:00
#SBATCH -e $outerr
#SBATCH -o $outlog
#SBATCH -p $queue

cd $fulldir

# Set number of OpenMP threads on target machine
export OMP_NUM_THREADS=$nthreads

echo Time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`

$RUNCMD $FDS $in 

# Run by $QSUB $QOPT $queue $scriptfile
EOF

echo Running `basename $FDS` $in 
if [ "$background" != "yes" ]; then
  chmod +x $scriptfile
  $QSUB $QOPT $queue $scriptfile
else
  cd $fulldir
  $QSUB $FDS $in
fi

rm $scriptfile
