#!/bin/bash
# $Date$ 
# $Revision$
# $Author$
#
PROG=$0

# setup default queue name

progname=qfds.sh
queue=batch
stopjob=0

nprocesses=1
nprocesses_per_node=1
nthreads=1

if [ $# -lt 1 ]
then
  echo "Usage: $progname [-d directory] [-f repository root] [-n processes per node] [-o nthreads]"
  echo "                 [-q queue] [-r] [-p nprocesses] [fds_command] casename.fds"
  echo ""
  echo "This script runs 64 bit serial or parallel versions of FDS using an executable"
  echo "specified on the command line or FDS from the respository if -r is specified."
  echo "The parallel FDS is invoked by using -p to specifying multiple processes."
  echo "Alternate queues (vis, fire70s) are set using the -q option."
  echo ""
  echo " -b use debug version"
  echo " -d directory [default: .]"
  echo " -n processes per node - maximum number of processes per node [default: 1]"
  echo "    (serial: 1, parallel: 8 for new cluster and fire70s, 4 for the vis queues)"
  echo " -o nthreads - run FDS (OpenMP) with a specified number of threads [default: $nthreads]"
  echo " -p nprocesses - number of processes used to run a case [default: 1] "
  echo " -q queue - name of the queue. choices: [default: $queue (other choices:"  
  echo "    vis and fire70s)"
  echo " -r - use FDS located in repository"
  echo " -s stop job"
  echo " -t - used for timing studies, run a job alone on a node"
  echo " -f repository root - name and location of repository where FDS is located"
  echo "    [default: ~/FDS-SMV]"
  echo " command - full path to command name (not used if either -f or -r"
  echo "    options are specified)"
  echo "input_file - input file"
  echo ""
  exit
fi

# default parameter settings

use_repository=0
use_debug=0
FDSROOT=~/FDS-SMV
MPIRUN=
dir=.
exe2=
ABORTRUN=n
IB=
DB=
SCRIPTFILE=
benchmark=no

if [ "$FDSNETWORK" == "infiniband" ] ; then
IB=ib
fi

# read in parameters from command line

while getopts 'bd:f:m:n:o:p:q:rst' OPTION
do
case $OPTION  in
  b)
   use_debug=1
   ;;
  d)
   dir="$OPTARG"
   ;;
  f)
   FDSROOT="$OPTARG"
   use_repository=1
   ;;
  m)
   SCRIPTFILE="$OPTARG"
   ;;
  n)
   nprocesses_per_node="$OPTARG"
   ;;
  o)
   nthreads="$OPTARG"
   ;;
  p)
   nprocesses="$OPTARG"
   ;;
  q)
   queue="$OPTARG"
   ;;
  r)
   use_repository=1
   ;;
  s)
   stopjob=1
   ;;
  t)
   benchmark="yes"
   ;;
esac
done
shift $(($OPTIND-1))

if [ "$use_debug" == "1" ] ; then
DB=_db
fi

if [ $use_repository -eq 0 ]
then
#set fds and the input file using the command line
  exe=$1
  in=$2
else
 if [ $nprocesses -gt 1 ]
# only set the input file using the command line, the fds exe is defined
# using the repository (serial if nprocesses==1 parallel otherwise)
  then
  exe=$FDSROOT/FDS_Compilation/mpi_intel_linux_64$IB$DB/fds_mpi_intel_linux_64$IB$DB
 else
  exe=$FDSROOT/FDS_Compilation/intel_linux_64$DB/fds_intel_linux_64$DB
 fi
 in=$1
fi

infile=${in%.*}

# if there is more than 1 process then use the mpirun command

TITLE="$infile"

if [ $nprocesses -gt 1 ] ; then
  MPIRUN="$MPIDIST/bin/mpirun -np $nprocesses"
  TITLE="$infile(MPI)"
  case $FDSNETWORK in
    "infiniband") TITLE="$infile(MPI_IB)"
  esac
fi

nnodes=$(echo "($nprocesses-1)/$nprocesses_per_node+1" | bc)
if test $nnodes -le 0
then
  nnodes=1
elif test $nnodes -gt 32
then
  nnodes=32
fi

# in benchmark mode run a case "alone" on one node
if [ "$benchmark" == "yes" ]; then
  nodes=1
  nprocesses_per_node=8
fi

cd $dir
fulldir=`pwd`

out=$fulldir/$infile.err
outlog=$fulldir/$infile.log
stopfile=$fulldir/$infile.stop


in_full_file=$fulldir/$in

# make sure files that are needed exist

if ! [ -e $in_full_file ]; then
  echo "The input file, $in_full_file, does not exist. Run aborted."
  ABORTRUN=y
fi
if ! [ -e $exe ]; then
  echo "The program, $exe, does not exist. Run aborted."
  ABORTRUN=y
fi
if [ -e $outlog ]; then
  echo "Removing log file: $outlog"
  rm $outlog
fi
if [ "$ABORTRUN" == "y" ] ; then
  exit
fi
if [ $STOPFDS ]; then
 echo "stopping case: $in"
 touch $stopfile
 exit
fi
if [ "$stopjob" == "1" ]; then
 echo "stopping case: $in"
 touch $stopfile
 exit
fi
if [ -e $stopfile ]; then
 rm $stopfile
fi

QSUB="qsub -q $queue"
if [ "$queue" == "terminal" ] ; then
  QSUB=
  MPIRUN=
fi

scriptfile=/tmp/script.$$
cat << EOF > $scriptfile
#!/bin/bash -f
#PBS -N $TITLE
#PBS -e $out
#PBS -o $outlog
#PBS -l nodes=$nnodes:ppn=$nprocesses_per_node
#\$ -N $TITLE
#\$ -e $out
#\$ -o $outlog
#\$ -l nodes=$nnodes:ppn=$nprocesses_per_node

export OMP_NUM_THREADS=$nthreads

cd $fulldir
echo Start time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`
$MPIRUN $exe $exe2 $in
EOF
echo "        Input file:$in"
echo "        Executable:$exe"
echo "             Queue:$queue"
echo "         Processes:$nprocesses"
if test $nprocesses -gt 1
then
echo "             Nodes:$nnodes"
echo "Processes per node:$nprocesses_per_node"
fi
chmod +x $scriptfile
$QSUB $scriptfile
rm $scriptfile
