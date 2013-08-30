#!/bin/bash
# $Date$ 
# $Revision$
# $Author$
#
PROG=$0

# setup default queue name

progname=qfds.sh
queue=batch

nprocesses=1
nprocesses_per_node=1

if [ $# -lt 1 ]
then
  echo "Usage: $progname [-d directory] [-f repository root] [-n processes per node] [-q queue]"
  echo "               [-r] [-p nprocesses] [fds_command] casename.fds"
  echo ""
  echo "This script runs 64 bit serial or parallel versions of FDS using an executable"
  echo "specified on the command line or FDS from the respository if -r is specified."
  echo "The parallel FDS is invoked by using -p to specifying multiple processes."
  echo "Alternate queues (vis, fire60s or fire70s) are set using the -q option."
  echo ""
  echo " -b use debug version"
  echo " -d directory [default: .]"
  echo " -n processes per node - maximum number of processes per node [default: "
  echo "    (serial: 1, parallel: 8 for new cluster and fire70s, 4 for the fire60s" 
  echo "                          and vis queues)]"
  echo " -p nprocesses - number of processes used to run a case [default: 1] "
  echo " -q queue - name of the queue. choices: [default: $queue (other choices:"  
  echo "    vis, fire60s and fire70s)"
  echo " -r - use FDS (or Smokeview if -s is specified) located in repository"
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
nprocesses_per_node_defined=0
dir=.
# parameters used by Smokeview
USE_SMOKEVIEW=
VOLRENDER=
SKIPFRAME=1
STARTFRAME=0
exe2=
ABORTRUN=n
IB=
DB=
if [ "$FDSNETWORK" == "infiniband" ] ; then
IB=ib
fi

# read in parameters from command line

while getopts 'bd:f:n:p:q:rsxy:z:' OPTION
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
  n)
   nprocesses_per_node="$OPTARG"
   nprocesses_per_node_defined=1
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
   USE_SMOKEVIEW="y"
   ;;
  x)
   VOLRENDER=y
   ;;
  y)
   STARTFRAME="$OPTARG"
   ;;
  z)
   SKIPFRAME="$OPTARG"
   ;;

esac
done
shift $(($OPTIND-1))

#  if smokeview is invoked then override various options that 
#  may have been specified  (ie can't use batch queue, must use
#  smokeview bash script in repository)

if [ "$USE_SMOKEVIEW" == "y" ] ; then
  use_repository=1
  if [ "$queue" == "batch" ] ; then
    queue=fire70s
  fi
fi
if [ "$use_debug" == "1" ] ; then
DB=_db
fi

# set number of processes per node  to 4 if the fire60s queue is being used
# (the fire60s only have 4 cores)

if [ "$queue" == "fire60s" ]
then
if test $nprocesses_per_node_defined -eq 0 
then
  nprocesses_per_node=4
fi
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
  if [ "$USE_SMOKEVIEW" == "y" ] ; then
# for now only one instance of smokeview can occur per node
    nprocesses_per_node=1
    nprocesses=1
    exe="$FDSROOT/Verification/scripts/runsmv_single.sh"
    exe2="-x -y $STARTFRAME -z $SKIPFRAME"
  else
    exe=$FDSROOT/FDS_Compilation/intel_linux_64$DB/fds_intel_linux_64$DB
  fi
 fi
 in=$1
fi

infile=${in%.*}

# define options used by smokeview (for computing volume rendering smoke frames) 

if [[ "$USE_SMOKEVIEW" == "y" && "$VOLRENDER" == "y" ]] ; then
  VOLRENDER="-volrender"
  STARTFRAME="-startframe $STARTFRAME"
  SKIPFRAME="-skipframe $SKIPFRAME"
else
  VOLRENDER=
  STARTFRAME=
  SKIPFRAME=
fi

# if there is more than 1 process then use the mpirun command
#  (which will never happen if smokeview is running)

TITLE="$infile"

if [ $nprocesses -gt 1 ] ; then
  MPIRUN="$MPIDIST/bin/mpirun -np $nprocesses"
  TITLE="$infile(MPI)"
  case $FDSNETWORK in
    "infiniband") TITLE="$infile(MPI_IB)"
  esac
fi
if [ "$USE_SMOKEVIEW" == "y" ] ; then
  TITLE="$infile(SMV)"
fi

nnodes=$(echo "($nprocesses-1)/$nprocesses_per_node+1" | bc)
if test $nnodes -le 0
then
  nnodes=1
fi

cd $dir
fulldir=`pwd`

out=$fulldir/$infile.err
outlog=$fulldir/$infile.log
stopfile=$fulldir/$infile.stop


if [ "$USE_SMOKEVIEW" = "y" ] ; then
  in_full_file=$fulldir/$infile.smv
else
  in_full_file=$fulldir/$in
fi

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
if [ -e $stopfile ]; then
 rm $stopfile
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
qsub -q $queue $scriptfile
rm $scriptfile
