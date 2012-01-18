#!/bin/bash
# $Date$ 
# $Revision$
# $Author$

EXPECTED_ARGS=1

if [ $# -lt $EXPECTED_ARGS ]
then
  echo "Usage: qfds.sh [-f repository root] [-n processes per node] [-q queue]"
  echo "               [-r] [-p nprocesses] [fds_command] casename.fds"
  echo ""
  echo "This script runs FDS on a Linux cluster using an executable specified" 
  echo "on the command line or if the -r option is specified an executable"
  echo "located in the repository. The parallel version of FDS may be run by" 
  echo "using -p to specifying multiple processes. Alternate queues (vis, "
  echo "fire60s or fire70s) are set using the -q option."
  echo ""
  echo " -n processes per node - maximum number of processes assigned per"
  echo "    node [default: 8, (4 for the fire60s queue)]"
  echo " -p nprocesses - number of processes used to run a case [default: 1] "
  echo " -q queue - name of the queue. choices: [default: batch (other choices:"  
  echo "    vis, fire60s and fire70s)"
  echo " -r - use FDS located in repository"
  echo " -f repository root - name and location of repository where FDS "
  echo "    is located  [default: ~/FDS-SMV]"
  echo " fds_command - full path to fds command name (not needed if either"
  echo "    -f or -r options are specified)"
  echo "casename.fds - FDS input file"
  echo ""
  echo "Examples:"
  echo "qfds.sh -r casename.fds"
  echo "  run casename.fds using FDS located in the repository at"
  echo "  ~/FDS-SMV/FDS_Compilation/intel_linux_64/fds_intel_linux_64"
  echo "qfds.sh ~/bin/fds_linux_64 casename.fds"
  echo "  run casename.fds using FDS located at ~/bin/fds_linux_64"
  echo "qfds.sh -r -p 8 -n 2 casename.fds"
  echo "  run casename.fds using 8 processes, 2 processes per node (4 nodes)"
  echo "  and FDS located in the repository at"
  echo "  ~/FDS-SMV/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64"
  echo ""
  exit
fi

# default parameter settings

queue=batch
nprocesses=1
use_repository=0
nprocesses_per_node=8
FDSROOT=~/FDS-SMV
MPIRUN=
nprocesses_per_node_defined=0

# read in parameters from command line

while getopts 'f:n:p:q:r' OPTION
do
case $OPTION  in
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
esac
done
shift $(($OPTIND-1))

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
  fdsexe=$1
  in=$2
else
 if [ $nprocesses -gt 1 ]
# only set the input file using the command line, the fds exe is defined
# using the repository (serial if nprocesses==1 parallel otherwise)
  then
  fdsexe=$FDSROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
 else
  fdsexe=$FDSROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
 fi
 in=$1
fi

infile=${in%.*}

# if there is more than 1 process then use the mpirun command

if [ $nprocesses -gt 1 ]
then
MPIRUN="mpirun -np $nprocesses"
TITLE="$infile(MPI)"
else
TITLE="$infile"
fi

nnodes=$(echo "($nprocesses-1)/$nprocesses_per_node+1" | bc)
if test $nnodes -le 0
then
nnodes=1
fi

fulldir=`pwd`

out=$fulldir/$infile.err
outlog=$fulldir/$infile.log

# make sure files that are needed exist

if ! [ -e $fulldir/$in ]; then
  echo "The FDS input file, $fulldir/$in, does not exist. Run aborted."
fi
if ! [ -e $fdsexe ]; then
  echo "The FDS program name, $fdsexe, does not exist. Run aborted."
fi
if [ -e $outlog ]; then
  echo "Removing log file: $outlog"
  rm $outlog
fi
if ! [ -e $fulldir/$in ]; then
  exit
fi
if ! [ -e $fdsexe ]; then
  exit
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
$MPIRUN $fdsexe $in
EOF
echo "        Input file:$in"
echo "        Executable:$fdsexe"
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
