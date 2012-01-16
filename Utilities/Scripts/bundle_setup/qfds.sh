#!/bin/bash
EXPECTED_ARGS=1

if [ $# -lt $EXPECTED_ARGS ]
then
  echo "Usage: qfds.sh [-f repository root] [-n processes per node] [-q queue]"
  echo "               [-r] [-t nthreads] [fds_command] casename.fds"
  echo ""
  echo "Runs an FDS case on a Linux cluster using the PBS/SGE qsub command"
  echo "This script can use the FDS executable specified on the command line"
  echo "or located in the users FDS repository (if the -r or -f option is "
  echo "specified.  Alternate queues may be specified using the -q option."
  echo ""
  echo " -n procs per node -  maximum processes assigned per node [default: 8]"
  echo " -t nthreads - number of threads used to run case [default: 1] "
  echo " -q queue - name of the queue. choices: batch (default), vis, fire60s,"
  echo "    fire70s"
  echo " -r - use repository executable"
  echo " -f repository root - repository where FDS executable"
  echo "    is located  [default: ~/FDS-SMV]"
  echo " fds_command - full path to fds command name (not used if -f or "
  echo "    -r options are specified)"
  echo "casename.fds - FDS input file"
  echo ""
  echo "Examples:"
  echo "qfds.sh -r casename.fds"
  echo "    run casename.fds using "
  echo "    ~/FDS-SMV/FDS_Compilation/intel_linux_64/fds_intel_linux_64"
  echo "qfds.sh ~/bin/fds_linux_64 casename.fds"
  echo "    run casename.fds using ~/bin/fds_linux_64"
  echo "qfds.sh -r -t 8 -n 2 casename.fds"
  echo "  run casename.fds using 8 processes, 2 processes per node (4 nodes)"
  echo "  and ~/FDS-SMV/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64"
  echo ""
  exit
fi
queue=batch
nthreads=1
use_repository=0
nprocs=8
FDSROOT=~/FDS-SMV
MPIRUN=
ndefined=0
while getopts 'f:n:q:rt:' OPTION
do
case $OPTION  in
  f)
   FDSROOT="$OPTARG"
   use_repository=1
   ;;
  n)
   nprocs="$OPTARG"
   ndefined=1
  ;;
  q)
   queue="$OPTARG"
  ;;
  r)
   use_repository=1
  ;;
  t)
   nthreads="$OPTARG"
  ;;
esac
done
shift $(($OPTIND-1))

if [ "$queue" == "fire60s" ]
then
  if [ $ndefined -eq 0 ]
  then
    nprocs=4
  fi
fi
if [ $use_repository -eq 0 ]
then
  fdsexe=$1
  in=$2
else
 if [ $nthreads -gt 1 ]
  then
  fdsexe=$FDSROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
 else
  fdsexe=$FDSROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
 fi
 in=$1
fi
if [ $nthreads -gt 1 ]
then
MPIRUN="mpirun -np $nthreads"
fi
echo MPIRUN=$MPIRUN

echo fdsexe=$fdsexe
echo input=$in
echo queue=$queue
echo nthreads=$nthreads
echo nprocs=$nprocs
echo use_repository=$use_repository


nprocs=8
nnodes=$(echo "($nthreads-1)/$nprocs+1" | bc)
if test $nnodes -le 0
then
nnodes=1
fi

infile=${in%.*}
fulldir=`pwd`

out=$fulldir/$infile.err
outlog=$fulldir/$infile.log

scriptfile=/tmp/script.$$
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

cat << EOF > $scriptfile
#!/bin/bash -f
#PBS -N $infile(MPI)
#PBS -e $out
#PBS -o $outlog
#PBS -l nodes=$nnodes:ppn=$nprocs
#\$ -N $infile(MPI)
#\$ -e $out
#\$ -o $outlog
#\$ -l nodes=$nnodes:ppn=$nprocs

cd $fulldir
echo Start time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`

$MPIRUN $fdsexe $in
EOF
chmod +x $scriptfile
echo Running $in using $queue queue
#qsub -q $queue $scriptfile
#rm $scriptfile
