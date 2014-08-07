#!/bin/bash
# $Date$ 
# $Revision$
# $Author$
#


if [ $# -lt 1 ]
then
  echo "Usage: qfds.sh [-d directory] [-f repository root] [-n mpi processes per node] [-o nopenmp_threads]"
  echo "                 [-q queue] [-p nmpi_processes] [-e fds_command] casename.fds"
  echo ""
  echo "qfds.sh runs FDS using an executable specified with the -e option or"
  echo "from the respository if -e is not specified (the -r option is no longer" 
  echo "used).  A parallel version of FDS is invoked by using -p to specify the"
  echo "number of MPI processes and/or -o to specify the number of OpenMP threads."
  echo ""
  echo " -b     - use debug version of FDS"
  echo " -d dir - specify directory where the case is found [default: .]"
  echo " -e exe - full path of FDS used to run case"
  echo " -f repository root - name and location of repository where FDS is located"
  echo "    [default: ~/FDS-SMV]"
  echo " -m m - reserve m processes per node [default: 1]"
  echo " -n n - number of MPI processes per node [default: 1]"
  echo " -o o - number of OpenMP threads per process [default: 1]"
  echo " -p p - number of MPI processes [default: 1] "
  echo " -q q - name of queue. [default: batch]"
  echo "        If queue is terminal then job is run in the foreground on local computer"
  echo " -r   - report bindings"
  echo " -s   - stop job"
  echo " -t   - used for timing studies, run a job alone on a node"
  echo " -v   - list script used to run case to standard output"
  echo "input_file - input file"
  echo ""
  exit
fi

# default parameter settings

ncores=8
if [ "`uname`" != "Darwin" ]; then
  ncores=`grep processor /proc/cpuinfo | wc -l`
fi
FDSROOT=~/FDS-SMV
MPIRUN=
ABORTRUN=n
IB=
DB=
JOBPREFIX=
if [ "$FDSNETWORK" == "infiniband" ] ; then
  IB=ib
fi

# --------------------------- parse options --------------------

# default parameter settings

queue=batch
stopjob=0

nmpi_processes=1
nmpi_processes_per_node=1
max_processes_per_node=1
nopenmp_threads=1

use_repository=0
use_debug=0
dir=.
benchmark=no
showinput=0
use_repository=1
strip_extension=0
REPORT_BINDINGS=

# read in parameters from command line

while getopts 'bcd:e:f:j:m:n:o:p:q:rstv' OPTION
do
case $OPTION  in
  b)
   use_debug=1
   ;;
  c)
   strip_extension=1
   ;;
  d)
   dir="$OPTARG"
   ;;
  e)
   exe="$OPTARG"
   use_repository=0
   ;;
  f)
   FDSROOT="$OPTARG"
   ;;
  j)
   JOBPREFIX="$OPTARG"
   ;;
  m)
   max_processes_per_node="$OPTARG"
   ;;
  n)
   nmpi_processes_per_node="$OPTARG"
   ;;
  o)
   nopenmp_threads="$OPTARG"
   ;;
  p)
   nmpi_processes="$OPTARG"
   ;;
  q)
   queue="$OPTARG"
   ;;
  r)
   REPORT_BINDINGS="--report-bindings"
   ;;
  s)
   stopjob=1
   ;;
  t)
   benchmark="yes"
   ;;
  v)
   showinput=1
   ;;
esac
done
shift $(($OPTIND-1))

# ^^^^^^^^^^^^^^^^^^^^^^^^parse options^^^^^^^^^^^^^^^^^^^^^^^^^

if [ "$use_debug" == "1" ] ; then
  DB=_db
fi

# define executables if the repository is used

if [ $use_repository -eq 1 ] ; then
# use fds from repository (-e was not specified)
 if [ $nmpi_processes -gt 1 ] ; then
# use mpi version of fds 
  exe=$FDSROOT/FDS_Compilation/mpi_intel_linux_64$IB$DB/fds_mpi_intel_linux_64$IB$DB
 else
# use non-mpi version of fds 
  exe=$FDSROOT/FDS_Compilation/intel_linux_64$DB/fds_intel_linux_64$DB
 fi
fi

#define input file

in=$1
infile=${in%.*}

# if there is more than 1 process then use the mpirun command

TITLE="$infile"

# define number of nodes

let "nodes=($nmpi_processes-1)/$nmpi_processes_per_node+1"
if test $nodes -lt 1 ; then
  nodes=1
fi

# define processes per node

let "ppn=($nopenmp_threads)*($nmpi_processes_per_node)"
if test $ppn -le $max_processes_per_node ; then
  ppn=$max_processes_per_node
fi

# bind to sockets if OpenMP is being used (number of threads > 1)

SOCKET_OPTION=" "
if test $nopenmp_threads -gt 1 ; then
  SOCKET_OPTION="--bind-to socket --map-by ppr:1:socket"
fi

# in benchmark mode run a case "alone" on one node

if [ "$benchmark" == "yes" ]; then
  nodes=1
  nmpi_processes_per_node=$ncores
fi

# use mpirun if there is more than 1 process

if [ $nmpi_processes -gt 1 ] ; then
  MPIRUN="$MPIDIST/bin/mpirun $REPORT_BINDINGS $SOCKET_OPTION -np $nmpi_processes"
  TITLE="$infile(MPI)"
  case $FDSNETWORK in
    "infiniband") TITLE="$infile(MPI_IB)"
  esac
fi

cd $dir
fulldir=`pwd`

# define files

outerr=$fulldir/$infile.err
outlog=$fulldir/$infile.log
stopfile=$fulldir/$infile.stop
in_full_file=$fulldir/$in

# make sure various files exist before running case

if ! [ -e $in_full_file ]; then
  if [ "$showinput" == "0" ] ; then
    echo "The input file, $in_full_file, does not exist. Run aborted."
    ABORTRUN=y
  fi
fi
if [ "$strip_extension" == "1" ] ; then
  in=$infile
fi
if ! [ -e $exe ]; then
  if [ "$showinput" == "0" ] ; then
    echo "The program, $exe, does not exist. Run aborted."
    ABORTRUN=y
  fi
fi
if [ -e $outlog ]; then
  echo "Removing log file: $outlog"
  rm $outlog
fi
if [ "$ABORTRUN" == "y" ] ; then
  if [ "$showinput" == "0" ] ; then
    exit
  fi
fi
if [ $STOPFDS ]; then
 echo "stopping case: $in"
 touch $stopfile
 exit
fi
if [ "$STOPFDSMAXITER" != "" ]; then
  echo "creating delayed stop file: $infile"
  echo $STOPFDSMAXITER > $stopfile
fi
if [ "$stopjob" == "1" ]; then
 echo "stopping case: $in"
 touch $stopfile
 exit
fi
if [ "$STOPFDSMAXITER" == "" ]; then
  if [ -e $stopfile ]; then
    rm $stopfile
  fi
fi

QSUB="qsub -q $queue"

if [ "$queue" == "terminal" ] ; then
  QSUB=
  MPIRUN=
fi

if [ "$queue" == "none" ]; then
  if [ "$BACKGROUND" == "" ]; then
    BACKGROUND=background
  fi
  QSUB="$BACKGROUND -u 75 -d 10 "
  background=yes;
  notfound=`$BACKGROUND -help 2>&1 | tail -1 | grep "not found" | wc -l`
  if [ "$showinput" == "0" ]; then
    if [ "$notfound" == "1" ];  then
      echo "The program $BACKGROUND (background) is not available. Run aborted"
      exit
    fi
  fi
fi

# create script file for queuing systems

scriptfile=/tmp/script.$$


cat << EOF > $scriptfile
#!/bin/bash
EOF

if [ "$queue" != "none" ] ; then
cat << EOF >> $scriptfile
#PBS -N $JOBPREFIX$TITLE
#PBS -e $outerr
#PBS -o $outlog
#PBS -l nodes=$nodes:ppn=$ppn
#SBATCH -J $JOBPREFIX$infile
#SBATCH -e $outerr
#SBATCH -o $outlog
#SBATCH -p $queue
EOF
fi

cat << EOF >> $scriptfile
export OMP_NUM_THREADS=$nopenmp_threads

cd $fulldir
echo Start time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`
$MPIRUN $exe $in
EOF

# if requested, output script file to screen

if [ "$showinput" == "1" ] ; then
  cat $scriptfile
  exit
fi

# output info to screen

if [ "$queue" != "none" ] ; then
  echo "         Input file:$in"
  echo "         Executable:$exe"
  echo "              Queue:$queue"
  echo "              Nodes:$nodes"
  echo "          Processes:$nmpi_processes"
  echo " Processes per node:$nmpi_processes_per_node"
  if test $nmpi_processes -gt 1 ; then
    echo "Threads per process:$nopenmp_threads"
  fi
fi

# run script

chmod +x $scriptfile
$QSUB $scriptfile
if [ "$queue" != "none" ] ; then
  rm $scriptfile
fi

