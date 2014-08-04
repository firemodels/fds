#!/bin/bash
# $Date: 2014-08-01 12:07:22 -0400 (Fri, 01 Aug 2014) $ 
# $Revision: 20080 $
# $Author: gforney $
#


if [ $# -lt 1 ]
then
  echo "Usage: qfds.sh [-d directory] [-f repository root] [-n mpi processes per node] [-o nopenmp_threads]"
  echo "                 [-q queue] [-p nmpi_processes] [-e fds_command] casename.fds"
  echo ""
  echo "qfds.sh runs FDS using an executable specified with the -e option or"
  echo "from the respository if -e is not specified (the -r option is no longer" 
  echo "used.  A parallel version of FDS is invoked by using -p to specify the"
  echo "number of MPI processes and/or -o to specify the number of OpenMP threads."
  echo ""
  echo " -b     - use debug version of FDS"
  echo " -d dir - specify directory where the case is found [default: .]"
  echo " -e exe - full path of FDS used to run case"
  echo " -i   - output script file used to run case (case is not run)"
  echo " -m m - reserve m processes per node [default: 1]"
  echo " -n n - number of MPI processes per node [default: 1]"
  echo " -o o - number of OpenMP threads per process [default: 1]"
  echo " -p p - number of MPI processes [default: 1] "
  echo " -q q - name of queue. [default: batch]"  
  echo " -s   - stop job"
  echo " -t   - used for timing studies, run a job alone on a node"
  echo "input_file - input file"
  echo ""
  exit
fi

# default parameter settings

if [ "$FDSROOT" == "" ] ; then
  FDSROOT=~/FDS-SMV
fi
MPIRUN=
ABORTRUN=n
IB=
DB=
if [ "$FDSNETWORK" == "infiniband" ] ; then
  IB=ib
fi

# --------------------------- parse options --------------------

# default parameter settings

queue=batch
stopjob=0

nmpi_processes=1
nmpi_processes_per_node=1
maxmpi_processes_per_node=1
nopenmp_threads=1

use_repository=0
use_debug=0
dir=.
benchmark=no
showinput=0
use_repository=1

# read in parameters from command line

while getopts 'bd:e:im:n:o:p:q:st' OPTION
do
case $OPTION  in
  b)
   use_debug=1
   ;;
  d)
   dir="$OPTARG"
   ;;
  e)
   exe="$OPTARG"
   use_repository=0
   ;;
  i)
   showinput=1
   ;;
  m)
   maxmpi_processes_per_node="$OPTARG"
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
  s)
   stopjob=1
   ;;
  t)
   benchmark="yes"
   ;;
esac
done
shift $(($OPTIND-1))

# ^^^^^^^^^^^^^^^^^^^^^^^^parse options^^^^^^^^^^^^^^^^^^^^^^^^^

if [ "$use_debug" == "1" ] ; then
  DB=_db
fi

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
in=$1

infile=${in%.*}

# if there is more than 1 process then use the mpirun command

TITLE="$infile"

if [ $nmpi_processes -gt 1 ] ; then
  MPIRUN="$MPIDIST/bin/mpirun --map-by ppr:$nmpi_processes_per_node:node -np $nmpi_processes"
  TITLE="$infile(MPI)"
  case $FDSNETWORK in
    "infiniband") TITLE="$infile(MPI_IB)"
  esac
fi

let "nodes=($nmpi_processes-1)/$nmpi_processes_per_node+1"
let "ppn=($nopenmp_threads)*($nmpi_processes_per_node)"
if test $maxmpi_processes_per_node -gt $ppn ; then
  ppn=$maxmpi_processes_per_node
fi

if test $nodes -le 0 ; then
  nodes=1
fi

# in benchmark mode run a case "alone" on one node

if [ "$benchmark" == "yes" ]; then
  nodes=1
# use 8 on blaze cluster
  nmpi_processes_per_node=8
# use 12 on burn cluster
#  nmpi_processes_per_node=12
fi

cd $dir
fulldir=`pwd`

out=$fulldir/$infile.err
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

# create script file

scriptfile=/tmp/script.$$
cat << EOF > $scriptfile
#!/bin/bash
#PBS -N $TITLE
#PBS -e $out
#PBS -o $outlog
#PBS -l nodes=$nodes:ppn=$ppn
#\$ -N $TITLE
#\$ -e $out
#\$ -o $outlog
#\$ -l nodes=$nodes:ppn=$ppn

export OMP_NUM_THREADS=$nopenmp_threads

cd $fulldir
echo Start time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`
$MPIRUN $exe $in
EOF

# output info to terminal

echo "         Input file:$in"
echo "         Executable:$exe"
echo "              Queue:$queue"
echo "              Nodes:$nodes"
echo "          Processes:$nmpi_processes"
echo " Processes per node:$nmpi_processes_per_node"
if test $nmpi_processes -gt 1 ; then
  echo "Threads per process:$nopenmp_threads"
fi

# output script file to terminal

if [ "$showinput" == "1" ] ; then
  cat $scriptfile
  exit
fi

# run script

chmod +x $scriptfile
$QSUB $scriptfile
rm $scriptfile
