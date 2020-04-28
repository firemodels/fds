#!/bin/bash

#*** environment varables

#*** environment variables used by the bots
# BACKGROUND_PROG  - defines location of background program
#                    ( if the 'none' queue is also specified)

# ---------------------------- stop_fds_if_requested ----------------------------------

function stop_fds_if_requested {
if [ "$OPENMPCASES" == "" ]; then
  if [ "$STOPFDS" != "" ]; then
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
else
  for i in `seq 1 $OPENMPCASES`; do
    stopfile=${filebase[$i]}.stop
    if [ "$STOPFDS" != "" ]; then
      echo "stopping case: ${files[$i]}"
      touch $stopfile
    fi

    if [ "$STOPFDSMAXITER" != "" ]; then
      echo "creating delayed stop file: $stopfile"
      echo $STOPFDSMAXITER > $stopfile
    fi

    if [ "$stopjob" == "1" ]; then
      echo "stopping case: ${files[$i]}"
      touch $stopfile
    fi

    if [ "$STOPFDSMAXITER" == "" ]; then
      if [ "$STOPFDS" == "" ]; then
        if [ -e $stopfile ]; then
          rm $stopfile
        fi
      fi
    fi
  done
  if [ "$STOPFDS" != "" ]; then
    exit
  fi
  if [ "$stopjob" == "1" ]; then
    exit
  fi
fi
}

# ---------------------------- usage ----------------------------------

function usage {
  if [ "$use_intel_mpi" == "1" ]; then
    MPI=impi
  else
    MPI=mpi
  fi
  echo "Usage: qfds.sh [-p n_mpi_processes] [-o nthreads] [-e fds_command] [-q queue]  casename.fds"
  echo ""
  echo "qfds.sh runs FDS using an executable from the repository or one specified with the -e option."
  echo "A parallel version of FDS is invoked by using -p to specify the number of MPI processes and/or"
  echo "-o to specify the number of OpenMP threads."
  echo ""
  echo "qfds.sh loads the modules that were loaded when fds was built unless the -C option is specified"
  echo "then the currently loaded modules are used."
  echo ""
  echo " -e exe - full path of FDS used to run case "
  echo "    [default: $FDSROOT/fds/Build/${MPI}_intel_${platform}_64$DB/fds_${MPI}_intel_${platform}_64$DB]"
  echo " -h   - show commonly used options"
  echo " -H   - show all options"
  echo " -o o - number of OpenMP threads per process [default: 1]"
  echo " -p p - number of MPI processes [default: 1] "
  echo " -q q - name of queue. [default: batch]"
  echo " -v   - output generated script to standard output"
  echo "input_file - input file"
  if [ "$HELP" == "" ]; then
    exit
  fi
  echo "Other options:"
  echo " -c file - loads Intel Trace Collector configuration file "
  echo " -C   - use modules currently loaded rather than modules loaded when fds was built."
  echo " -d dir - specify directory where the case is found [default: .]"
  echo " -D n - delay the case submission by n seconds"
  echo " -E - use tcp transport (only available with the Intel compiled versions of fds)"
  echo "      This options adds export I_MPI_FABRICS=shm:tcp to the run script"
  echo " -f repository root - name and location of repository where FDS is located"
  echo "    [default: $FDSROOT]"
  echo " -i use installed fds"
  echo " -I use Intel MPI version of fds"
  echo " -j prefix - specify a job prefix"
  echo " -L use Open MPI version of fds"
  echo " -m m - reserve m processes per node [default: 1]"
  echo " -M   -  add --mca plm_rsh_agent /usr/bin/ssh to mpirun command "
  echo " -n n - number of MPI processes per node [default: 1]"
  echo " -N   - run as many MPI processes on a node as possible"
  echo "        MIN ( number of cores, number of mpi processes)"
  echo " -O n - run cases casea.fds, caseb.fds, ... using 1, ..., N OpenMP threads"
  echo "        where case is specified on the command line. N can be at most 9."
  echo " -s   - stop job"
  echo " -S   - use startup files to set the environment, do not load modules"
  echo " -r   - append trace flag to the mpiexec call generated"
  echo " -t   - used for timing studies, run a job alone on a node (reserving $NCORES_COMPUTENODE cores)"
  echo " -T type - run dv (development), db (debug), inspect, advise, or vtune version of fds"
  echo "           if -T is not specified then the release version of fds is used"
  echo " -V   - show command line used to invoke qfds.sh"
  echo " -w time - walltime, where time is hh:mm for PBS and dd-hh:mm:ss for SLURM. [default: $walltime]"
  echo ""
  echo " Resource manager: $RESOURCE_MANAGER"
  exit
}

#*** get directory containing qfds.sh

QFDS_PATH=$(dirname `which $0`)
CURDIR=`pwd`
cd $QFDS_PATH
QFDS_DIR=`pwd`
cd $CURDIR
SLEEP=

#*** define toplevel of the repos

FDSROOT=~/FDS-SMV
if [ "$FIREMODELS" != "" ]; then
  FDSROOT=$FIREMODELS
fi

#*** determine platform

platform="linux"
if [ "$WINDIR" != "" ]; then
  platform="win"
  if [ "$I_MPI_ROOT" != "" ]; then
    echo $I_MPI_ROOT | sed 's/\\/\//g' -  | sed 's/C:/\/C/g' -  > var.out.$$
    I_MPI_ROOT=`head -1 var.out.$$`
    rm var.out.$$
  fi
fi
if [ "`uname`" == "Darwin" ] ; then
  platform="osx"
fi

#*** determine number of cores and default queue

if [ "$platform" == "osx" ]; then
  queue=none
  ncores=`system_profiler SPHardwareDataType|grep Cores|awk -F' ' '{print $5}'`
else
  if [ "$platform" == "win" ]; then
    queue=terminal
  else
    queue=batch
  fi
  ncores=`grep processor /proc/cpuinfo | wc -l`
fi
if [ "$NCORES_COMPUTENODE" == "" ]; then
  NCORES_COMPUTENODE=$ncores
else
  ncores=$NCORES_COMPUTENODE
fi

#*** set default parameter values

showcommandline=
HELP=
FDS_MODULE_OPTION=1
MPIRUN=
ABORTRUN=n
DB=
OUT2ERROR=
stopjob=0
MCA=
OPENMPCASES=
OPENMPTEST=
if [ "$MPIRUN_MCA" != "" ]; then
  MCA=$MPIRUN_MCA
fi

n_mpi_processes=1
n_mpi_processes_per_node=1
if [ "$platform" == "linux" ]; then
max_processes_per_node=`cat /proc/cpuinfo | grep cores | wc -l`
else
max_processes_per_node=1
fi
n_openmp_threads=1
trace=
use_installed=
use_debug=
use_devel=
use_inspect=
use_advise=
use_vtune=
use_intel_mpi=1
iinspectresdir=
iinspectargs=
vtuneresdir=
vtuneargs=
use_config=""

# determine which resource manager is running (or none)

STATUS_FILE=status_file.$$
srun -V >& $STATUS_FILE
missing_slurm=`cat $STATUS_FILE | tail -1 | grep "not found" | wc -l`
rm -f $STATUS_FILE

RESOURCE_MANAGER="NONE"
if [ $missing_slurm -eq 0 ]; then
  RESOURCE_MANAGER="SLURM"
else
  echo | qmgr -n >& $STATUS_FILE
  missing_torque=`cat $STATUS_FILE | tail -1 | grep "not found" | wc -l`
  rm -f $STATUS_FILE
  if [ $missing_torque -eq 0 ]; then
    RESOURCE_MANAGER="TORQUE"
  fi
fi
if [ "$RESOURCE_MANAGER" == "SLURM" ]; then
  if [ "$SLURM_MEM" != "" ]; then
   SLURM_MEM="#SBATCH --mem=$SLURM_MEM"
  fi
  if [ "$SLURM_MEMPERCPU" != "" ]; then
   SLURM_MEM="#SBATCH --mem-per-cpu=$SLURM_MEMPERCPU"
  fi
fi

# the mac doesn't have Intel MPI
if [ "`uname`" == "Darwin" ]; then
  use_intel_mpi=
fi
dir=.
benchmark=no
showinput=0
exe=
STARTUP=
SET_MPI_PROCESSES_PER_NODE=
if [ "$QFDS_STARTUP" != "" ]; then
  STARTUP=$QFDS_STARTUP
fi

if [ $# -lt 1 ]; then
  usage
fi

if [ "$BACKGROUND_PROG" == "" ]; then
   BACKGROUND_PROG=background
fi
if [ "$BACKGROUND_DELAY" == "" ]; then
   BACKGROUND_DELAY=10
fi
if [ "$BACKGROUND_LOAD" == "" ]; then
   BACKGROUND_LOAD=75
fi

commandline=`echo $* | sed 's/-V//' | sed 's/-v//'`

#*** read in parameters from command line

while getopts 'Aa:c:Cd:D:e:Ef:hHiIj:Lm:Mn:No:O:p:Pq:rsStT:vVw:x:' OPTION
do
case $OPTION  in
  A) # used by timing scripts to identify benchmark cases
   DUMMY=1
   ;;
  a)
   vtuneresdir="$OPTARG"
   use_vtune=1
   ;;
  c)
   use_config="$OPTARG"
   ;;
  C)
   FDS_MODULE_OPTION=
   ;;
  d)
   dir="$OPTARG"
   ;;
  D)
   SLEEP="sleep $OPTARG"
   ;;
  e)
   exe="$OPTARG"
   ;;
  E)
   TCP=1
   ;;
  f)
   FDSROOT="$OPTARG"
   ;;
  h)
   usage
   exit
   ;;
  H)
   HELP=ALL
   usage
   exit
   ;;
  i)
   use_installed=1
   ;;
  I)
   use_intel_mpi=1
   ;;
  j)
   JOBPREFIX="$OPTARG"
   ;;
  L)
   use_intel_mpi=
   ;;
  M)
   MCA="--mca plm_rsh_agent /usr/bin/ssh "
   ;;
  m)
   max_processes_per_node="$OPTARG"
   ;;
  n)
   n_mpi_processes_per_node="$OPTARG"
   ;;
  N)
   SET_MPI_PROCESSES_PER_NODE=1
   ;;
  o)
   n_openmp_threads="$OPTARG"
   ;;
  O)
   OPENMPCASES="$OPTARG"
   if [ $OPENMPCASES -gt 9 ]; then
     OPENMPCASES=9
   fi
   n_mpi_process=1
   benchmark="yes"
   if [ "$RESOURCE_MANAGER" == "TORQUE" ]; then
     if [ "$NCORES_COMPUTENODE" != "" ]; then
       n_mpi_processes_per_node="$NCORES_COMPUTENODE"
     fi
   fi
   ;;
  p)
   n_mpi_processes="$OPTARG"
   ;;
  P)
   OPENMPCASES="2"
   OPENMPTEST="1"
   benchmark="yes"
   n_mpi_process=1
   if [ "$RESOURCE_MANAGER" == "TORQUE" ]; then
     if [ "$NCORES_COMPUTENODE" != "" ]; then
       n_mpi_processes_per_node="$NCORES_COMPUTENODE"
     fi
   fi  
   ;;
  q)
   queue="$OPTARG"
   ;;
  r)
   trace="-trace"
   ;;
  s)
   stopjob=1
   ;;
  S)
   STARTUP=1
   ;;
  t)
   benchmark="yes"
   if [ "$RESOURCE_MANAGER" == "TORQUE" ]; then
     if [ "$NCORES_COMPUTENODE" != "" ]; then
       n_mpi_processes_per_node="$NCORES_COMPUTENODE"
     fi
   fi
   ;;
  T)
   TYPE="$OPTARG"
   use_devel=
   use_debug=
   if [ "$TYPE" == "dv" ]; then
     use_devel=1
   fi
   if [ "$TYPE" == "db" ]; then
     use_debug=1
   fi
   if [ "$TYPE" == "inspect" ]; then
     use_inspect=1
   fi
   if [ "$TYPE" == "advise" ]; then
     use_advise=1
   fi
   if [ "$TYPE" == "vtune" ]; then
     use_vtune=1
   fi
   ;;
  v)
   showinput=1
   ;;
  V)
   showcommandline=1
   ;;
  w)
   walltime="$OPTARG"
   ;;
  x)
  iinspectresdir="$OPTARG"
   use_inspect=1
   ;;
   
esac
done
shift $(($OPTIND-1))

if [ "$showcommandline" == "1" ]; then
  echo $0 $commandline
  exit
fi

if [ "$SET_MPI_PROCESSES_PER_NODE" == "1" ]; then
   n_mpi_processes_per_node=$n_mpi_processes
   if test $n_mpi_processes_per_node -gt $ncores ; then
     n_mpi_processes_per_node=$ncores
   fi
fi

#*** define input file

in=$1
infile=${in%.*}

if [[ "$TCP" != "" ]] && [[ "$use_intel_mpi" == "" ]]; then
  echo "***error: -The E option for specifying tcp transport is only available"
  echo "          with Intel the compiled versions of fds"
  exit
fi

if [ "$OPENMPCASES" == "" ]; then
  files[1]=$in
  filebase[1]=$infile
  nthreads[1]=$n_openmp_threads
else
  for i in `seq 1 $OPENMPCASES`; do
    nthreads[$i]=$i
    if [[ "$OPENMPTEST" == "1" ]] && [[ "$i" == "2" ]]; then
      nthreads[$i]=4
    fi
    arg=`echo ${nthreads[$i]} | tr 123456789 abcdefghi`
    filebase[$i]=$in$arg
    files[$i]=$in$arg.fds
  done
fi

#*** parse options

if [ "$walltime" == "" ]; then
    if [ "$RESOURCE_MANAGER" == "SLURM" ]; then
	walltime=99-99:99:99
    else
	walltime=999:0:0
    fi
fi

#*** define executable

if [ "$use_installed" == "1" ]; then
  notfound=`echo | fds 2>&1 >/dev/null | tail -1 | grep "not found" | wc -l`
  if [ $notfound -eq 1 ]; then
    echo "fds is not installed. Run aborted."
    ABORTRUN=y
    exe=
  else
    fdspath=`which fds`
    fdsdir=$(dirname "${fdspath}")
    curdir=`pwd`
    cd $fdsdir
    exe=`pwd`/fds
    cd $curdir
  fi
else
  if [ "$use_debug" == "1" ]; then
    DB=_db
  fi
  if [ "$use_devel" == "1" ]; then
    DB=_dv
  fi
  if [ "$use_inspect" == "1" ]; then
    DB=_inspect
    if [ "$iinspectresdir" != "" ]; then
    iinspectargs="inspxe-cl -collect ti2 -knob stack-depth=32 -s-f $FDSROOT/fds/Build/impi_intel_linux_64_inspect/suppressions/default.sup  -result-dir $iinspectresdir --"
    fi
  fi
  if [ "$use_advise" == "1" ]; then
    DB=_advise
  fi
  if [ "$use_vtune" == "1" ]; then
    DB=_vtune
    if [ "$vtuneresdir" != "" ]; then
    vtuneargs="amplxe-cl -collect hpc-performance -result-dir $vtuneresdir --"
    fi
  fi
  if [ "$use_intel_mpi" == "1" ]; then
    if [ "$exe" == "" ]; then
      exe=$FDSROOT/fds/Build/impi_intel_${platform}_64$DB/fds_impi_intel_${platform}_64$DB
    fi
  fi
  if [ "$exe" == "" ]; then
    exe=$FDSROOT/fds/Build/mpi_intel_${platform}_64$DB/fds_mpi_intel_${platform}_64$DB
  fi
fi

#*** check to see if fds was built using Intel MPI

if [ -e $exe ]; then
  if [ "$use_intel_mpi" == "" ]; then
    is_intel_mpi=`echo "" | $exe 2>&1 >/dev/null | grep MPI | grep library | grep Intel | wc -l`
    if [ "$is_intel_mpi" == "1" ]; then
         use_intel_mpi=1
    fi
  fi
fi

#*** modules loaded currently

if [ "$STARTUP" == "" ]; then

  CURRENT_LOADED_MODULES=`echo $LOADEDMODULES | tr ':' ' '`

# modules loaded when fds was built

  if [ "$exe" != "" ]; then  # first look for file that contains the list
    FDSDIR=$(dirname "$exe")
    if [ -e $FDSDIR/.fdsinfo ]; then
      FDS_LOADED_MODULES=`tail -1 $FDSDIR/.fdsinfo`
      OPENMPI_PATH=`head -1 $FDSDIR/.fdsinfo`
    fi
  fi

  if [[ "$FDS_MODULE_OPTION" == "1" ]] && [[ "$FDS_LOADED_MODULES" != "" ]]; then
    MODULES=$FDS_LOADED_MODULES               # modules loaded when fds was built
  else
    MODULES=$CURRENT_LOADED_MODULES
  fi
fi

#*** define number of nodes

if test $n_mpi_processes_per_node -gt $ncores ; then
  n_mpi_processes_per_node=$ncores
fi

if test $n_mpi_processes_per_node = -1 ; then
  if test $n_mpi_processes -gt 1 ; then
    n_mpi_processes_per_node=2
  else
    n_mpi_processes_per_node=1
  fi
fi

let "nodes=($n_mpi_processes-1)/$n_mpi_processes_per_node+1"
if test $nodes -lt 1 ; then
  nodes=1
fi

#*** define processes per node

let ppn="$n_mpi_processes_per_node*n_openmp_threads"
if [ $ppn -gt $max_processes_per_node ]; then
  ppn=$max_processes_per_node
fi

if [[ $n_openmp_threads -gt 1 ]] && [[ "$use_intel_mpi" == "1" ]]; then
  ppn=2
fi

# Socket or node bindings

if test $n_mpi_processes -gt 1 ; then
 if test $n_openmp_threads -gt 1 ; then
  SOCKET_OPTION="--map-by socket:PE=$n_openmp_threads"
 else
  SOCKET_OPTION=" "
 fi
else
 SOCKET_OPTION="--map-by node:PE=$n_openmp_threads"
fi

if [ "$use_intel_mpi" == "1" ]; then
  SOCKET_OPTION=" "
fi

# Report bindings in the .err or .log file

if [ "$use_intel_mpi" == "1" ]; then
 REPORT_BINDINGS=" "
else
 REPORT_BINDINGS="--report-bindings"
fi

# The "none" queue does not use the queing system

if [ "$queue" == "none" ]; then
 SOCKET_OPTION=
 REPORT_BINDINGS=
fi

#*** define MPIRUNEXE and do some error checking

if [ "$use_intel_mpi" == "1" ]; then
  if [ "$use_installed" == "1" ]; then
    MPIRUNEXE=$fdsdir/INTEL/mpi/intel64/bin/mpiexec
    if [ ! -e $MPIRUNEXE ]; then
      MPIRUNEXE=$fdsdir/INTEL/bin/mpiexec
      if [ ! -e $MPIRUNEXE ]; then
        echo "Intel mpiexec not found at:"
        echo "$fdsdir/INTEL/mpi/intel64/bin/mpiexec or"
        echo "$fdsdir/bin/mpiexec"
        echo "Run aborted"
        ABORT=y
      fi
    fi
  else
    if [ "$I_MPI_ROOT" == "" ]; then
      echo "Intel MPI environment not setup. Run aborted."
      ABORTRUN=y
    else
      MPIRUNEXE=$I_MPI_ROOT/intel64/bin/mpiexec
      if [ ! -e "$MPIRUNEXE" ]; then
        MPIRUNEXE="$I_MPI_ROOT/bin/mpiexec"
        if [ ! -e "$MPIRUNEXE" ]; then
          echo "Intel mpiexec not found at:"
          echo "$I_MPI_ROOT/intel64/bin/mpiexec or"
          echo "$I_MPI_ROOT/bin/mpiexec"
          ABORTRUN=y
          echo "Run aborted"
        fi
      fi
    fi
  fi
else                                 # using OpenMPI
  if [ "$OPENMPI_PATH" != "" ]; then
    if [ -e $OPENMPI_PATH/bin ]; then
      mpibindir=$OPENMPI_PATH/bin/
    fi
  fi
  MPIRUNEXE=${mpibindir}mpirun
  if [ "$mpibindir" == "" ]; then  # OPENMPI_PATH blank so mpirun needs to be in path
    notfound=`$MPIRUNEXE -h 2>&1 >/dev/null | head -1 | grep "not found" | wc -l`
    if [ $notfound -eq 1 ]; then
      echo "*** error: $MPIRUNEXE not in PATH"
      ABORTRUN=y
    fi
  else                             # use full path to mpirun
    if [ ! -e $MPIRUNEXE ]; then
      echo "*** error: $MPIRUNEXE does not exist"
      ABORTRUN=y
    fi
  fi
fi

TITLE="$infile"
MPIRUN="$MPIRUNEXE $REPORT_BINDINGS $SOCKET_OPTION $MCA -np $n_mpi_processes $trace"

cd $dir
fulldir=`pwd`

#*** define files

outerr=$fulldir/$infile.err
outlog=$fulldir/$infile.log
qlog=$fulldir/$infile.qlog
stopfile=$fulldir/$infile.stop
scriptlog=$fulldir/$infile.slog
in_full_file=$fulldir/$in

#*** make sure various files exist before running the case

if [ "$OPENMPCASES" == "" ]; then
  if ! [ -e $in_full_file ]; then
    if [ "$showinput" == "0" ]; then
      echo "The input file, $in_full_file, does not exist. Run aborted."
      ABORTRUN=y
    fi
  fi
else
for i in `seq 1 $OPENMPCASES`; do
  in_full_file=$fulldir/${files[$i]}
  if ! [ -e $in_full_file ]; then
    if [ "$showinput" == "0" ]; then
      echo "The input file, $in_full_file, does not exist."
      ABORTRUN=y
    fi
  fi
done
if [ "$ABORTRUN" == "y" ]; then
  echo "Run aborted."
fi
fi

if [ "$STOPFDS" == "" ]; then
  if [ "$exe" != "" ]; then
    if ! [ -e "$exe" ]; then
      if [ "$showinput" == "0" ]; then
        echo "The program, $exe, does not exist. Run aborted."
        ABORTRUN=y
      fi
    fi
  fi

  if [ -e $outlog ]; then
    echo "Removing log file: $outlog"
    rm $outlog
  fi

  if [ "$ABORTRUN" == "y" ]; then
    if [ "$showinput" == "0" ]; then
      exit
    fi
  fi
fi

stop_fds_if_requested

QSUB="qsub -q $queue"

#*** use the queue none and the program background on systems
#    without a queing system

if [ "$queue" == "none" ]; then
  OUT2ERROR=" 2> $outerr"
  notfound=`$BACKGROUND_PROG -help 2>&1 | tail -1 | grep "not found" | wc -l`
  if [ "$showinput" == "0" ]; then
    if [ "$notfound" == "1" ];  then
      echo "The program $BACKGROUND_PROG was not found."
      echo "Install FDS which has the background utility."
      echo "Run aborted"
      exit
    fi
  fi
  MPIRUN=
  QSUB="$BACKGROUND_PROG -u $BACKGROUND_LOAD -d $BACKGROUND_DELAY "
  USE_BACKGROUND=1
else

  if [ "$queue" == "terminal" ]; then
    QSUB=
    MPIRUN=
  fi

#*** setup for SLURM (alternative to torque)

  if [ "$RESOURCE_MANAGER" == "SLURM" ]; then
    QSUB="sbatch -p $queue --ignore-pbs"
    MPIRUN="srun -N $nodes -n $n_mpi_processes --ntasks-per-node $n_mpi_processes_per_node"
  fi
fi

#*** Set walltime parameter only if walltime is specified as input argument

walltimestring_pbs=
walltimestring_slurm=
if [ "$walltime" != "" ]; then
  walltimestring_pbs="-l walltime=$walltime"
  walltimestring_slurm="-t $walltime"
fi

#*** create a random script file for submitting jobs

scriptfile=`mktemp /tmp/script.$$.XXXXXX`

cat << EOF > $scriptfile
#!/bin/bash
# $0 $commandline
EOF

if [ "$queue" != "none" ]; then
  if [ "$RESOURCE_MANAGER" == "SLURM" ]; then
    cat << EOF >> $scriptfile
#SBATCH -J $JOBPREFIX$infile
#SBATCH -e $outerr
#SBATCH -o $outlog
#SBATCH -p $queue
#SBATCH -n $n_mpi_processes
#SBATCH --nodes=$nodes
#SBATCH --cpus-per-task=$n_openmp_threads
#SBATCH --ntasks-per-node=$n_mpi_processes_per_node
EOF

if [ "$benchmark" == "yes" ]; then
cat << EOF >> $scriptfile
#SBATCH --exclusive
EOF
fi

cat << EOF >> $scriptfile
$SLURM_MEM
EOF
    if [ "$walltimestring_slurm" != "" ]; then
      cat << EOF >> $scriptfile
#SBATCH $walltimestring_slurm
EOF
    fi

  else
    cat << EOF >> $scriptfile
#PBS -N $JOBPREFIX$TITLE
#PBS -W umask=0022
#PBS -e $outerr
#PBS -o $outlog
#PBS -l nodes=$nodes:ppn=$ppn
EOF
    if [ "$walltimestring_pbs" != "" ]; then
      cat << EOF >> $scriptfile
#PBS $walltimestring_pbs
EOF
    fi
    if [[ $n_openmp_threads -gt 1 ]] && [[ "$use_intel_mpi" == "1" ]]; then
      cat << EOF >> $scriptfile
#PBS -l naccesspolicy=SINGLEJOB -n
EOF
    fi
  fi
fi

if [[ "$MODULES" != "" ]]; then
  cat << EOF >> $scriptfile
export MODULEPATH=$MODULEPATH
module purge
module load $MODULES
EOF
fi

if [ "$OPENMPCASES" == "" ]; then
cat << EOF >> $scriptfile
export OMP_NUM_THREADS=$n_openmp_threads
EOF
fi

if [ "$use_vtune" == "1" ]; then
cat << EOF >> $scriptfile
source /opt/intel19/vtune_amplifier_2019/amplxe-vars.sh quiet
EOF
fi

if [ "$use_inspect" == "1" ]; then
cat << EOF >> $scriptfile
source /opt/intel19/inspector_2019/inspxe-vars.sh quiet
EOF
fi


if [ "$use_intel_mpi" == "1" ]; then
cat << EOF >> $scriptfile
export I_MPI_DEBUG=5
EOF
fi

if [[ $n_openmp_threads -gt 1 ]] && [[ "$use_intel_mpi" == "1" ]]; then
cat << EOF >> $scriptfile
export I_MPI_PIN_DOMAIN=omp
EOF
fi

if [ "$TCP" != "" ]; then
  cat << EOF >> $scriptfile
export I_MPI_FABRICS=shm:tcp
EOF
fi

if [ "$use_config" != "" ]; then
  cat << EOF >> $scriptfile
export VT_CONFIG=$use_config
EOF
fi

cat << EOF >> $scriptfile
cd $fulldir
echo
echo \`date\`
EOF

if [ "$OPENMPCASES" == "" ]; then
cat << EOF >> $scriptfile
echo "    Input file: $in"
EOF
else
cat << EOF >> $scriptfile
echo "    Input files: "
EOF
for i in `seq 1 $OPENMPCASES`; do
cat << EOF >> $scriptfile
echo "       ${files[$i]}"
EOF
done
fi
cat << EOF >> $scriptfile
echo "     Directory: \`pwd\`"
echo "          Host: \`hostname\`"
echo "----------------" >> $qlog
echo "started running at \`date\`" >> $qlog
EOF
if [ "$OPENMPCASES" == "" ]; then
if [ "$vtuneresdir" == "" ]; then
if [ "$iinspectresdir" == "" ]; then
cat << EOF >> $scriptfile
$MPIRUN $exe $in $OUT2ERROR
EOF
else
cat << EOF >> $scriptfile
$MPIRUN $iinspectargs $exe $in $OUT2ERROR
EOF
fi
else
cat << EOF >> $scriptfile
$MPIRUN $vtuneargs $exe $in $OUT2ERROR
EOF
fi
else
for i in `seq 1 $OPENMPCASES`; do
if [ "$vtuneresdir" == "" ]; then
if [ "$iinspectresdir" == "" ]; then
cat << EOF >> $scriptfile

export OMP_NUM_THREADS=${nthreads[$i]}
$MPIRUN $exe ${files[$i]} $OUT2ERROR
EOF
else
cat << EOF >> $scriptfile

export OMP_NUM_THREADS=${nthreads[$i]}
$MPIRUN $iinspectargs $exe ${files[$i]} $OUT2ERROR
EOF
fi
else
cat << EOF >> $scriptfile

export OMP_NUM_THREADS=${nthreads[$i]}
$MPIRUN $vtuneargs $exe ${files[$i]} $OUT2ERROR
EOF
fi
done
fi
cat << EOF >> $scriptfile
echo "finished running at \`date\`" >> $qlog
EOF
if [ "$queue" == "none" ]; then
cat << EOF >> $scriptfile
rm -f $scriptfile
EOF
fi

#*** output script file to screen if -v option was selected

if [ "$showinput" == "1" ]; then
  cat $scriptfile
  echo
  exit
fi

#*** output info to screen
echo "submitted at `date`"                          > $qlog
if [ "$queue" != "none" ]; then
if [ "$OPENMPCASES" == "" ]; then
  echo "         Input file:$in"             | tee -a $qlog
else
  echo "         Input files:"               | tee -a $qlog
for i in `seq 1 $OPENMPCASES`; do
  echo "            ${files[$i]}"            | tee -a $qlog
done
fi
  echo "         Executable:$exe"            | tee -a $qlog
  if [ "$OPENMPI_PATH" != "" ]; then
    echo "            OpenMPI:$OPENMPI_PATH" | tee -a $qlog
  fi
  if [ "$use_intel_mpi" != "" ]; then
    echo "           Intel MPI"              | tee -a $qlog
  fi

#*** output currently loaded modules and modules when fds was built if the
#    1) -C option was selected and
#    2) currently loaded modules and fds loaded modules are diffent

  if [ "$FDS_MODULE_OPTION" == "" ]; then
    if [[ "$FDS_LOADED_MODULES" != "" ]] && [[ "$CURRENT_LOADED_MODULES" != "" ]]; then
      if [ "$FDS_LOADED_MODULES" != "$CURRENT_LOADED_MODULES" ]; then
        echo "  Modules(when run):$CURRENT_LOADED_MODULES" | tee -a $qlog
        echo "Modules(when built):$FDS_LOADED_MODULES"     | tee -a $qlog
        MODULES_OUT=1
      fi
    fi
  fi

#*** otherwise output modules used when fds is run
  if [[ "$MODULES" != "" ]] && [[ "$MODULES_OUT" == "" ]]; then
    echo "            Modules:$MODULES"                    | tee -a $qlog
  fi
  echo "   Resource Manager:$RESOURCE_MANAGER"             | tee -a $qlog
  echo "              Queue:$queue"                        | tee -a $qlog
  echo "              Nodes:$nodes"                        | tee -a $qlog
  echo "          Processes:$n_mpi_processes"              | tee -a $qlog
  echo " Processes per node:$n_mpi_processes_per_node"     | tee -a $qlog
  if test $n_openmp_threads -gt 1 ; then
    echo "Threads per process:$n_openmp_threads"           | tee -a $qlog
  fi
fi

#*** run script

$SLEEP
echo 
chmod +x $scriptfile
if [ "$queue" != "none" ]; then
  $QSUB $scriptfile | tee -a $qlog
else
  $QSUB $scriptfile
fi
if [ "$queue" != "none" ]; then
  cat $scriptfile > $scriptlog
  echo "#$QSUB $scriptfile" >> $scriptlog
  rm $scriptfile
fi
