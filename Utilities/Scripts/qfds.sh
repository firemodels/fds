#!/bin/bash

#*** environment varables

# OMP_PLACES       - cores, sockets or threads
# OMP_PROC_BIND    - false, true, master, close or spread
# RESOURCE_MANAGER - SLURM or TORQUE (default TORQUE)

#*** environment variables used by the bots
# JOBPREFIX        - prefix job title with $JOBPREFIX eg. FB_ or  SB_
# BACKGROUND_PROG  - defines location of background program
#                    ( if the 'none' queue is also specified)
# SCRIPTFILES      - outputs the name of the script file to $SCRIPTFILES
#                    ( used to kill jobs )

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
  echo "Usage: qfds.sh [-p nmpi_processes] [-o nthreads] [-e fds_command] [-q queue]  casename.fds"
  echo ""
  echo "qfds.sh runs FDS using an executable from the repository or one specified with the -e option."
  echo "A parallel version of FDS is invoked by using -p to specify the number of MPI processes and/or"
  echo "-o to specify the number of OpenMP threads."
  echo ""
  echo "qfds.sh loads the modules that were loaded when fds was built unless the -C option is specified"
  echo "then the currently loaded modules are used."
  echo ""
  echo " -e exe - full path of FDS used to run case "
  echo "    [default: $FDSROOT/fds/Build/mpi_intel_${platform}_64$DB/fds_mpi_intel_${platform}_64$DB]"
  echo " -h   - show commony used options"
  echo " -H   - show all options"
  echo " -o o - number of OpenMP threads per process [default: 1]"
  echo " -p p - number of MPI processes [default: 1] "
  echo " -q q - name of queue. [default: batch]"
  echo "        If q is terminal then casename.fds is run in the foreground on the local computer"
  echo " -v   - output generated script to standard output"
  echo "input_file - input file"
  if [ "$HELP" == "" ]; then
    exit
  fi
  echo "Other options:"
  echo " -C   - use modules currently loaded rather than modules loaded when fds was built."
  echo " -d dir - specify directory where the case is found [default: .]"
  echo " -E - use tcp transport (only available with the Intel compiled versions of fds)"
  echo "      This options adds export I_MPI_FABRICS=shm:tcp to the run script"
  echo " -f repository root - name and location of repository where FDS is located"
  echo "    [default: $FDSROOT]"
  echo " -i use installed fds"
  echo " -I use Intel mpi version of fds"
  echo " -m m - reserve m processes per node [default: 1]"
  echo " -M   -  add --mca plm_rsh_agent /usr/bin/ssh to mpirun command "
  echo " -n n - number of MPI processes per node [default: 1]"
  echo " -N   - do not use socket or report binding options"
  echo " -O n - run cases casea.fds, caseb.fds, ... using 1, ..., N OpenMP threads"
  echo "        where case is specified on the command line. N can be at most 9."
  echo " -r   - report bindings"
  echo " -s   - stop job"
  echo " -S   - use startup files to set the environment, do not load modules"
  echo " -t   - used for timing studies, run a job alone on a node (reserving $NCORES_COMPUTENODE cores)"
  echo " -T type - run dv (development) or db (debug) version of fds"
  echo "           if -T is not specified then the release version of fds is used"
  echo " -w time - walltime, where time is hh:mm for PBS and dd-hh:mm:ss for SLURM. [default: $walltime]"
  echo ""
  exit
}

#*** get directory containing qfds.sh

QFDS_PATH=$(dirname `which $0`)
CURDIR=`pwd`
cd $QFDS_PATH
QFDS_DIR=`pwd`
cd $CURDIR

#*** define toplevel of the repos

FDSROOT=~/FDS-SMV
if [ "$FIREMODELS" != "" ]; then
  FDSROOT=$FIREMODELS
fi

#*** define resource manager that is used

if [ "$RESOURCE_MANAGER" == "SLURM" ]; then
  if [ "$SLURM_MEM" != "" ]; then
    SLURM_MEM="#SBATCH --mem=$SLURM_MEM"
  fi
  if [ "$SLURM_MEMPERCPU" != "" ]; then
    SLURM_MEM="#SBATCH --mem-per-cpu=$SLURM_MEMPERCPU"
  fi
else
  RESOURCE_MANAGER="TORQUE"
fi

#*** determine platform

platform="linux"
if [ "`uname`" == "Darwin" ] ; then
  platform="osx"
fi

#*** determine number of cores and default queue

if [ "$platform" == "osx" ]; then
  queue=none
  ncores=`system_profiler SPHardwareDataType|grep Cores|awk -F' ' '{print $5}'`
else
  queue=batch
  ncores=`grep processor /proc/cpuinfo | wc -l`
fi
if [ "$NCORES_COMPUTENODE" == "" ]; then
  NCORES_COMPUTENODE=$ncores
else
  ncores=$NCORES_COMPUTENODE
fi

#*** set default parameter values

OMPPLACES=$OMP_PLACES
OMPPROCBIND=$OMP_PROCBIND
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

nmpi_processes=1
nmpi_processes_per_node=-1
max_processes_per_node=1
nopenmp_threads=1
use_installed=
use_debug=
use_devel=
use_intel_mpi=
dir=.
benchmark=no
showinput=0
REPORT_BINDINGS="--report-bindings"
nosocket=
exe=
STARTUP=
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

#*** read in parameters from command line

while getopts 'ACd:e:Ef:hHiIm:MNn:o:O:p:Pq:rsStT:vw:' OPTION
do
case $OPTION  in
  A) # used by timing scripts to identify benchmark cases
   DUMMY=1
   ;;
  C)
   FDS_MODULE_OPTION=
   ;;
  d)
   dir="$OPTARG"
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
   nosocket="1"
   ;;
  M)
   MCA="--mca plm_rsh_agent /usr/bin/ssh "
   ;;
  m)
   max_processes_per_node="$OPTARG"
   ;;
  n)
   nmpi_processes_per_node="$OPTARG"
   ;;
  N)
   nosocket="1"
   ;;
  o)
   nopenmp_threads="$OPTARG"
   ;;
  O)
   OPENMPCASES="$OPTARG"
   if [ $OPENMPCASES -gt 9 ]; then
     OPENMPCASES=9
   fi
   nmpi_process=1
   benchmark="yes"
   if [ "$NCORES_COMPUTENODE" != "" ]; then
     nmpi_processes_per_node="$NCORES_COMPUTENODE"
   fi
   ;;
  p)
   nmpi_processes="$OPTARG"
   ;;
  P)
   OPENMPCASES="2"
   OPENMPTEST="1"
   benchmark="yes"
   nmpi_process=1
   if [ "$NCORES_COMPUTENODE" != "" ]; then
     nmpi_processes_per_node="$NCORES_COMPUTENODE"
   fi
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
  S)
   STARTUP=1
   ;;
  t)
   benchmark="yes"
   if [ "$NCORES_COMPUTENODE" != "" ]; then
     nmpi_processes_per_node="$NCORES_COMPUTENODE"
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
   ;;
  v)
   showinput=1
   ;;
  w)
   walltime="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

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
  nthreads[1]=$nopenmp_threads
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

if [[ "$OMPPLACES" != "" ]]; then
  if [[ "$OMPPLACES" != "cores" ]] &&  [[ "$OMPPLACES" != "sockets" ]] &&  [[ "$OMPPLACES" == "threads" ]]; then
    echo "*** error: OMP_PLACES can only be cores, sockets or threads"
    exit
  fi
  OMPPLACES="OMP_PLACES=$OMPPLACES"
fi

if [ "$OMPPROCBIND" != "" ]; then
  if [[ "$OMPPROCBIND" != "false" ]] &&  [[ "$OMPPROCBIND" != "true" ]] &&  [[ "$OMPPROCBIND" != "master" ]] &&  [[ "$OMPPROCBIND" == "close" ]] &&  [[ "$OMPPROCBIND" == "spread" ]]; then
    echo "*** error: OMP_PROCBIND can only be false, true, master, close or spread"
    exit
  fi
  OMPPROCBIND="OMP_PROC_BIND=$OMPPROCBIND"
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
  if [ "$use_mpi_intel" == "" ]; then
    is_intel_mpi=`echo "" | $exe 2>&1 >/dev/null | grep MPI | grep library | grep Intel | wc -l`
    if [ "$is_intel_mpi" == "1" ]; then
         use_intel_mpi=1
         nosocket="1"
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

if test $nmpi_processes_per_node -gt $ncores ; then
  nmpi_processes_per_node=$ncores
fi

if test $nmpi_processes_per_node = -1 ; then
  if test $nmpi_processes -gt 1 ; then
    nmpi_processes_per_node=2
  else
    nmpi_processes_per_node=1
  fi
fi

let "nodes=($nmpi_processes-1)/$nmpi_processes_per_node+1"
if test $nodes -lt 1 ; then
  nodes=1
fi
if [ "$RESOURCE_MANAGER" == "SLURM" ]; then
    nodes=""
fi

#*** define processes per node

let "ppn=($nmpi_processes_per_node)"
if test $ppn -le $max_processes_per_node ; then
  ppn=$max_processes_per_node
fi

# default: Use mpirun option to bind processes to socket (for MPI).
# Or, bind processs to and map processes by socket if
# OpenMP is being used (number of OpenMP threads > 1).

if test $nmpi_processes -gt 1 ; then
 if test $nopenmp_threads -gt 1 ; then
  SOCKET_OPTION="--map-by socket:PE=$nopenmp_threads"
 else
  SOCKET_OPTION=" "
 fi
else
 SOCKET_OPTION="--map-by node:PE=$nopenmp_threads"
fi

#*** the "none" queue does not use the queing system,
#    so blank out SOCKET_OPTIONS and REPORT_BINDINGS

if [[ "$queue" == "none" ]] || [[ "$nosocket" == "1" ]]; then
 SOCKET_OPTION=
 REPORT_BINDINGS=
fi

#*** define MPIRUNEXE and do some error checking

if [ "$use_intel_mpi" == "1" ]; then # using Intel MPI
  if [ "$use_installed" == "1" ]; then
    MPIRUNEXE=$fdsdir/mpiexec
    if [ ! -e $MPIRUNEXE ]; then
      echo "$MPIRUNEXE not found"
      echo "Run aborted"
      ABORT=y
    fi
    MPILABEL="IMPI"
  else
    if [ "$I_MPI_ROOT" == "" ]; then
      echo "Intel MPI environment not setup. Run aborted."
      ABORTRUN=y
    else
      MPIRUNEXE=$I_MPI_ROOT/bin64/mpiexec
      if [ ! -e $MPIRUNEXE ]; then
        echo "Intel mpiexec, $MPIRUNEXE, not found at:"
        echo "$MPIRUNEXE"
        ABORTRUN=y
        echo "Run aborted"
      fi
      MPILABEL="IMPI"
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
  MPILABEL="MPI"
fi

TITLE="$infile($MPILABEL)"
MPIRUN="$MPIRUNEXE $REPORT_BINDINGS $SOCKET_OPTION $MCA -np $nmpi_processes"

cd $dir
fulldir=`pwd`

#*** define files

outerr=$fulldir/$infile.err
outlog=$fulldir/$infile.log
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

#QSUB="qsub -k eo -q $queue"
QSUB="qsub -q $queue"

if [ "$queue" == "terminal" ]; then
  QSUB=
  MPIRUN=
fi

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

#*** setup for SLURM (alternative to torque)

  if [ "$RESOURCE_MANAGER" == "SLURM" ]; then
    QSUB="sbatch -p $queue --ignore-pbs"
    MPIRUN='srun'
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
EOF

if [ "$queue" != "none" ]; then
  if [ "$RESOURCE_MANAGER" == "SLURM" ]; then
    cat << EOF >> $scriptfile
#SBATCH -J $JOBPREFIX$infile
#SBATCH -e $outerr
#SBATCH -o $outlog
#SBATCH -p $queue
#SBATCH -n $nmpi_processes
####SBATCH --nodes=$nodes
#SBATCH --cpus-per-task=$nopenmp_threads
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
export OMP_NUM_THREADS=$nopenmp_threads
EOF
fi

if [ "$use_intel_mpi" == "1" ]; then
cat << EOF >> $scriptfile
export I_MPI_DEBUG=5
EOF
fi

if [ "$TCP" != "" ]; then
  cat << EOF >> $scriptfile
export I_MPI_FABRICS=shm:tcp
EOF
fi

if test $nopenmp_threads -gt 1 ; then
  if [ "$OMPPLACES" != "" ]; then
    cat << EOF >> $scriptfile
export $OMPPLACES
EOF
  fi

  if [ "$OMPPROCBIND" != "" ]; then
    cat << EOF >> $scriptfile
export $OMPPROCBIND
EOF
  fi
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
EOF
if [ "$OPENMPCASES" == "" ]; then
cat << EOF >> $scriptfile
$MPIRUN $exe $in $OUT2ERROR
EOF
else
for i in `seq 1 $OPENMPCASES`; do
cat << EOF >> $scriptfile

export OMP_NUM_THREADS=${nthreads[$i]}
$MPIRUN $exe ${files[$i]} $OUT2ERROR
EOF
done
fi
if [ "$queue" == "none" ]; then
cat << EOF >> $scriptfile
rm -f $scriptfile
EOF
fi

#*** output script file to screen if -v option was selected

if [ "$showinput" == "1" ]; then
  cat $scriptfile
  exit
fi

#*** output info to screen

if [ "$queue" != "none" ]; then
if [ "$OPENMPCASES" == "" ]; then
  echo "         Input file:$in"
else
  echo "         Input files:"
for i in `seq 1 $OPENMPCASES`; do
  echo "            ${files[$i]}"
done
fi
  echo "         Executable:$exe"
  if [ "$OPENMPI_PATH" != "" ]; then
    echo "            OpenMPI:$OPENMPI_PATH"
  fi
  if [ "$use_intel_mpi" != "" ]; then
    echo "           Intel MPI"
  fi

#*** output currently loaded modules and modules when fds was built if the
#    1) -C option was selected and
#    2) currently loaded modules and fds loaded modules are diffent

  if [ "$FDS_MODULE_OPTION" == "" ]; then
    if [[ "$FDS_LOADED_MODULES" != "" ]] && [[ "$CURRENT_LOADED_MODULES" != "" ]]; then
      if [ "$FDS_LOADED_MODULES" != "$CURRENT_LOADED_MODULES" ]; then
        echo "  Modules(when run):$CURRENT_LOADED_MODULES"
        echo "Modules(when built):$FDS_LOADED_MODULES"
        MODULES_OUT=1
      fi
    fi
  fi

#*** otherwise output modules used when fds is run

  if [[ "$MODULES" != "" ]] && [[ "$MODULES_OUT" == "" ]]; then
    echo "            Modules:$MODULES"
  fi
  echo "              Queue:$queue"
  echo "              Nodes:$nodes"
  echo "          Processes:$nmpi_processes"
  echo " Processes per node:$nmpi_processes_per_node"
  if test $nopenmp_threads -gt 1 ; then
    echo "Threads per process:$nopenmp_threads"
  fi
fi

#*** run script

chmod +x $scriptfile
if [ "$SCRIPTFILES" != "" ]; then
  echo $(basename "$scriptfile") >> $SCRIPTFILES
fi
$QSUB $scriptfile
if [ "$queue" != "none" ]; then
  cat $scriptfile > $scriptlog
  echo "#$QSUB $scriptfile" >> $scriptlog
  rm $scriptfile
fi
