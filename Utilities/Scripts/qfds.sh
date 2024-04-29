#!/bin/bash

# ---------------------------- stop_fds_if_requested ----------------------------------

function stop_fds_if_requested {
if [ "$STOPFDS" != "" ]; then
  echo "stopping case: $in"
 touch $stopfile $stopcatfile
 exit
fi

if [ "$STOPFDSMAXITER" != "" ]; then
  echo "creating delayed stop file: $infile"
  echo $STOPFDSMAXITER > $stopfile
  echo $STOPFDSMAXITER > $stopcatfile
fi

if [ "$stopjob" == "1" ]; then
  echo "stopping case: $in"
  touch $stopfile $stopcatfile
  exit
fi

if [ "$STOPFDSMAXITER" == "" ]; then
  rm -f $stopfile $stopcatfile
fi
}

# ---------------------------- usage ----------------------------------

function usage {
  if [ "$use_intel_mpi" == "1" ]; then
    MPI=impi
  else
    MPI=ompi
  fi
  echo "Usage: qfds.sh [-p n_mpi_processes] [-o nthreads] [-e fds_command] [-q queue]  casename.fds"
  echo ""
  echo "qfds.sh runs FDS using an executable from the repository or one specified with the -e option."
  echo ""
  echo " -e exe - full path of FDS used to run case "
  echo "    [default: $FDSROOT/fds/Build/${MPI}_intel_linux$DB/fds_${MPI}_intel_linux$DB]"
  echo " -h   - show commonly used options"
  echo " -H   - show all options"
  echo " -o o - number of OpenMP threads per process [default: 1]"
  echo " -p p - number of MPI processes [default: 1] "
  echo " -q q - name of queue. [default: batch]"
  echo " -v   - output generated script to standard output"
  if [ "$HELP" == "" ]; then
    return
  fi
  echo "Other options:"
  echo " -b email_address - send an email to email_address when jobs starts, aborts and finishes"
  echo " -d dir - specify directory where the case is found [default: .]"
  echo " -g   - only run if input file and executable are not dirty"
  echo " -I use Intel MPI version of fds"
  echo " -j prefix - specify a job prefix"
  echo " -L use Open MPI version of fds"
  echo " -n n - number of MPI processes per node [default: 1]"
  echo " -O n - run cases casea.fds, caseb.fds, ... using 1, ..., N OpenMP threads"
  echo "        where case is specified on the command line. N can be at most 9."
  echo " -s   - stop job"
  echo " -t   - used for timing studies, run a job alone on a node (reserving $NCORES_COMPUTENODE cores)"
  echo " -T type - run dv (development) or db (debug) version of fds"
  echo "           if -T is not specified then the release version of fds is used"
  echo " -U n - only allow n jobs owned by `whoami` to run at a time"
  echo " -w time - maximum run time, where time is dd-hh:mm:ss [default: $walltime]"
  echo " -y dir - run case in directory dir"
  echo " -Y   - run case in directory casename where casename.fds is the case being run"
  echo ""
  echo " Resource manager: $RESOURCE_MANAGER"
}

#*** get directory containing qfds.sh

QFDS_PATH=$(dirname `which $0`)
CURDIR=`pwd`
cd $QFDS_PATH
QFDS_DIR=`pwd`

#*** define toplevel of the repos

FDSROOT=$QFDS_DIR/../../..
cd $FDSROOT
FDSROOT=`pwd`
cd $CURDIR

#*** qfds.sh only works on a linux OS

if [ "$WINDIR" != "" ]; then
  echo "***Error: only Linux platforms are supported"
  exit
fi
if [ "`uname`" == "Darwin" ] ; then
  echo "***Error: only Linux platforms are supported"
  exit
fi

#*** determine number of cores and default queue

queue=batch
if [ "$QFDS_NCORES" == "" ]; then
  ncores=`grep processor /proc/cpuinfo | wc -l`
else
  ncores=$QFDS_NCORES
fi
if [ "$NCORES_COMPUTENODE" == "" ]; then
  NCORES_COMPUTENODE=$ncores
else
  ncores=$NCORES_COMPUTENODE
fi

#*** set default parameter values

HELP=
MPIRUN=
ABORTRUN=n
DB=
OUT2ERROR=
stopjob=0

n_mpi_processes=1
max_mpi_processes_per_node=1000
n_openmp_threads=1
use_debug=
use_devel=
use_intel_mpi=1
EMAIL=
CHECK_DIRTY=
casedir=
use_default_casedir=
USERMAX=

# determine which resource manager is running (or none)

STATUS_FILE=status_file.$$
srun -V >& $STATUS_FILE
missing_slurm=`cat $STATUS_FILE | tail -1 | grep "not found" | wc -l`
rm -f $STATUS_FILE

RESOURCE_MANAGER="NONE"
if [ $missing_slurm -eq 0 ]; then
  RESOURCE_MANAGER="SLURM"
else
  echo "***error: The slurm resource manager was not found and is required."
  exit
fi
if [ "$SLURM_MEM" != "" ]; then
 SLURM_MEM="#SBATCH --mem=$SLURM_MEM"
fi
if [ "$SLURM_MEMPERCPU" != "" ]; then
 SLURM_MEM="#SBATCH --mem-per-cpu=$SLURM_MEMPERCPU"
fi

dir=.
benchmark=no
showinput=0
exe=
walltime=99-99:99:99

if [ $# -lt 1 ]; then
  usage
  exit
fi

commandline=`echo $* | sed 's/-V//' | sed 's/-v//'`

#*** read in parameters from command line

while getopts 'b:d:e:ghHIj:Ln:o:p:q:stT:U:vw:y:Y' OPTION
do
case $OPTION  in
  b)
   EMAIL="$OPTARG"
   ;;
  d)
   dir="$OPTARG"
   ;;
  e)
   exe="$OPTARG"
   ;;
  g)
   CHECK_DIRTY=1
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
  I)
   use_intel_mpi=1
   ;;
  j)
   JOBPREFIX="$OPTARG"
   ;;
  L)
   use_intel_mpi=
   ;;
  n)
   max_mpi_processes_per_node="$OPTARG"
   ;;
  o)
   n_openmp_threads="$OPTARG"
   ;;
  p)
   n_mpi_processes="$OPTARG"
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
  U)
   USERMAX="$OPTARG"
   ;;
  v)
   showinput=1
   ;;
  w)
   walltime="$OPTARG"
   ;;
  y)
   casedir="$OPTARG"
   ;;
  Y)
   use_default_casedir=1
   ;;
  *)
   usage
   exit 1
   ;;
esac
done
shift $(($OPTIND-1))

#*** define input file

in=$1
infile=${in%.*}

#*** run case in a sub-directory

if [ "$use_default_casedir" != "" ]; then
  casedir=$infile
fi
if [ "$casedir" != "" ]; then
  if [ ! -d $casedir ]; then
    mkdir $casedir
  fi
  cp $in $casedir/.
  cd $casedir
fi

#*** define executable

if [ "$use_debug" == "1" ]; then
  DB=_db
fi
if [ "$use_devel" == "1" ]; then
  DB=_dv
fi
if [ "$use_intel_mpi" == "1" ]; then
  if [ "$exe" == "" ]; then
    exe=$FDSROOT/fds/Build/impi_intel_linux$DB/fds_impi_intel_linux$DB
  fi
  if [[ $n_openmp_threads > 1 ]]; then
    exe=$FDSROOT/fds/Build/impi_intel_linux_openmp$DB/fds_impi_intel_linux_openmp$DB
  fi
fi
if [ "$exe" == "" ]; then
  exe=$FDSROOT/fds/Build/ompi_intel_linux$DB/fds_ompi_intel_linux$DB
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

CURRENT_LOADED_MODULES=`echo $LOADEDMODULES | tr ':' ' '`

MODULES=$CURRENT_LOADED_MODULES

#*** define number of nodes

let "n_mpi_processes_per_node=($ncores)/($n_openmp_threads)"
if [ $n_mpi_processes_per_node -gt $max_mpi_processes_per_node ]; then
  n_mpi_processes_per_node=$max_mpi_processes_per_node
fi
let nodes="($n_mpi_processes+$n_mpi_processes_per_node-1)/$n_mpi_processes_per_node"

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
MPIRUN="$MPIRUNEXE $REPORT_BINDINGS $SOCKET_OPTION -np $n_mpi_processes"

cd $dir
fulldir=`pwd`

#*** check if exe and/or input file is dirty before running

if [[ "$CHECK_DIRTY" == "1" ]] && [[ "$exe" != "" ]]; then
  if [ -e $exe ]; then
    is_dirty_exe=`echo "" | $exe |& grep dirty |& wc -l`
    dirty_exe=`   echo "" | $exe |& grep dirty |& awk '{print $3}'`
    is_dirty_input=`git diff $in   |& wc -l`

    is_dirty=
    if [ $is_dirty_exe -gt 0 ]; then
      is_dirty=1
    fi
    if [ $is_dirty_input -gt 0 ]; then
      is_dirty=1
    fi

    if [ "$is_dirty" == "1" ]; then
      echo ""
      if [ $is_dirty_exe -gt 0 ]; then
        echo "***error: source used to build FDS is dirty."
      fi
      echo "executable: $exe"
      echo "          $dirty_exe"
      if [ $is_dirty_input -gt 0 ]; then
        echo "***error: input file $in is dirty."
      else
        echo "input file: $in"
      fi
    fi
    if [ "$is_dirty" == "1" ]; then
      echo "Use the -g option to ignore this error"
      echo "Exiting."
      exit 1
    fi
  fi
fi

#*** define files

outerr=$fulldir/$infile.err
outlog=$fulldir/$infile.log
qlog=$fulldir/$infile.qlog
stopfile=$fulldir/$infile.stop
stopcatfile=$fulldir/${infile}_cat.stop
scriptlog=$fulldir/$infile.slog
in_full_file=$fulldir/$in

#*** make sure various files exist before running the case

if ! [ -e $in_full_file ]; then
  if [ "$showinput" == "0" ]; then
    echo "The input file, $in_full_file, does not exist. Run aborted."
    ABORTRUN=y
  fi
fi

if [ "$STOPFDS" == "" ]; then
  if [ "$exe" != "" ]; then
    if ! [ -e $exe ]; then
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

#*** setup for SLURM

QSUB="sbatch -p $queue"
if [ "$use_intel_mpi" == "1" ]; then
   MPIRUN="srun --mpi=pmi2"
else
   MPIRUN="srun "
fi

#*** Set walltime parameter only if walltime is specified as input argument

walltimestring_slurm=
if [ "$walltime" != "" ]; then
  walltimestring_slurm="--time=$walltime"
fi

#*** create a random script filename for submitting jobs

scriptfile=`mktemp /tmp/script.$$.XXXXXX`

cat << EOF > $scriptfile
#!/bin/bash
# $0 $commandline
EOF

cat << EOF >> $scriptfile
#SBATCH -J $JOBPREFIX$infile
#SBATCH -e $outerr
#SBATCH -o $outlog
#SBATCH --partition=$queue
#SBATCH --ntasks=$n_mpi_processes
#SBATCH --nodes=$nodes
#SBATCH --cpus-per-task=$n_openmp_threads
#SBATCH --ntasks-per-node=$n_mpi_processes_per_node
EOF
if [ "$EMAIL" != "" ]; then
    cat << EOF >> $scriptfile
#SBATCH --mail-user=$EMAIL
#SBATCH --mail-type=ALL
EOF
fi

if [ "$benchmark" == "yes" ]; then
cat << EOF >> $scriptfile
#SBATCH --exclusive
#SBATCH --cpu-freq=Performance
EOF
fi

if [ "$SLURM_MEM" != "" ]; then
cat << EOF >> $scriptfile
$SLURM_MEM
EOF
fi

if [ "$SLURM_PSM" != "" ]; then
cat << EOF >> $scriptfile
$SLURM_PSM
EOF
fi

if [ "$walltimestring_slurm" != "" ]; then
      cat << EOF >> $scriptfile
#SBATCH $walltimestring_slurm

EOF
fi

cat << EOF >> $scriptfile
export OMP_NUM_THREADS=$n_openmp_threads
EOF

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

if [ "$PROVIDER" != "" ]; then
cat << EOF >> $scriptfile
$PROVIDER
EOF
fi

cat << EOF >> $scriptfile

cd $fulldir
echo
echo \`date\`
EOF

cat << EOF >> $scriptfile
echo "    Input file: $in"
EOF
if [ "$casedir" != "" ]; then
cat << EOF >> $scriptfile
echo "     Input dir: $casedir"
EOF
fi

cat << EOF >> $scriptfile
echo "     Directory: \`pwd\`"
echo "          Host: \`hostname\`"
echo "----------------" >> $qlog
echo "started running at \`date\`" >> $qlog
EOF

cat << EOF >> $scriptfile
$MPIRUN $exe $in $OUT2ERROR
EOF

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

# wait until number of jobs running alread by user is less than USERMAX
if [ "$USERMAX" != "" ]; then
  nuser=`squeue | grep -v JOBID | awk '{print $4}' | grep $USER | wc -l`
  while [ $nuser -gt $USERMAX ]
  do
    nuser=`squeue | grep -v JOBID | awk '{print $4}' | grep $USER | wc -l`
    sleep 10
  done
fi

#*** output info to screen
echo "submitted at `date`"                          > $qlog
if [ "$queue" != "none" ]; then
  echo "         Input file:$in"             | tee -a $qlog
if [ "$casedir" != "" ]; then
  echo "          Input dir:$casedir"             | tee -a $qlog
fi

  echo "         Executable:$exe"            | tee -a $qlog
  if [ "$OPENMPI_PATH" != "" ]; then
    echo "            OpenMPI:$OPENMPI_PATH" | tee -a $qlog
  fi
  if [ "$use_intel_mpi" != "" ]; then
    echo "           Intel MPI"              | tee -a $qlog
  fi

#*** output modules used when fds is run
  if [[ "$MODULES" != "" ]] && [[ "$MODULES_OUT" == "" ]]; then
    echo "            Modules:$MODULES"                    | tee -a $qlog
  fi
  echo "   Resource Manager:$RESOURCE_MANAGER"             | tee -a $qlog
  echo "              Queue:$queue"                        | tee -a $qlog
  echo "              Nodes:$nodes"                        | tee -a $qlog
  echo "          Processes:$n_mpi_processes"              | tee -a $qlog
  if test $n_openmp_threads -gt 1 ; then
    echo "Threads per process:$n_openmp_threads"           | tee -a $qlog
  fi
fi

#*** run script

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
