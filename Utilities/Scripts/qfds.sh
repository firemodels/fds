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
  echo " -I use Intel MPI version of fds"
  echo " -j prefix - specify a job prefix"
  echo " -L use Open MPI version of fds"
  echo " -n n - number of MPI processes per node [default: 1]"
  echo " -P use PBS/Torque"
  echo " -s   - stop job"
  echo " -t   - used for timing studies, run a job alone on a node (reserving $NCORES_COMPUTENODE cores)"
  echo " -T type - run dv (development) or db (debug) version of fds"
  echo "           if -T is not specified then the release version of fds is used"
  echo " -U n - only allow n jobs owned by `whoami` to run at a time"
  echo " -w time - maximum run time: SLURM dd-hh:mm:ss ; PBS/Torque hh:mm:ss "
  echo " -y dir - run case in directory dir"
  echo " -Y   - run case in directory casename where casename.fds is the case being run"
  echo ""
  echo " Resource manager: $RESOURCE_MANAGER"
}

#*** Define the default resource manager

RESOURCE_MANAGER="SLURM"

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

#*** determine number of cores

if [ "$QFDS_NCORES" == "" ]; then
  n_cores=`grep processor /proc/cpuinfo | wc -l`
else
  n_cores=$QFDS_NCORES
fi
if [ "$NCORES_COMPUTENODE" == "" ]; then
  NCORES_COMPUTENODE=$n_cores
else
  n_cores=$NCORES_COMPUTENODE
fi

#*** set default parameter values

queue=batch
HELP=
ABORT_RUN=n
DB=
stopjob=0
n_mpi_processes=1
max_mpi_processes_per_node=1000
n_openmp_threads=1
use_debug=
use_devel=
use_intel_mpi=1
EMAIL=
casedir=
use_default_casedir=
USERMAX=
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

while getopts 'b:d:e:hHIj:Ln:o:Pp:q:stT:U:vw:y:Y' OPTION
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
  P)
   RESOURCE_MANAGER="PBS/Torque"
   if [ "$walltime" == "99-99:99:99" ]; then
     walltime=999:99:99
   fi
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

#*** define number of nodes

let "n_mpi_processes_per_node=($n_cores)/($n_openmp_threads)"
if [ $n_mpi_processes_per_node -gt $max_mpi_processes_per_node ]; then
  n_mpi_processes_per_node=$max_mpi_processes_per_node
fi
if [ $n_mpi_processes_per_node -gt $n_mpi_processes ]; then
  n_mpi_processes_per_node=$n_mpi_processes
fi
let nodes="($n_mpi_processes+$n_mpi_processes_per_node-1)/$n_mpi_processes_per_node"

cd $dir
fulldir=`pwd`

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
    ABORT_RUN=y
  fi
fi

if [ "$STOPFDS" == "" ]; then
  if [ "$exe" != "" ]; then
    if ! [ -e $exe ]; then
      if [ "$showinput" == "0" ]; then
        echo "The program, $exe, does not exist. Run aborted."
        ABORT_RUN=y
      fi
    fi
  fi

  if [ -e $outlog ]; then
    echo "Removing log file: $outlog"
    rm $outlog
  fi

  if [ "$ABORT_RUN" == "y" ]; then
    if [ "$showinput" == "0" ]; then
      exit
    fi
  fi
fi

stop_fds_if_requested

#*** define the execute command and the submit command

if [ "$RESOURCE_MANAGER" == "SLURM" ]; then
  QSUB="sbatch -p $queue"
  if [ "$USE_MPIRUN" == "" ]; then
# use the native Slurm process allocator
     MPIRUN="srun --mpi=pmi2 "
  else
# use on spark ( USE_MPIRUN variable is set to 1 in /etc/profile )
     MPIRUN="mpirun "
  fi
else
  QSUB="qsub -q $queue"
  MPIRUN="mpiexec -n $n_mpi_processes"
fi

#*** create a random script filename for submitting jobs

scriptfile=`mktemp /tmp/script.$$.XXXXXX`

cat << EOF > $scriptfile
#!/bin/bash
# $0 $commandline
EOF

if [ "$RESOURCE_MANAGER" == "SLURM" ]; then

cat << EOF >> $scriptfile
#SBATCH -J $JOBPREFIX$infile
#SBATCH -e $outerr
#SBATCH -o $outlog
#SBATCH --partition=$queue
#SBATCH --ntasks=$n_mpi_processes
#SBATCH --cpus-per-task=$n_openmp_threads
#SBATCH --time=$walltime
EOF

if [[ $n_openmp_threads -gt 1 ]] || [[ $max_mpi_processes_per_node -lt 1000 ]] ; then
cat << EOF >> $scriptfile
#SBATCH --nodes=$nodes
#SBATCH --ntasks-per-node=$n_mpi_processes_per_node
EOF
fi

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

else # PBS/Torque

cat << EOF >> $scriptfile
#PBS -N $JOBPREFIX$infile
#PBS -e $outerr
#PBS -o $outlog
#PBS -q $queue
#PBS -l nodes=$nodes:ppn=$n_mpi_processes
#PBS -l walltime=$walltime
EOF
if [ "$EMAIL" != "" ]; then
    cat << EOF >> $scriptfile
#PBS -M $EMAIL
#PBS -m abe
EOF
fi
if [ "$benchmark" == "yes" ]; then
cat << EOF >> $scriptfile
#PBS l naccesspolicy=singlejob
EOF
fi
fi

cat << EOF >> $scriptfile
export OMP_NUM_THREADS=$n_openmp_threads
EOF

if [ "$use_intel_mpi" == "1" ]; then
cat << EOF >> $scriptfile
export I_MPI_DEBUG=5
export I_MPI_PMI_VALUE_LENGTH_MAX=512
EOF
fi

if [[ $n_openmp_threads -gt 1 ]] && [[ "$use_intel_mpi" == "1" ]]; then
cat << EOF >> $scriptfile
export I_MPI_PIN_DOMAIN=omp
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
$MPIRUN $exe $in 
EOF

cat << EOF >> $scriptfile
echo "finished running at \`date\`" >> $qlog
EOF

#*** output script file to screen if -v option was selected

if [ "$showinput" == "1" ]; then
  cat $scriptfile
  echo
  exit
fi

#*** wait until number of jobs running already by user is less than USERMAX

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

echo "         Input file:$in"             | tee -a $qlog
if [ "$casedir" != "" ]; then
  echo "         Input dir:$casedir"             | tee -a $qlog
fi

echo "         Executable:$exe"            | tee -a $qlog
if [ "$use_intel_mpi" != "" ]; then
  echo "         Intel MPI"              | tee -a $qlog
fi

#*** output modules used when fds is run

CURRENT_LOADED_MODULES=`echo $LOADEDMODULES | tr ':' ' '`
MODULES=$CURRENT_LOADED_MODULES

if [ "$MODULES" != "" ]; then
    echo "            Modules:$MODULES"                    | tee -a $qlog
fi
echo "   Resource Manager:$RESOURCE_MANAGER"             | tee -a $qlog
echo "              Queue:$queue"                        | tee -a $qlog
echo "              Nodes:$nodes"                        | tee -a $qlog
echo "          Processes:$n_mpi_processes"              | tee -a $qlog
if test $n_openmp_threads -gt 1 ; then
  echo "Threads per process:$n_openmp_threads"           | tee -a $qlog
fi

#*** run script

echo
chmod +x $scriptfile

$QSUB $scriptfile | tee -a $qlog

cat $scriptfile > $scriptlog
echo "#$QSUB $scriptfile" >> $scriptlog
rm $scriptfile

