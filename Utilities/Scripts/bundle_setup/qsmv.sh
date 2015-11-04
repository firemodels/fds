#!/bin/bash
# $Date$ 
# $Revision$
# $Author$
#
PROG=$0

# setup default queue name

progname=qsmv.sh
queue=batch
haltjob=0

nthreads=1

if [ $# -lt 1 ]
then
  echo "Usage: $progname [-d directory] [-r] [-f repository root] [-q queue] [smokeview_command] casename"
  echo ""
  echo "This script runs Smokeiew using an executable"
  echo "specified on the command line or from the respository if -r is specified."
  echo "Alternate queues (vis, fire70s) are set using the -q option."
  echo ""
  echo " -b use debug version"
  echo " -d directory [default: .]"
  echo " -q queue - name of the queue. choices: [default: $queue (other choices:"  
  echo "    vis and fire70s)"
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
dir=.
# parameters used by Smokeview
VOLRENDER=
SKIPFRAME=1
STARTFRAME=0
exe2=
ABORTRUN=n
DB=
SCRIPTFILE=

# read in parameters from command line

while getopts 'bd:f:hm:n:o:p:q:rsxy:z:t' OPTION
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
  q)
   queue="$OPTARG"
   ;;
  r)
   use_repository=1
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

use_repository=1
if [ "$queue" == "batch" ] ; then
  queue=fire70s
fi
if [ "$SCRIPTFILE" != "" ] ; then
  SCRIPTFILE = "-m $SCRIPTFILE"
fi

if [ "$use_debug" == "1" ] ; then
DB=_db
fi

if [ $use_repository -eq 0 ]
then
#set fds and the input file using the command line
  exe=$1
  in=$2
else
# for now only one instance of smokeview can occur per node
    exe="$FDSROOT/Verification/scripts/runsmv_single.sh"
    exe2="-x -y $STARTFRAME -z $SKIPFRAME $SCRIPTFILE"
 in=$1
fi

infile=${in%.*}

# define options used by smokeview (for computing volume rendering smoke frames) 

if [[ "$VOLRENDER" == "y" ]] ; then
  VOLRENDER="-volrender"
  STARTFRAME="-startframe $STARTFRAME"
  SKIPFRAME="-skipframe $SKIPFRAME"
fi

# if there is more than 1 process then use the mpirun command
#  (which will never happen if smokeview is running)

TITLE="$infile"

TITLE="$infile(SMV)"

cd $dir
fulldir=`pwd`

out=$fulldir/$infile.err
outlog=$fulldir/$infile.log
stopfile=$fulldir/$infile.stop


in_full_file=$fulldir/$infile.smv

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

QSUB="qsub -q $queue"
if [ "$queue" == "terminal" ] ; then
  QSUB=
  MPIRUN=
fi

scriptfile=/tmp/script.$$
cat << EOF > $scriptfile
#!/bin/bash
#PBS -N $TITLE
#PBS -e $out
#PBS -o $outlog
#PBS -l nodes=1:ppn=1
#\$ -N $TITLE
#\$ -e $out
#\$ -o $outlog
#\$ -l nodes=1:ppn=1

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
echo "         Processes:1"
chmod +x $scriptfile
$QSUB $scriptfile
rm $scriptfile
