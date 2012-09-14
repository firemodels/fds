#!/bin/bash -f

# defaults

queue=
background=no
QSUB=qsub

if [ "$JOBPREFIX" == "" ]; then
  JOBPREFIX=VV_
fi

while getopts 'q:' OPTION
do
case $OPTION in
  q)
   queue="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

# If queue is "none" then use "background" to submit jobs
# instead of qsub

if [ "$queue" == "none" ]; then
  queue=
  QSUB="$BACKGROUND -u 75 -d 10 "
  background=yes;
  if ! [ -e $BACKGROUND ];  then
    echo "The file $BACKGROUND does not exist. Run aborted"
    exit
  fi
fi
if [ "$queue" != "" ]; then
   queue="-q $queue"
fi

# setup parameters

scratchdir=$SVNROOT/Utilities/Scripts/tmp
dir=$1
infile=$2

fulldir=$BASEDIR/$dir
in=$infile
outerr=$fulldir/$infile.err
outlog=$fulldir/$infile.log

scriptfile=$scratchdir/script.$$

# Do not run CFAST if -s option is specified

if [ $STOPFDS ]; then
 echo "skipping CFAST case: $infile"
 exit
fi

# ensure that various files and directories exist

if ! [ -e $CFAST ];  then
  echo "The file $CFAST does not exist. Run aborted"
  exit
fi
if ! [ -d $fulldir ]; then
  echo "The directory $fulldir does not exist. Run aborted."
  exit
fi
if ! [ -e $fulldir/$in.in ]; then
  echo "The cfast input file, $fulldir/$in.in, does not exist. Run aborted."
  exit
fi
if [ -e $outlog ]; then
 rm $outlog
fi

# create run script

cat << EOF > $scriptfile
#!/bin/bash -f
#\$ -S /bin/bash
#\$ -N $JOBPREFIX$infile -e $outerr -o $outlog
#PBS -N $JOBPREFIX$infile -e $outerr -o $outlog
cd $fulldir

echo Time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`

$CFAST $in 
EOF

echo Running `basename $CFAST` $in 
if [ "$background" != "yes" ]; then
  chmod +x $scriptfile
  $QSUB $queue $scriptfile
else
  cd $fulldir
  $QSUB $CFAST $in
fi

rm $scriptfile
