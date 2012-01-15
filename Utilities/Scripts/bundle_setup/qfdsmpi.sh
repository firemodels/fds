#!/bin/bash -f
EXPECTED_ARGS=3

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: qfdsmpi.sh [-q queue] nthreads fds_command casename.fds"
  echo ""
  echo "Runs a parallel FDS case on a Linux cluster using the "
  echo "PBS/SGE qsub batch queuing command"
  echo ""
  echo " -q queue     - optional parameter used to specify the name of the queue"
  echo                  choices: batch (default), vis,fire60s, fire70s"
  echo "    nthreads - number of threads (usually number of &mesh lines)"
  echo " fds_command - full path to fds command name"
  echo "casename.fds - FDS input file"
  echo 
  exit
fi
queue=batch
ARGS=("$@")
for ((argi=0,i=0;i<$ARGC;i++)) ; do
arg=${ARGS[$i]}
case $arg in
 -q)
   let i=$i+1
   queue=${ARGS[$i]}
  ;;
  *)
   if [ $argi -eq 0 ] ; then
     nthreads=$arg
     let argi=$argi+1
   elif [ $argi -eq 1 ] ; then
     fds=$arg
     let argi=$argi+1
   elif [ $argi -eq 1 ] ; then
     in=$arg
     let argi=$argi+1
   fi
  ;;
esac
done


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
  echo "The fds input file, $fulldir/$in, does not exit. Run aborted."
  exit
fi
if ! [ -e $fds ]; then
  echo "The FDS program name, $fds, does not exit. Run aborted."
  exit
fi
if [ -e $outlog ]; then
  echo "Removing log file: $outlog"
  rm $outlog
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

mpirun -np $nthreads $fds $in
EOF
chmod +x $scriptfile
echo Running $in 
qsub -q $queue $scriptfile
rm $scriptfile
