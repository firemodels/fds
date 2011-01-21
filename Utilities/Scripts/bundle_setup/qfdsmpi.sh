#!/bin/bash -f
EXPECTED_ARGS=3

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: qfdsmpi.sh nthreads fds_command casename.fds"
  echo "    nthreads - number of threads (usually number of &mesh lines)"
  echo " fds_command - full path to fds command name"
  echo "casename.fds - fds case"
  echo 
  exit
fi
nthreads=$1
fds=$2
in=$3
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
echo Processors: \`cat \$PBS_NODEFILE\`

mpirun -np $nthreads $fds $in
EOF
chmod +x $scriptfile
echo Running $in 
qsub $scriptfile
#rm $scriptfile
