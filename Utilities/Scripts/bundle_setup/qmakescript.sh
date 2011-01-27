#!/bin/bash -f
EXPECTED_ARGS=4

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: qmakescript.sh nthreads nprocs fds_command casename.fds"
  echo ""
  echo "Creates a script used by a batch queuing system to run"
  echo "FDS on a Linux cluster"
  echo ""
  echo "    nthreads - number of threads (usually number of &mesh lines)"
  echo "      nprocs - number of processes/threads per node"
  echo " fds_command - full path to fds command name"
  echo "casename.fds - fds case"
  echo 
  exit
fi
nthreads=$1
nprocs=$2
fds=$3
in=$4
if test $nprocs -gt 8
then
nprocs=8
fi
nnodes=$(echo "($nthreads-1)/$nprocs+1" | bc)
if test $nnodes -le 0
then
nnodes=1
fi

infile=${in%.*}
fulldir=`pwd`

out=$fulldir/$infile.err
outlog=$fulldir/$infile.log

cat << EOF 
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
