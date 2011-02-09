#!/bin/bash -f
EXPECTED_ARGS=4

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: qfdsmpi.sh nthreads platform fds_command casename.fds"
  echo ""
  echo "Runs a parallel FDS case on a Linux cluster using the "
  echo "PBS/SGE qsub batch queuing command"
  echo ""
  echo "    nthreads - number of threads (usually number of &mesh lines)"
  echo "    platform - ia32 or intel64"
  echo " fds_command - full path to fds command name"
  echo "casename.fds - fds case"
  echo 
  exit
fi
nthreads=$1
platform=$2
fds=$3
in=$4
if test $nthreads -le 0
then
echo "Number of threads specified is $nthreads . Must be bigger than 0."
exit
fi
nnodes=$(echo "($nthreads-1)/8+1" | bc)
nprocs=$(echo "($nthreads-1)/$nnodes+1" | bc)

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

source ~/.bashrc_fds $platform
mpirun -np $nthreads $fds $in
EOF
chmod +x $scriptfile
echo Running $in 
qsub $scriptfile
rm $scriptfile
