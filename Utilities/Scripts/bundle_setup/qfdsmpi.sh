#!/bin/bash -f
EXPECTED_ARGS=4

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: qfdsmpi.sh nthreads nnodes fds_command casename.fds"
  echo "    nthreads - number of threads"
  echo "      nnodes - number of nodes"
  echo " fds_command - full path to fds command name"
  echo "casename.fds - fds case"
  echo 
  exit
fi
nthreads=$1
nnodes=$2
fds=$3
in=$4
nprocnodes=$(echo "($nthreads-1)/$nnodes+1" | bc)
if test $nprocnodes -le 0
then
nprocnodes=1
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
#PBS -l nodes=$nnodes:ppn=$nprocnodes
cd $fulldir
mpirun -np $nthreads $fds $in
EOF
chmod +x $scriptfile
echo Running $in 
qsub $scriptfile
rm $scriptfile
