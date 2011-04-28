#!/bin/bash -f
EXPECTED_ARGS=2

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: qfdsmpi.sh nthreads casename.fds"
  echo ""
  echo "Runs a 64 bit parallel FDS case on a Linux cluster using the "
  echo "PBS/SGE qsub batch queuing command"
  echo ""
  echo "    nthreads - number of threads (usually number of &mesh lines)"
  echo "casename.fds - fds case"
  echo 
  exit
fi
nthreads=$1
in=$2
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

source ~/.bashrc_fds intel64
$FDSBINDIR/qfdsmpi.sh $nthreads intel64 $FDSBINDIR/fds_mpi_linux_64 $in
