#!/bin/bash -f
EXPECTED_ARGS=1

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: qfds.sh casename.fds"
  echo ""
  echo "Runs a serial 32 bit FDS case on a Linux cluster using the "
  echo "PBS/SGE qsub batch queuing command"
  echo ""
  echo " casename.fds - fds case"
  echo 
  exit
fi
in=$1

infile=${in%.*}
fulldir=`pwd`

out=$fulldir/$infile.err
outlog=$fulldir/$infile.log

scriptfile=/tmp/script.$$
if ! [ -e $fulldir/$in ]; then
  echo "The fds input file, $fulldir/$in, does not exit. Run aborted."
  exit
fi
source ~/.bashrc_fds ia32
$FDSBINDIR/qfds.sh $FDSBINDIR/fds_linux_32 $in
