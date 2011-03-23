#!/bin/bash -f
EXPECTED_ARGS=2

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: qfds.sh fds_command casename.fds"
  echo ""
  echo "Runs a serial FDS case on a Linux cluster using the "
  echo "PBS/SGE qsub batch queuing command (SVN $Revision$')"
  echo ""
  echo "  fds_command - full path to fds command name"
  echo " casename.fds - fds case"
  echo 
  exit
fi
fds=$1
in=$2

infile=${in%.*}
fulldir=`pwd`

err=$fulldir/$infile.err
err2=$fulldir/$infile.err2
outlog=$fulldir/$infile.log

stopfile=${in%.*}.stop
if [ $STOPFDS ]; then
 echo "stopping case: $infile"
 touch $stopfile
 exit
fi
if [ -e $stopfile ]; then
  rm $stopfile
fi
scriptfile=/tmp/script.$$
if ! [ -e $fulldir/$in ]; then
  echo "The fds input file, $fulldir/$in, does not exist. Run aborted."
  exit
fi
if ! [ -e $fds ]; then
  echo "The FDS program name, $fds, does not exist. Run aborted."
  exit
fi
if [ -e $outlog ]; then
  echo "Removing log file: $outlog"
  rm $outlog
fi

cat << EOF > $scriptfile
#!/bin/bash -f
#PBS -N $infile
#PBS -e $err2
#PBS -o $outlog
#\$ -N $infile
#\$ -e $err2
#\$ -o $outlog

cd $fulldir

echo Start time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`

$fds $in 2>$err
EOF
chmod +x $scriptfile
echo Running $in 
qsub $scriptfile
rm $scriptfile
