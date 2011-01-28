#!/bin/bash -f
scratchdir=$SVNROOT/Utilities/Scripts/tmp
nthreads=$1
dir=$2
infile=$3
nprocs=8
nnodes=$(echo "($nthreads-1)/$nprocs+1" | bc)
if test $nnodes -le 0
then
nnodes=1
fi


fulldir=$BASEDIR/$dir
in=$infile.fds
outerr=$fulldir/$infile.err
outlog=$fulldir/$infile.log
stopfile=$infile.stop

scriptfile=$scratchdir/script.$$
if ! [ -e $FDS ];  then
  echo "The file $FDS does not exit. Run aborted"
  exit
fi
if ! [ -d $fulldir ]; then
  echo "The directory $fulldir does not exit. Run aborted."
  exit
fi
if ! [ -e $fulldir/$in ]; then
  echo "The fds input file, $fulldir/$in, does not exit. Run aborted."
  exit
fi
if [ $STOPFDS ]; then
 echo "stopping case: $infile"
 touch $fulldir/$stopfile
 exit
fi
if [ -e $fulldir/$stopfile ]; then
 rm $fulldir/$stopfile
fi
if [ -e $outlog ]; then
 rm $outlog
fi
cat << EOF > $scriptfile
#!/bin/bash -f
#PBS -N VV_$infile(MPI)
#PBS -l nodes=$nnodes:ppn=$nprocs
#PBS -S /bin/bash
#PBS -e $outerr
#PBS -o $outlog
#\$ -N VV_$infile(MPI)
#\$ -l nodes=$nnodes:ppn=$nprocs
#\$ -S /bin/bash
#\$ -e $outerr
#\$ -o $outlog

cd $fulldir

echo Time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`

mpirun -np $nthreads $FDS $in 
EOF
chmod +x $scriptfile
echo Running $in 
qsub $scriptfile
rm $scriptfile
