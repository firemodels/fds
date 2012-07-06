#!/bin/bash -f
scratchdir=$SVNROOT/Utilities/Scripts/tmp
dir=$1
infile=$2

fulldir=$BASEDIR/$dir
in=$infile
outerr=$fulldir/$infile.err
outlog=$fulldir/$infile.log

scriptfile=$scratchdir/script.$$
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
cat << EOF > $scriptfile
#!/bin/bash -f
#\$ -S /bin/bash
#\$ -N VV_$infile -e $outerr -o $outlog
#PBS -N VV_$infile -e $outerr -o $outlog
cd $fulldir

echo Time: \`date\`
echo Running $infile on \`hostname\`
echo Directory: \`pwd\`

$CFAST $in 
EOF
chmod +x $scriptfile
echo Running $in 
qsub -q fire70s $scriptfile
rm $scriptfile
