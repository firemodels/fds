#!/bin/csh -f
set scratchdir=$SVNROOT/Utilities/Scripts/tmp
set dir=$1
set infile=$2

setenv fulldir $BASEDIR/$dir
set in=$infile.fds
set out=$infile.err

set scriptfile=$scratchdir/script.$$

if(! -e $FDS) then
  echo "The file $FDS does not exist. Run aborted"
  exit
endif
if(! -d $fulldir) then
  echo "The directory $fulldir does not exist. Run aborted."
  exit
endif
if(! -e $fulldir/$in) then
  echo "The fds input  files, $fulldir/$in, does not exist. Run aborted."
  exit
endif

#cd $fulldir
#qsub sge-fds.sh

cat << EOF > $scriptfile
#!/bin/bash
#\$ -S /bin/bash
#\$ -cwd -N $infile -V -e /dev/null -o /dev/null
cd $fulldir
$FDS $in >& $out
EOF
chmod +x $scriptfile
qsub $scriptfile
rm $scriptfile

