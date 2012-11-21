#!/bin/csh -f
set scratchdir=$SVNROOT/Utilities/Scripts/tmp
set dir=$1
set infile=$2

set fulldir=$BASEDIR/$dir
set in=$infile.fds
set out=$infile.err
set stopfile=$infile.stop

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
setenv BACKG ""
if($?BACKGROUND) then
 unsetenv BACKG
 setenv BACKG "$BACKGROUND"
endif
if($?STOPFDS) then
 echo "stopping case: $infile"
 touch $fulldir/$stopfile
 exit
endif
if(-e $fulldir/$stopfile) then
 rm $fulldir/$stopfile
endif
cat << EOF > $scriptfile
#!/bin/csh -f
cd $fulldir
$BACKG $FDS $in >& $out
EOF
chmod +x $scriptfile
echo submitted $in 
$scriptfile ; rm $scriptfile 
