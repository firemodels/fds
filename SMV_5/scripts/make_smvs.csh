#!/bin/csh -f

set SVNROOT=~/FDS-SMV
set OSXHOST=tiger.cfr.nist.gov

echo
echo building Linux smokeview
cd $SVNROOT
svn update >& svn.out
cd $SVNROOT/SMV_5/INTEL/sv5p0
make clean >& /dev/null
make >& $SVNROOT/SMV_5/scripts/make_linux.out
cd $SVNROOT/SMV_5/bin
ls -l smv5_linux

echo 
echo building OSX smokeview
ssh $OSXHOST FDS-SMV/SMV_5/scripts/make_smv_osx.csh
scp $OSXHOST\:FDS-SMV/SMV_5/bin/make_osx.out $SVNROOT/SMV_5/scripts/.
