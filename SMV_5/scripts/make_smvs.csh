#!/bin/csh -f

set SVNROOT=~/FDS-SMV
set OSXHOST=tiger.cfr.nist.gov

echo
echo building Linux smokeview
cd $SVNROOT
cd $SVNROOT/SMV_5/INTEL/sv5p0
make clean >& /dev/null
date >& $SVNROOT/SMV_5/scripts/make_linux.out
make >>& $SVNROOT/SMV_5/scripts/make_linux.out
cd $SVNROOT/SMV_5/bin

echo 
echo building OSX smokeview
ssh $OSXHOST FDS-SMV/SMV_5/scripts/make_smv_osx.csh
date > $SVNROOT/SMV_5/scripts/make_osx.out
scp $OSXHOST\:FDS-SMV/SMV_5/bin/make_osx.out $SVNROOT/SMV_5/scripts/.
