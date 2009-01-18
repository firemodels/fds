#!/bin/csh -f

set SVNROOT=~/FDS-SMV
set OSXHOST=tiger.cfr.nist.gov
set LINUX64HOST=fire79

echo
echo building 32 bit Linux Smokeview
#cd $SVNROOT/SMV_5/INTEL/sv5p0
cd $SVNROOT/SMV_5/Build/INTEL_LINUX_32
make clean >& /dev/null
date >& $SVNROOT/SMV_5/scripts/make_linux.out
make >>& $SVNROOT/SMV_5/scripts/make_linux.out
cd $SVNROOT/SMV_5/bin

echo 
echo building 64 bit Linux Smokeview
ssh $LINUX64HOST FDS-SMV/SMV_5/scripts/make_smv_64linux.csh
date > $SVNROOT/SMV_5/scripts/make_linux64.out
scp $OSXHOST\:FDS-SMV/SMV_5/bin/make_linux64.out $SVNROOT/SMV_5/scripts/.

echo 
echo building OSX Smokeview
ssh $OSXHOST FDS-SMV/SMV_5/scripts/make_smv_osx.csh
date > $SVNROOT/SMV_5/scripts/make_osx.out
scp $OSXHOST\:FDS-SMV/SMV_5/bin/make_osx.out $SVNROOT/SMV_5/scripts/.
