#!/bin/csh -f

set SVNROOT=~/FDS-SMV
set OSXHOST=tiger.cfr.nist.gov
set LINUX64HOST=fire79

echo
echo building 32 bit Linux Smokeview
#cd $SVNROOT/SMV_5/INTEL/sv5p0
cd $SVNROOT/SMV_5/Build/INTEL_LINUX_32
make -f ../Makefile clean >& /dev/null
date >& $SVNROOT/SMV_5/scripts/make_intel_linux_32.out
./make_smv.csh >>& $SVNROOT/SMV_5/scripts/make_intel_linux_32.out
cd $SVNROOT/SMV_5/bin

echo 
echo building 64 bit Linux Smokeview
ssh $LINUX64HOST FDS-SMV/SMV_5/scripts/make_smv_intel_linux_64.csh
scp $LINUX64HOST\:FDS-SMV/SMV_5/bin/make_intel_linux_64.out $SVNROOT/SMV_5/scripts/.

echo 
echo building 32 bit OSX Smokeview
ssh $OSXHOST FDS-SMV/SMV_5/scripts/make_smv_osx.csh
scp $OSXHOST\:FDS-SMV/SMV_5/bin/make_intel_osx_32.out $SVNROOT/SMV_5/scripts/.
