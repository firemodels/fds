#!/bin/csh -f
set SVNROOT=~/FDS-SMV
#set MACHOST=tiger.cfr.nist.gov
set MACHOST=bluesky.cfr.nist.gov

cd $SVNROOT
echo updating LINUX  repository
svn update

echo updating FIRE72 repository
ssh fire72 \(cd FIRE72/FDS-SMV \; svn update \)

echo updating MAC ( $MACHOST )  repository
ssh $MACHOST \(cd FDS-SMV \; /usr/local/bin/svn update \)

