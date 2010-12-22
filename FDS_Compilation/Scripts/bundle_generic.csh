#!/bin/csh -f

# this script is called by bundle_platform_size.csh
# where platform may be linux or osx and size may be 32 or 64

set scp_fds_smvroot=$fds_smvroot
set fds_smvroot=~/$fds_smvroot
set fdsroot=$scp_fds_smvroot/FDS_Compilation
set googledir=$fds_smvroot/Utilities/to_google
set bundledir=$bundlebase/FDS/FDS5

cd $googledir
rm -rf $bundlebase
mkdir $bundlebase
mkdir $bundlebase/FDS
mkdir $bundledir
mkdir $bundledir/bin

# FDS 

echo copying $fds5 from $fds5dir on $fdshost
scp -q $fdshost\:$fdsroot/$fds5dir/$fds5 $bundledir/bin/$fds5out

echo copying $fds5mpi from $fds5dir on $fdshost
scp -q $fdshost\:$fdsroot/$fds5mpidir/$fds5mpi $bundledir/bin/$fds5mpiout

echo Building archive
rm -rf $googledir/$bundlebase.tar
rm -rf $googledir/$bundlebase.tar.gz
cd $googledir/$bundlebase
tar cf ../$bundlebase.tar .
echo Compressing archive
gzip    ../$bundlebase.tar
