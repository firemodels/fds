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

echo copying $fds from $fdsdir on $fdshost
scp -q $fdshost\:$fdsroot/$fdsdir/$fds $bundledir/bin/$fdsout

echo copying $fdsmpi from $fdsdir on $fdshost
scp -q $fdshost\:$fdsroot/$fdsmpidir/$fdsmpi $bundledir/bin/$fdsmpiout

echo Building archive
rm -rf $googledir/$bundlebase.tar
rm -rf $googledir/$bundlebase.tar.gz
cd $googledir/$bundlebase
tar cf ../$bundlebase.tar .
echo Compressing archive
gzip    ../$bundlebase.tar
