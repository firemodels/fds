#!/bin/csh -f

# this script is called by bundle_platform_size.csh
# where platform may be linux or osx and size may be 32 or 64

cd $googledir
rm -rf $bundlebase
mkdir $bundlebase
mkdir $bundlebase/FDS
mkdir $bundledir
mkdir $bundledir/bin
mkdir $bundledir/Documentation
mkdir $bundledir/Examples

echo Copying program files
if $?fdshost then
scp $fdshost\:$makedir/$fds5dir/$fds5 $bundledir/bin/.
scp $fdshost\:$makedir/$fds5mpidir/$fds5mpi $bundledir/bin/.
scp $fdshost\:$smvbindir/$smokeview $bundledir/bin/.
scp $fdshost\:$smvbindir/$smokezip $bundledir/bin/.
else
cp $makedir/$fds5dir/$fds5 $bundledir/bin/.
cp $makedir/$fds5mpidir/$fds5mpi $bundledir/bin/.
cp $smvbindir/$smokeview $bundledir/bin/.
cp $smvbindir/$smokezip $bundledir/bin/.
endif
cp $smvbindir/smokeview.ini $bundledir/bin/.
cp $fds2asciidir/$fds2ascii $bundledir/bin/.

echo Copying documentation
cp $bundle_setup/readme_docs.html $bundledir/Documentation/.
cp $mandir/FDS_5_User_Guide.pdf $bundledir/Documentation/.
cp $mandir/SMV_5_User_Guide.pdf $bundledir/Documentation/.

echo Obtaining example files from the repository
cp $bundle_setup/readme_examples.html $bundledir/Examples/.
svn export -q --force https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Verification $bundledir/Examples/.

echo Building archive
rm -rf $googledir/$bundlebase.tar
rm -rf $googledir/$bundlebase.tar.gz
cd $googledir/$bundlebase
tar cf ../$bundlebase.tar --exclude='*.csv' .
echo Compressing archive
gzip    ../$bundlebase.tar
