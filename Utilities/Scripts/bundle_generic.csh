#!/bin/csh -f

# this script is called by bundle_platform_size.csh
# where platform may be linux or osx and size may be 32 or 64

if $?fdshost then
set scp_fds_smvroot=$fds_smvroot
else
set scp_fds_smvroot=~/$fds_smvroot
endif
set fds_smvroot=~/$fds_smvroot
set makedir=$scp_fds_smvroot/FDS_Compilation
set googledir=$fds_smvroot/Utilities/to_google
set bundledir=$bundlebase/FDS/FDS5
set bundle_setup=$fds_smvroot/Utilities/Scripts/bundle_setup
set mandir=$fds_smvroot/Manuals/All_PDF_Files
set smvbindir=$scp_fds_smvroot/SMV_5/bin
set forbundle=$fds_smvroot/SMV_5/for_bundle
set fds2asciidir=$fds_smvroot/Utilities/fds2ascii

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
cp $forbundle/smokeview.ini $bundledir/bin/.
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
