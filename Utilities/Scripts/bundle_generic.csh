#!/bin/csh -f

# this script is called by bundle_platform_size.csh
# where platform may be linux or osx and size may be 32 or 64

if $?fdshost then
set scp_fds_smvroot=$fds_smvroot
else
set scp_fds_smvroot=~/$fds_smvroot
endif
set fds_smvroot=~/$fds_smvroot
set fdsroot=$scp_fds_smvroot/FDS_Compilation
set smokediffroot=$scp_fds_smvroot/Utilities/smokediff
set smokeziproot=$scp_fds_smvroot/Utilities/smokezip
set googledir=$fds_smvroot/Utilities/to_google
set bundledir=$bundlebase/FDS/FDS5
set bundle_setup=$fds_smvroot/Utilities/Scripts/bundle_setup
set mandir=$fds_smvroot/Manuals/All_PDF_Files
set smvbindir=$scp_fds_smvroot/SMV_5/bin
set texturedir=$scp_fds_smvroot/SMV_5/bin/textures
set forbundle=$fds_smvroot/SMV_5/for_bundle
set fds2asciidir=$fds_smvroot/Utilities/fds2ascii
set wikify=$fds_smvroot/Utilities/Scripts/wikify.py

cd $googledir
rm -rf $bundlebase
mkdir $bundlebase
mkdir $bundlebase/FDS
mkdir $bundledir
mkdir $bundledir/bin
mkdir $bundledir/Documentation
mkdir $bundledir/Examples
mkdir $bundledir/bin/textures

echo Copying program files

# share libraries for INTEL build
if $?INTELLIB then
cp $bundle_setup/README_LINUX.html $bundledir/bin/.
cp -r $bundle_setup/$INTELLIB $bundledir/bin/.
endif

# smokeview

if $?smvhost then
scp $smvhost\:$smvbindir/$smokeview $bundledir/bin/.
else
cp $smvbindir/$smokeview $bundledir/bin/.
endif
cp $texturedir/*.png $bundledir/bin/textures/.
cp $texturedir/*.jpg $bundledir/bin/textures/.

# smokediff

if $?fdshost then
scp $fdshost\:$smokediffroot/$smokediffdir/$smokediff $bundledir/bin/.
else
cp $smokediffroot/$smokediffdir/$smokediff $bundledir/bin/.
endif

# smokezip

if $?fdshost then
scp $fdshost\:$smokeziproot/$smokezipdir/$smokezip $bundledir/bin/.
else
cp $smokeziproot/$smokezipdir/$smokezip $bundledir/bin/.
endif

# FDS 

if $?fdshost then
scp $fdshost\:$fdsroot/$fds5dir/$fds5 $bundledir/bin/$fds5out
scp $fdshost\:$fdsroot/$fds5mpidir/$fds5mpi $bundledir/bin/$fds5mpiout
cp $bundle_setup/FDS-SMV_5_OSX_Launcher.app.zip $bundledir/bin/.
cp $bundle_setup/README_OSX.html $bundledir/bin/.
else
cp $fdsroot/$fds5dir/$fds5 $bundledir/bin/$fds5out
cp $fdsroot/$fds5mpidir/$fds5mpi $bundledir/bin/$fds5mpiout
endif

cp $forbundle/smokeview.ini $bundledir/bin/.
cp $fds2asciidir/$fds2ascii $bundledir/bin/.

echo Copying documentation
cp $bundle_setup/Overview.html $bundledir/Documentation/.
cp $mandir/FDS_5_User_Guide.pdf $bundledir/Documentation/.
cp $mandir/SMV_5_User_Guide.pdf $bundledir/Documentation/.

echo Copy objects.svo
cp $forbundle/objects.svo $bundledir/bin/objects.svo

echo
echo Getting the FDS release notes from the repository
svn export --quiet --force http://fds-smv.googlecode.com/svn/wiki/FDS_Release_Notes.wiki $bundle_setup/FDS_Release_Notes.wiki

echo 
echo Converting the FDS release notes from wiki to html format
$wikify -r $bundle_setup/FDS_Release_Notes.wiki > $bundledir/Documentation/FDS_Release_Notes.html

echo
echo Copying Smokeview release from  the repository
cp $forbundle/readme.html $bundledir/Documentation/SMV_Release_Notes.html


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
