#!/bin/csh -f

# this script is called by bundle_platform_size.csh
# where platform may be linux or osx and size may be 32 or 64

setenv manifest manifest$FDSOS.html

setenv smokeviewdir intel$FDSOS
setenv smokeview smokeview$FDSOS
setenv smokeviewout smokeview$MAJOR$FDSOS
echo smokeview=$smokeview
echo smokeviewout=$smokeviewout

setenv smokezipdir intel$FDSOS
setenv smokezip smokezip$FDSOS
setenv smokezipout smokezip$MAJOR$FDSOS

setenv smokediffdir intel$FDSOS
setenv smokediff smokediff$FDSOS
setenv smokediffout smokediff$MAJOR$FDSOS

setenv fdsdir intel$FDSOS
setenv fds fds_intel$FDSOS
setenv fdsout fds$MAJOR$FDSOS

setenv fdsmpidir mpi_intel$FDSOS
setenv fdsmpi fds_mpi_intel$FDSOS
setenv fdsmpiout fds$MAJOR\_mpi$FDSOS

setenv fds2asciidir intel$FDSOS
setenv fds2ascii fds2ascii$FDSOS
setenv fds2asciiout fds2ascii$MAJOR$FDSOS


set scp_fds_smvroot=$fds_smvroot
set fds_smvroot=~/$fds_smvroot
set fdsroot=$scp_fds_smvroot/FDS_Compilation
set smokediffroot=$scp_fds_smvroot/Utilities/smokediff
set smokeziproot=$scp_fds_smvroot/Utilities/smokezip
set googledir=$fds_smvroot/Utilities/to_google
set bundledir=$bundlebase
set bundle_setup=$fds_smvroot/Utilities/Scripts/bundle_setup
set mandir=$fds_smvroot/Manuals/All_PDF_Files
set smvbindir=$scp_fds_smvroot/SMV/Build/$smokeviewdir
set forbundle=$fds_smvroot/SMV/for_bundle
set texturedir=$forbundle/textures
set fds2asciiroot=$scp_fds_smvroot/Utilities/fds2ascii
set wikify=$fds_smvroot/Utilities/Scripts/wikify.py
set fullmanifest=$googledir/$bundledir/bin/$manifest
set makeinstaller=$fds_smvroot/Utilities/Scripts/make_installer.sh

cd $googledir
rm -rf $bundlebase
mkdir $bundledir
mkdir $bundledir/bin
mkdir $bundledir/Documentation
mkdir $bundledir/Examples
mkdir $bundledir/bin/textures

echo Copying program files

# share libraries for INTEL build
if $?INTELLIB then
cp $bundle_setup/README_LINUX.html $bundledir/bin/.
cp -r $INTELLIB $bundledir/bin/.
set SETLDPATH="LD_LIBRARY_PATH=$INTELLIB;" 
else
set SETLDPATH=
endif

# smokeview

echo copying $smokeview from $smvbindir on $smvhost
scp -q $smvhost\:$smvbindir/$smokeview $bundledir/bin/$smokeviewout

echo copying textures
cp $texturedir/*.png $bundledir/bin/textures/.
cp $texturedir/*.jpg $bundledir/bin/textures/.

# smokediff

echo copying $smokediff from $smokediffdir on $fdshost
scp -q $fdshost\:$smokediffroot/$smokediffdir/$smokediff $bundledir/bin/$smokediffout

# smokezip

echo copying $smokezip from $smokezipdir on $fdshost
scp -q $fdshost\:$smokeziproot/$smokezipdir/$smokezip $bundledir/bin/$smokezipout

# FDS 

echo copying $fds from $fdsdir on $fdshost
scp -q $fdshost\:$fdsroot/$fdsdir/$fds $bundledir/bin/$fdsout

echo copying $fdsmpi from $fdsdir on $fdshost
scp -q $fdshost\:$fdsroot/$fdsmpidir/$fdsmpi $bundledir/bin/$fdsmpiout

if ($PLATFORM == "LINUX32" || $PLATFORM == "LINUX64") then
set ostype=LINUX
endif
if ($PLATFORM == "OSX32" || $PLATFORM == "OSX64") then
set ostype=OSX
endif
if ($PLATFORM == "LINUX32" || $PLATFORM == "OSX32") then
set ossize=ia32
endif
if ($PLATFORM == "LINUX64" || $PLATFORM == "OSX64") then
set ossize=intel64
endif

cat <<EOF > $fullmanifest
<html>
<head>
<TITLE>FDS-SMV Bundle Manifest</TITLE>
</HEAD>
<BODY BGCOLOR="#FFFFFF" >
<pre>
EOF
echo
echo Creating Manifest
echo $PLATFORM FDS-Smokeview bundle created >> $fullmanifest
date >> $fullmanifest
echo  >> $fullmanifest
echo Versions:>> $fullmanifest
echo  >> $fullmanifest
echo ------fds-------------------- >> $fullmanifest
ssh -q $runhost " echo 0 | $fdsroot/$fdsdir/$fds" >>& $fullmanifest

echo  >> $fullmanifest
echo ------fds2ascii-------------------- >> $fullmanifest
ssh -q $runhost $fds2asciiroot/$fds2asciidir/$fds2ascii -v >> $fullmanifest

echo  >> $fullmanifest
echo ------smokeview-------------------- >> $fullmanifest
ssh -q $runhost $smvbindir/$smokeview -v  >> $fullmanifest
echo  >> $fullmanifest
echo ------smokediff-------------------- >> $fullmanifest
ssh -q $runhost $smokediffroot/$smokediffdir/$smokediff -v >> $fullmanifest
echo  >> $fullmanifest
echo ------smokezip-------------------- >> $fullmanifest
ssh -q $runhost $smokeziproot/$smokezipdir/$smokezip -v >> $fullmanifest

if ($?OSXBUNDLE) then
echo copying OSX launcher script
cp $bundle_setup/FDS-SMV_OSX_Launcher.app.zip $bundledir/bin/.
cp $bundle_setup/README_OSX.html $bundledir/bin/.
endif

echo copying smokeview.ini from $forbundle
cp $forbundle/smokeview.ini $bundledir/bin/.

echo copying $fds2ascii from $fds2asciiroot on $fdshost
scp -q $fdshost\:$fds2asciiroot/$fds2asciidir/$fds2ascii $bundledir/bin/$fds2asciiout

echo Copying documentation
cp $bundle_setup/Overview.html $bundledir/Documentation/.
cp $mandir/FDS_User_Guide.pdf $bundledir/Documentation/.
cp $mandir/SMV_User_Guide.pdf $bundledir/Documentation/.
cp $mandir/SMV_Technical_Reference_Guide.pdf $bundledir/Documentation/.
cp $mandir/FDS_Technical_Reference_Guide.pdf $bundledir/Documentation/.

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

echo >> $fullmanifest
echo ------file listing---------------------------------- >> $fullmanifest
set curdir=`pwd`
cd $googledir/$bundledir
find . -print >> $fullmanifest
cd $curdir

cat <<EOF>>$fullmanifest
</pre>
</body>
</html>
EOF

cp $fullmanifest ~/$manifest
cp $fullmanifest $googledir/$manifest
cat $fullmanifest | Mail -s " $PLATFORM" `whoami`

echo Building archive
rm -rf $googledir/$bundlebase.tar
rm -rf $googledir/$bundlebase.tar.gz
cd $googledir/$bundlebase
tar cf ../$bundlebase.tar --exclude='*.csv' .
echo Compressing archive
gzip    ../$bundlebase.tar
echo Creating installer
cd ..
$makeinstaller $ostype $ossize $bundlebase.tar.gz $bundlebase.sh \~/FDS/$FDSEDITION
