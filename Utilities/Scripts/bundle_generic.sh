#!/bin/bash

# this script is called by bundle_platform_size.csh
# where platform may be linux or osx and size may be 32 or 64

SCP ()
{
  HOST=$1
  FROMDIR=$2
  FROMFILE=$3
  TODIR=$4
  TOFILE=$5

  scp $HOST\:$FROMDIR/$FROMFILE $TODIR/$TOFILE 2>/dev/null
  if [ -e $TODIR/$TOFILE ]; then
    echo "$TOFILE copied from $HOST"
  else
    echo "***error: the file $TOFILE failed to copy from $HOST"
  fi
}

CP ()
{
  FROMDIR=$1
  FROMFILE=$2
  TODIR=$3
  TOFILE=$4
  if [ ! -e $FROMDIR/$FROMFILE ]; then
    echo "***error: the file $FROMFILE does not exist"
  else
    cp $FROMDIR/$FROMFILE $TODIR/$TOFILE
  fi
  if [ -e $TODIR/$TOFILE ]; then
    echo "$TOFILE copied"
  else
    echo "***error: the file $TOFILE failed to copy"
  fi
}

CP2 ()
{
  FROMDIR=$1
  FROMFILE=$2
  TODIR=$3
  TOFILE=$FROMDIR
  if [ ! -e $FROMDIR/$FROMFILE ]; then
    echo "***error: the file $FROMFILE does not exist"
  else
    cp $FROMDIR/$FROMFILE $TODIR/$TOFILE
  fi
  if [ -e $TODIR/$TOFILE ]; then
    echo "$TOFILE copied"
  else
    echo "***error: the file $TOFILE failed to copy"
  fi
}

CPDIR ()
{
  FROMDIR=$1
  TODIR=$2
  if [ ! -e $FROMDIR ]; then
    echo "***error: the directory $FROMDIR does not exist"
  else
    cp -r $FROMDIR $TODIR
  fi
  if [ -e $TODIR ]; then
    echo "$TODIR copied"
  else
    echo "***error: the directory $TODIR failed to copy"
  fi
}


manifest=manifest$FDSOS.html
OUT=$MAJOR$FDSOS
OUT=

smokeviewdir=intel$FDSOS
smokeview=smokeview$FDSOS
smokeviewout=smokeview$OUT

smokezipdir=intel$FDSOS
smokezip=smokezip$FDSOS
smokezipout=smokezip$OUT

wind2fdsdir=intel$FDSOS
wind2fds=wind2fds$FDSOS
wind2fdsout=wind2fds$OUT

smokediffdir=intel$FDSOS
smokediff=smokediff$FDSOS
smokediffout=smokediff$OUT

backgrounddir=intel$FDSOS
background=background
backgroundout=background

fdsdir=intel$FDSOS
fds=fds_intel$FDSOS
fdsout=fds$OUT

fdsmpidir=mpi_intel$FDSOS
fdsmpi=fds_mpi_intel$FDSOS
fdsmpiout=fds$MAJOR\_mpi$FDSOS
fdsmpiout=fds_mpi
fdsmpiout=fds$OUT

fds2asciidir=intel$FDSOS
fds2ascii=fds2ascii$FDSOS
fds2asciiout=fds2ascii$OUT


scp_fds_smvroot=$fds_smvroot
fds_smvroot=~/$fds_smvroot
fdsroot=$scp_fds_smvroot/FDS_Compilation
backgroundroot=$scp_fds_smvroot/SMV/Build/background
smokediffroot=$scp_fds_smvroot/SMV/Build/smokediff
smokeziproot=$scp_fds_smvroot/SMV/Build/smokezip
wind2fdsroot=$scp_fds_smvroot/SMV/Build/wind2fds
uploaddir=$fds_smvroot/Utilities/uploads
bundledir=$bundlebase
bundle_setup=$fds_smvroot/Utilities/Scripts/bundle_setup
mandir=~/FIRE-LOCAL/reports/fds_manuals
smvbindir=$scp_fds_smvroot/SMV/Build/smokeview/$smokeviewdir
forbundle=$fds_smvroot/SMV/for_bundle
texturedir=$forbundle/textures
fds2asciiroot=$scp_fds_smvroot/Utilities/fds2ascii
wikify=$fds_smvroot/Utilities/Scripts/wikify.py
fullmanifest=$uploaddir/$bundledir/bin/$manifest
makeinstaller=$fds_smvroot/Utilities/Scripts/make_installer.sh

fds_cases=$fds_smvroot/Verification/FDS_Cases.sh
smv_cases=$fds_smvroot/Verification/scripts/SMV_Cases.sh
wui_cases=$fds_smvroot/Verification/scripts/WUI_Cases.sh
copyfdscase=$fds_smvroot/Utilities/Scripts/copyfdscase.sh
copycfastcase=$fds_smvroot/Utilities/Scripts/copycfastcase.sh
ExamplesDirectory=$fds_smvroot/Verification

cd $uploaddir
rm -rf $bundlebase
mkdir $bundledir
mkdir $bundledir/bin
mkdir $bundledir/Documentation
mkdir $bundledir/Examples
mkdir $bundledir/bin/textures

# background

SCP $fdshost $backgroundroot/$backgrounddir $background $bundledir/bin $backgroundout

# smokeview

SCP $smvhost $smvbindir $smokeview $bundledir/bin $smokeviewout

# textures

CPDIR $texturedir $bundledir/bin/textures

# smokediff

SCP $fdshost $smokediffroot/$smokediffdir $smokediff $bundledir/bin $smokediffout

# smokezip

SCP $fdshost $smokeziproot/$smokezipdir $smokezip $bundledir/bin $smokezipout

# wind2fds

SCP $fdshost $wind2fdsroot/$wind2fdsdir $wind2fds $bundledir/bin $wind2fdsout

# FDS 

SCP $fdshost $fdsroot/$fdsmpidir $fdsmpi $bundledir/bin $fdsmpiout

if [ "$PLATFORM" == "LINUX64" ]; then
   ostype=LINUX
   ossize=intel64
fi
if [ "$PLATFORM" == "OSX64" ]; then
   ostype=OSX
   ossize=intel64
fi

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
ssh -q $runhost " echo 0 | $fdsroot/$fdsmpidir/$fdsmpi" >> $fullmanifest

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

if [ "$OSXBUNDLE" == "yes" ]; then
  CP $bundle_setup FDS-SMV_OSX_Launcher.app.zip $bundledir/bin FDS-SMV_OSX_Launcher.app.zip
  CP $bundle_setup README_OSX.html $bundledir/bin README_OSX.html
fi

CP $forbundle smokeview.ini $bundledir/bin smokeview.ini

CP $forbundle volrender.ssf $bundledir/bin volrender.ssf

CP $forbundle objects.svo $bundledir/bin objects.svo

SCP $fdshost $fds2asciiroot/$fds2asciidir $fds2ascii $bundledir/bin $fds2asciiout

echo Copying documentation
CP $bundle_setup Overview_linux_osx.html $bundledir/Documentation Overview.html
CP2 $mandir FDS_Configuration_Management_Plan.pdf $bundledir/Documentation
CP2 $mandir FDS_Technical_Reference_Guide.pdf $bundledir/Documentation
CP2 $mandir FDS_User_Guide.pdf $bundledir/Documentation
CP2 $mandir FDS_Validation_Guide.pdf $bundledir/Documentation
CP2 $mandir FDS_Verification_Guide.pdf $bundledir/Documentation
CP2 $mandir SMV_User_Guide.pdf $bundledir/Documentation
CP2 $mandir SMV_Technical_Reference_Guide.pdf $bundledir/Documentation
CP2 $mandir SMV_Verification_Guide.pdf $bundledir/Documentation

if [ ! "$INTELLIB" == "" ]; then
if [ -d $INTELLIB ]; then
echo copying  run time libraries
CPDIR $INTELLIB $bundledir/bin/$DESTLIB
fi
fi

CP $bundle_setup FDS_Release_Notes.htm $bundledir/Documentation FDS_Release_Notes.html

CP ~/FDS-SMVwebpages smv_readme.html $bundledir/Documentation SMV_Release_Notes.html


CP2 $bundle_setup readme_examples.html $bundledir/Examples

cd $ExamplesDirectory
export OUTDIR=$uploaddir/$bundledir/Examples
export QFDS=$copyfdscase
export RUNTFDS=$copyfdscase
export RUNCFAST=$copycfastcase
echo Copying example files to bundle directory
$fds_cases
$smv_cases
$wui_cases
rm -rf $OUTDIR/Immersed_Boundary_Method

echo >> $fullmanifest
echo ------file listing---------------------------------- >> $fullmanifest
set curdir=`pwd`
cd $uploaddir/$bundledir
find . -print >> $fullmanifest
cd $curdir

cat <<EOF>>$fullmanifest
</pre>
</body>
</html>
EOF

cp $fullmanifest ~/$manifest
cp $fullmanifest $uploaddir/$manifest
cat $fullmanifest | Mail -s " $PLATFORM" `whoami`

echo Building archive
rm -rf $uploaddir/$bundlebase.tar
rm -rf $uploaddir/$bundlebase.tar.gz
cd $uploaddir/$bundlebase
tar cf ../$bundlebase.tar --exclude='*.csv' .
echo Compressing archive
gzip    ../$bundlebase.tar
echo Creating installer
cd ..
$makeinstaller -o $ostype -i $bundlebase.tar.gz -d $INSTALLDIR $bundlebase.sh 

