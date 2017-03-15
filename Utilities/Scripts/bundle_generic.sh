#!/bin/bash

# this script is called by BUNDLE_linux64.sh and BUNDLE_osx64.sh

errlog=/tmp/errlog.$$
OPENMPI_VERSION=2.0.2

# ---------------- MDHASH -----------------

MD5HASH ()
{
local DIR=$1
local FILE=$2

local curdir=`pwd`
cd $DIR

if [ "`uname`" == "Darwin" ] ; then
  hash=`cat $FILE | md5`
  echo "$hash   $FILE"
else
  md5sum $FILE
fi

cd $curdir
}


# -------------------- SCP -------------------

SCP ()
{
  HOST=$1
  FROMDIR=$2
  FROMFILE=$3
  TODIR=$4
  TOFILE=$5

  scp $HOST\:$FROMDIR/$FROMFILE $TODIR/$TOFILE 2>/dev/null
  if [ -e $TODIR/$TOFILE ]; then
    echo "$FROMFILE copied from host:$HOST"
  else
    echo "***error: $TOFILE not copied to bundle from $HOST at $FROMDIR/$FROMFILE " >> $errlog
  fi
}

# -------------------- CP -------------------

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
    echo "$FROMFILE copied"
  else
    echo "***error: $FROMFILE not copied to bundle" >> $errlog
  fi
}

# -------------------- UNTAR -------------------

UNTAR ()
{
  FROMDIR=$1
  FROMFILE=$2
  TODIR=$3
  TODIR2=$4
  if [ ! -e $FROMDIR/$FROMFILE ]; then
    echo "***error: the compressed file $FROMFILE does not exist"
  else
    curdir=`pwd`
    cd $TODIR
    tar xvf $FROMDIR/$FROMFILE
    cd $curdir
  fi
  if [ -e $TODIR/$TODIR2 ]; then
    echo "$FROMFILE untar'd"
  else
    echo "***error: $FROMFILE not untar'd to bundle" >> $errlog
  fi
}

# -------------------- CP2 -------------------

CP2 ()
{
  FROMDIR=$1
  FROMFILE=$2
  TODIR=$3
  TOFILE=$FROMFILE
  if [ ! -e $FROMDIR/$FROMFILE ]; then
    echo "***error: the file $FROMFILE does not exist"
  else
    cp $FROMDIR/$FROMFILE $TODIR/$TOFILE
  fi
  if [ -e $TODIR/$TOFILE ]; then
    echo "$FROMFILE copied"
  else
    echo "***error: $FROMFILE not copied to bundle" >> $errlog
  fi
}

# -------------------- CPDIR -------------------

CPDIR ()
{
  FROMDIR=$1
  TODIR=$2
  if [ ! -e $FROMDIR ]; then
    echo "***error: the directory $FROMDIR does not exist"
  else
    echo "*******************************"
    echo copying directory from $FROMDIR to $TODIR
    echo "*******************************"
    cp -r $FROMDIR $TODIR
  fi
  if [ -e $TODIR ]; then
    echo "$FROMDIR copied"
  else
    echo "***error: the directory $FROMDIR not copied to bundle" >> $errlog
  fi
}

# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

manifest=manifest$FDSOS.html
OUT=$MAJOR$FDSOS
OUT=

smokeviewdir=intel$FDSOS
smokeview=smokeview$FDSOS
smokeviewout=smokeview$OUT

smokezipdir=intel$FDSOS
smokezip=smokezip$FDSOS
smokezipout=smokezip$OUT

dem2fdsdir=intel$FDSOS
dem2fds=dem2fds$FDSOS
dem2fdsout=dem2fds$OUT

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

openmpidir=~/FDS_Guides

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
fdsroot=$scp_fds_smvroot/fds/Build
backgroundroot=$scp_fds_smvroot/smv/Build/background
smokediffroot=$scp_fds_smvroot/smv/Build/smokediff
smokeziproot=$scp_fds_smvroot/smv/Build/smokezip
dem2fdsroot=$scp_fds_smvroot/smv/Build/dem2fds
smvscriptdir=$scp_fds_smvroot/smv/scripts
wind2fdsroot=$scp_fds_smvroot/smv/Build/wind2fds
uploaddir=$fds_smvroot/fds/Utilities/uploads
bundledir=$bundlebase
webpagesdir=$fds_smvroot/webpages
mandir=~/FDS_Guides
smvbindir=$scp_fds_smvroot/smv/Build/smokeview/$smokeviewdir
fds_bundle=$fds_smvroot/fds/Utilities/Scripts/for_bundle
smv_bundle=$fds_smvroot/smv/for_bundle
texturedir=$smv_bundle/textures
fds2asciiroot=$scp_fds_smvroot/fds/Utilities/fds2ascii
makeinstaller=$fds_smvroot/fds/Utilities/Scripts/make_installer.sh

fds_cases=$fds_smvroot/fds/Verification/FDS_Cases.sh
fds_benchamrk_cases=$fds_smvroot/fds/Verification/FDS_Benchmark_Cases.sh
smv_cases=$fds_smvroot/smv/Verification/scripts/SMV_Cases.sh
wui_cases=$fds_smvroot/smv/Verification/scripts/WUI_Cases.sh
copyfdscase=$fds_smvroot/fds/Utilities/Scripts/copyfdscase.sh
copycfastcase=$fds_smvroot/fds/Utilities/Scripts/copycfastcase.sh
FDSExamplesDirectory=$fds_smvroot/fds/Verification
SMVExamplesDirectory=$fds_smvroot/smv/Verification

cd $uploaddir
rm -rf $bundlebase
mkdir $bundledir
mkdir $bundledir/bin
mkdir $bundledir/bin/MD5
mkdir $bundledir/Documentation
mkdir $bundledir/Examples
mkdir $bundledir/bin/textures

echo ""
echo "--- copying programs ---"
echo ""

# smokeview

SCP $fdshost $backgroundroot/$backgrounddir $background $bundledir/bin $backgroundout
SCP $smvhost $smvbindir                     $smokeview  $bundledir/bin $smokeviewout
SCP $fdshost $smokediffroot/$smokediffdir   $smokediff  $bundledir/bin $smokediffout
SCP $fdshost $smokeziproot/$smokezipdir     $smokezip   $bundledir/bin $smokezipout
SCP $fdshost $dem2fdsroot/$dem2fdsdir       $dem2fds    $bundledir/bin $dem2fdsout
SCP $fdshost $wind2fdsroot/$wind2fdsdir     $wind2fds   $bundledir/bin $wind2fdsout

mdversion=_$SMVVERSION.md5
MD5DIR=$bundledir/bin/MD5
MD5HASH $bundledir/bin $backgroundout > $MD5DIR/$backgroundout$mdversion
MD5HASH $bundledir/bin $smokeviewout  > $MD5DIR/$smokeviewout$mdversion
MD5HASH $bundledir/bin $smokediffout  > $MD5DIR/$smokediffout$mdversion
MD5HASH $bundledir/bin $smokezipout   > $MD5DIR/$smokezipout$mdversion
MD5HASH $bundledir/bin $dem2fdsout    > $MD5DIR/$dem2fdsout$mdversion
MD5HASH $bundledir/bin $wind2fdsout   > $MD5DIR/$wind2fdsout$mdversion

SCP $fdshost $smvscriptdir jp2conv.sh $bundledir/bin jp2conv.sh
CPDIR $texturedir $bundledir/bin/textures

# FDS 

SCP $fdshost $fdsroot/$fdsmpidir          $fdsmpi    $bundledir/bin $fdsmpiout
SCP $fdshost $fds2asciiroot/$fds2asciidir $fds2ascii $bundledir/bin $fds2asciiout

mdversion=_$FDSVERSION.md5
MD5HASH $bundledir/bin $fdsmpiout    > $MD5DIR/$fdsmpiout$mdversion
MD5HASH $bundledir/bin $fds2asciiout > $MD5DIR/$fds2asciiout$mdversion

if [ "$PLATFORM" == "LINUX64" ]; then
   ostype=LINUX
   ossize=intel64
   openmpifile=openmpi_${OPENMPI_VERSION}_linux_64.tar.gz
   MD5SUMMARY=$bundledir/bin/MD5/fds_${FDSVERSION}_linux_bundle.md5s
fi
if [ "$PLATFORM" == "OSX64" ]; then
   ostype=OSX
   ossize=intel64
   openmpifile=openmpi_${OPENMPI_VERSION}_osx_64.tar.gz
   MD5SUMMARY=$bundledir/bin/MD5/fds_${FDSVERSION}_osx_bundle.md5s
fi

echo ""
echo "--- copying configuration files ---"
echo ""
if [ "$OSXBUNDLE" == "yes" ]; then
  CP $fds_bundle FDS-SMV_OSX_Launcher.app.zip $bundledir/bin FDS-SMV_OSX_Launcher.app.zip
fi

CP $fds_bundle README.html   $bundledir/bin README.html

CP $smv_bundle smokeview.ini $bundledir/bin smokeview.ini

CP $smv_bundle volrender.ssf $bundledir/bin volrender.ssf

CP $smv_bundle objects.svo   $bundledir/bin objects.svo

CP $openmpidir $openmpifile  $bundledir/bin $openmpifile

echo ""
echo "--- copying documentation ---"
echo ""
CP2 $mandir FDS_Config_Management_Plan.pdf $bundledir/Documentation
CP2 $mandir FDS_Technical_Reference_Guide.pdf $bundledir/Documentation
CP2 $mandir FDS_User_Guide.pdf $bundledir/Documentation
CP2 $mandir FDS_Validation_Guide.pdf $bundledir/Documentation
CP2 $mandir FDS_Verification_Guide.pdf $bundledir/Documentation
CP2 $mandir SMV_User_Guide.pdf $bundledir/Documentation
CP2 $mandir SMV_Technical_Reference_Guide.pdf $bundledir/Documentation
CP2 $mandir SMV_Verification_Guide.pdf $bundledir/Documentation


if [ ! "$COMPFROM" == "" ]; then
  COMPLIBFROM=~/$COMPFROM
  if [ -d $COMPLIBFROM ]; then

    echo ""
    echo "--- copying compiler run time libraries ---"
    echo ""
    CPDIR $COMPLIBFROM $bundledir/bin/$COMPTO
  fi
fi
if [ ! "$MISCFROM" == "" ]; then
  MISCLIBFROM=~/$MISCFROM
  if [ -d $MISCLIBFROM ]; then

    echo ""
    echo "--- copying miscellaneous run time libraries ---"
    echo ""
    CPDIR $MISCLIBFROM $bundledir/bin/$MISCTO
  fi
fi

echo ""
echo "--- copying release notes ---"
echo ""
CP $fds_bundle FDS_Release_Notes.htm $bundledir/Documentation FDS_Release_Notes.html

CP $webpagesdir smv_readme.html $bundledir/Documentation SMV_Release_Notes.html


# CP2 $fds_bundle readme_examples.html $bundledir/Examples

export OUTDIR=$uploaddir/$bundledir/Examples
export QFDS=$copyfdscase
export RUNTFDS=$copyfdscase
export RUNCFAST=$copycfastcase

echo ""
echo "--- copying example files ---"
echo ""
cd $FDSExamplesDirectory
$fds_cases
$fds_benchmark_cases
cd $SMVExamplesDirectory
$wui_cases
$smv_cases
rm -rf $OUTDIR/Immersed_Boundary_Method

cd $curdir

echo ""
echo "--- building archive ---"
echo ""
rm -rf $uploaddir/$bundlebase.tar
rm -rf $uploaddir/$bundlebase.tar.gz
cd $uploaddir/$bundlebase
tar cf ../$bundlebase.tar --exclude='*.csv' .
echo Compressing archive
gzip    ../$bundlebase.tar
echo Creating installer
cd ..
$makeinstaller -o $ostype -i $bundlebase.tar.gz -d $INSTALLDIR $bundlebase.sh 
MD5HASH . $bundlebase.sh > $bundledir/bin/MD5/$bundlebase.sh.md5
MD5HASH . $bundlebase.sh > $uploaddir/$bundlebase.sh.md5
cat $bundledir/bin/MD5/*.md5 > $MD5SUMMARY

if [ -e $errlog ]; then
  numerrs=`cat $errlog | wc -l `
  if [ $numerrs -gt 0 ]; then
    echo ""
    echo "----------------------------------------------------------------"
    echo "---------------- bundle generation errors ----------------------"
    cat $errlog
    echo "----------------------------------------------------------------"
    echo "----------------------------------------------------------------"
    echo ""
  fi
  rm $errlog
fi
