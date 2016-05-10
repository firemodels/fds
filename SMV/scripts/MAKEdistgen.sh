#!/bin/bash
version=$1
platform=$2
HOST=$3
FDS_EDITION=$4
SVNROOT=$5

size=64

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

platformsize=${platform}_$size
BACKGROUNDDIR=$SVNROOT/SMV/Build/background/intel_${platform}_64
SMOKEVIEWDIR=$SVNROOT/SMV/Build/smokeview/intel_$platformsize
SMOKEZIPDIR=$SVNROOT/SMV/Build/smokezip/intel_$platformsize
DEM2FDSDIR=$SVNROOT/SMV/Build/dem2fds/intel_$platformsize
SMOKEDIFFDIR=$SVNROOT/SMV/Build/smokediff/intel_$platformsize
WINDDIR=$SVNROOT/SMV/Build/wind2fds/intel_$platformsize
FORBUNDLE=~/$SVNROOT/SMV/for_bundle
DIR=smv_${version}_$platform$size
UPDATER=~/$SVNROOT/Utilities/Scripts/make_updater.sh

cd ~/$SVNROOT/SMV/uploads

rm -rf $DIR
mkdir -p $DIR
mkdir -p $DIR/bin
mkdir -p $DIR/Documentation
cp ~/FDS-SMVwebpages/smv_readme.html $DIR/Documentation/release_notes.html

echo ""
echo "--- copying files ---"
echo ""
CPDIR $FORBUNDLE/textures $DIR/bin/textures
CP $FORBUNDLE objects.svo $DIR/bin objects.svo
CP $FORBUNDLE smokeview.ini $DIR/bin smokeview.ini
CP $FORBUNDLE volrender.ssf $DIR/bin volrender.ssf
SCP $HOST $BACKGROUNDDIR background $DIR/bin background
SCP $HOST $SMOKEDIFFDIR smokediff_$platformsize $DIR/bin smokediff
SCP $HOST $SMOKEVIEWDIR smokeview_$platformsize $DIR/bin smokeview
SCP $HOST $SMOKEZIPDIR smokezip_$platformsize $DIR/bin smokezip
SCP $HOST $DEM2FDSDIR dem2fds_$platformsize $DIR/bin dem2fds
SCP $HOST $WINDDIR wind2fds_$platformsize $DIR/bin wind2fds

echo ""
echo "--- building installer ---"
echo ""
rm -f $DIR.tar $DIR.tar.gz
cd $DIR
tar cvf ../$DIR.tar .
cd ..
gzip $DIR.tar

platform2=$platform
if [ "$platform" == "linux" ]; then
platform2=Linux
fi
if [ "$platform" == "osx" ]; then
platform2=OSX
fi

$UPDATER $platform2 $size $version $DIR.tar.gz $DIR.sh FDS/$FDS_EDITION release
