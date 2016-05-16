#!/bin/bash
revision=$1
REMOTESVNROOT=FDS-SMV
OSXHOST=$2
SVNROOT=~/$3

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


BACKGROUNDDIR=$REMOTESVNROOT/SMV/Build/background/intel_osx_64
SMVDIR=$REMOTESVNROOT/SMV/Build/smokeview/intel_osx_64
SMZDIR=$REMOTESVNROOT/SMV/Build/smokezip/intel_osx_64
DEM2FDSDIR=$REMOTESVNROOT/SMV/Build/dem2fds/intel_osx_64
SMDDIR=$REMOTESVNROOT/SMV/Build/smokediff/intel_osx_64
WINDDIR=$REMOTESVNROOT/SMV/Build/wind2fds/intel_osx_64
FORBUNDLE=$SVNROOT/SMV/for_bundle
OSXDIR=smv_test\_$revision\_osx64
UPDATER=$SVNROOT/Utilities/Scripts/make_updater.sh

cd $SVNROOT/SMV/uploads

rm -rf $OSXDIR
mkdir -p $OSXDIR
mkdir -p $OSXDIR/bin
mkdir -p $OSXDIR/Documentation

echo ""
echo "---- copying files ----"
echo ""
CP $FORBUNDLE objects.svo $OSXDIR/bin objects.svo
CP $FORBUNDLE smokeview.ini $OSXDIR/bin smokeview.ini
CPDIR $FORBUNDLE/textures $OSXDIR/bin/textures
cp $FORBUNDLE/*.po $OSXDIR/bin/.
CP $FORBUNDLE volrender.ssf $OSXDIR/bin volrender.ssf
SCP $OSXHOST $BACKGROUNDDIR background $OSXDIR/bin background
SCP $OSXHOST $SMVDIR smokeview_osx_test_64 $OSXDIR/bin smokeview
SCP $OSXHOST $DEM2FDSDIR dem2fds_osx_64 $OSXDIR/bin dem2fds
SCP $OSXHOST $SMDDIR smokediff_osx_64 $OSXDIR/bin smokediff
SCP $OSXHOST $WINDDIR wind2fds_osx_64 $OSXDIR/bin wind2fds
rm -f $OSXDIR.tar $OSXDIR.tar.gz
cd $OSXDIR
echo ""
echo "---- building installer ----"
echo ""
tar cvf ../$OSXDIR.tar .
cd ..
gzip $OSXDIR.tar
$UPDATER OSX $revision $OSXDIR.tar.gz $OSXDIR.sh FDS/FDS6 test
