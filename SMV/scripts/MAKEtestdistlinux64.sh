#!/bin/bash
revision=$1
SVNROOT=~/$2

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


BACKGROUNDDIR=$SVNROOT/SMV/Build/background/intel_linux_64
SMVDIR=$SVNROOT/SMV/Build/smokeview/intel_linux_64
SMZDIR=$SVNROOT/SMV/Build/smokezip/intel_linux_64
DEM2FDSDIR=$SVNROOT/SMV/Build/dem2fds/intel_linux_64
SMDDIR=$SVNROOT/SMV/Build/smokediff/intel_linux_64
WINDDIR=$SVNROOT/SMV/Build/wind2fds/intel_linux_64
FORBUNDLE=$SVNROOT/SMV/for_bundle
LINUXDIR=smv_test\_$revision\_linux64
UPDATER=$SVNROOT/Utilities/Scripts/make_updater.sh

cd $SVNROOT/SMV/uploads

echo ""
echo "--- copying files ---"
echo ""
rm -rf $LINUXDIR
mkdir -p $LINUXDIR
mkdir -p $LINUXDIR/bin
mkdir -p $LINUXDIR/Documentation
CPDIR $FORBUNDLE/textures $LINUXDIR/bin/textures
CP $FORBUNDLE objects.svo $LINUXDIR/bin objects.svo
CP $FORBUNDLE smokeview.ini $LINUXDIR/bin smokeview.ini
CP $FORBUNDLE volrender.ssf $LINUXDIR/bin volrender.ssf
cp $FORBUNDLE/*.po $LINUXDIR/bin/.
CP $BACKGROUNDDIR background $LINUXDIR/bin background
CP $SMVDIR smokeview_linux_test_64 $LINUXDIR/bin smokeview
CP $DEM2FDSDIR dem2fds_linux_64 $LINUXDIR/bin dem2fds
CP $SMDDIR smokediff_linux_64 $LINUXDIR/bin smokediff
CP $WINDDIR wind2fds_linux_64 $LINUXDIR/bin wind2fds
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
cd $LINUXDIR
echo ""
echo "--- building installer ---"
echo ""
tar cvf ../$LINUXDIR.tar .
cd ..
gzip $LINUXDIR.tar
$UPDATER Linux $revision $LINUXDIR.tar.gz $LINUXDIR.sh FDS/FDS6 test
