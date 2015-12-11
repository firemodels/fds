#!/bin/csh -f
set version=$1
set platform=$2
set size=$3
set HOST=$4
set FDS_EDITION=$5
set SVNROOT=$6

set platformsize=${platform}_$size
set BACKGROUNDDIR=$SVNROOT/Utilities/background/intel_${platform}_32
set SMOKEVIEWDIR=$SVNROOT/SMV/Build/intel_$platformsize
set SMOKEZIPDIR=$SVNROOT/Utilities/smokezip/intel_$platformsize
set SMOKEDIFFDIR=$SVNROOT/Utilities/smokediff/intel_$platformsize
set WINDDIR=$SVNROOT/Utilities/wind2fds/intel_$platformsize
set FORBUNDLE=~/$SVNROOT/SMV/for_bundle
set DIR=smv_${version}_$platform$size
set UPDATER=~/$SVNROOT/Utilities/Scripts/make_updater.sh

cd ~/$SVNROOT/SMV/uploads

rm -rf $DIR
mkdir -p $DIR
mkdir -p $DIR/bin
mkdir -p $DIR/Documentation
cp $FORBUNDLE/readme.html $DIR/Documentation/release_notes.html

cp -r $FORBUNDLE/textures $DIR/bin/.
cp $FORBUNDLE/objects.svo $DIR/bin/.
cp $FORBUNDLE/smokeview.ini $DIR/bin/.
cp $FORBUNDLE/volrender.ssf $DIR/bin/.
scp $HOST\:$BACKGROUNDDIR/background $DIR/bin/.
scp $HOST\:$SMOKEDIFFDIR/smokediff_$platformsize $DIR/bin/smokediff
scp $HOST\:$SMOKEVIEWDIR/smokeview_$platformsize $DIR/bin/smokeview
scp $HOST\:$SMOKEZIPDIR/smokezip_$platformsize $DIR/bin/smokezip
scp $HOST\:$WINDDIR/wind2fds_$platformsize $DIR/bin/wind2fds
rm -f $DIR.tar $DIR.tar.gz
cd $DIR
tar cvf ../$DIR.tar .
cd ..
gzip $DIR.tar

set platform2=$platform
if ( "$platform" == "linux" ) then
set platform2=Linux
endif
if ( "$platform" == "osx" ) then
set platform2=OSX
endif

$UPDATER $platform2 $size $version $DIR.tar.gz $DIR.sh FDS/$FDS_EDITION release
