#!/bin/bash
version=$1
platform=$2
size=$3
HOST=$4
FDS_EDITION=$5
SVNROOT=$6

platformsize=${platform}_$size
BACKGROUNDDIR=$SVNROOT/SMV/Build/background/intel_${platform}_64
SMOKEVIEWDIR=$SVNROOT/SMV/Build/intel_$platformsize
SMOKEZIPDIR=$SVNROOT/SMV/Build/smokezip/intel_$platformsize
SMOKEDIFFDIR=$SVNROOT/SMV/Buiild/smokediff/intel_$platformsize
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
if [ "$platform" == "linux" ]; then
platform2=Linux
fi
if [ "$platform" == "osx" ]; then
platform2=OSX
fi

$UPDATER $platform2 $size $version $DIR.tar.gz $DIR.sh FDS/$FDS_EDITION release
