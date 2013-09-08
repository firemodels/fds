#!/bin/csh -f
set version=$1
set platform=$2
set size=$3
set SVNROOT=FDS-SMV
set HOST=$4

set platformsize=${platform}_$size
set BACKGROUNDDIR=$SVNROOT/Utilities/background/intel_${platform}_32
set SMOKEVIEWDIR=$SVNROOT/SMV/Build/intel_$platformsize
set SMOKEZIPDIR=$SVNROOT/Utilities/smokezip/intel_$platformsize
set SMOKEDIFFDIR=$SVNROOT/Utilities/smokediff/intel_$platformsize
set FORBUNDLE=~/$SVNROOT/SMV/for_bundle
set DIR=smv_${version}_$platform$size

cd $FORBUNDLE/uploads

rm -rf $DIR
mkdir -p $DIR
mkdir -p $DIR/Documentation
cp $FORBUNDLE/readme.html $DIR/Documentation/release_notes.html

cp $FORBUNDLE/objects.svo $DIR/.
scp $HOST\:$BACKGROUNDDIR/background $DIR/.
scp $HOST\:$SMOKEDIFFDIR/smokediff_$platformsize $DIR/.
scp $HOST\:$SMOKEVIEWDIR/smokeview_$platformsize $DIR/.
scp $HOST\:$SMOKEZIPDIR/smokezip_$platformsize $DIR/.
rm -f $DIR.tar $DIR.tar.gz
tar cvf $DIR.tar $DIR/.
gzip $DIR.tar
