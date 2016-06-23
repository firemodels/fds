#!/bin/bash

curdir=`pwd`
cd ../..
reporoot=`pwd`
dir1=$reporoot/Utilities/Scripts
dir2=$reporoot/FDS/Build/Scripts
dir3=$reporoot/SMV/scripts
dir4=$reporoot/Verification/scripts
dir5=$reporoot/Utilities/Firebot
h1=$dir1/../build_bundle.html
h2=$dir2/../build_fds.html
h3=$dir3/../build_smokeview.html
h4=$dir3/../build_guides.html

cd $curdir

#FILES="$dir3/*.sh $dir3/*.bat"
FILES="$dir1/*.sh $dir1/*.bat"
for f in $FILES
do
file="${f##*/}"
file2=$file
file="${file%.*}"
numhits=`grep -i $file $dir1/*.sh $dir1/*.bat $dir2/*.bat $dir2/*.sh $dir3/*.sh $dir3/*.bat $dir4/*.sh $dir4/*.bat $dir5/*.bat $dir5/*.sh $h1 $h2 $h3 $h4 | wc -l`
if [ "$numhits" == "0" ]; then
echo file=$file2 hits=$numhits
fi
done
