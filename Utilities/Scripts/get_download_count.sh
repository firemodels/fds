#!/bin/bash
FDSURL=https://api.github.com/repos/firemodels/fds/releases
SMVURL=https://api.github.com/repos/firemodels/smv/releases
BUNDLE=FDS_6.6.0-SMV_6.6.0
SMV=SMV6.6.4
LINUX=_linux64.sh
OSX=_osx64.sh
WIN=_win64.exe

get_downloads()
{
FILE=$1
URL=$2
wget -q -O-  $URL | grep -E '\"name\":|download_count'   | grep --no-group-separator -B 1 download_count  | grep $FILE -A 1 | grep download_count | head -1 | awk '{print $2}' | awk -F"," '{print $1}'
}
echo bundle downloads
echo FILE=$BUNDLE$LINUX downloads=`get_downloads $BUNDLE$LINUX $FDSURL`
echo FILE=$BUNDLE$OSX   downloads=`get_downloads $BUNDLE$OSX $FDSURL`
echo FILE=$BUNDLE$WIN   downloads=`get_downloads $BUNDLE$WIN $FDSURL`

echo ""
echo smokeview downloads
echo FILE=$SMV$LINUX downloads=`get_downloads $SMV$LINUX $SMVURL`
echo FILE=$SMV$OSX   downloads=`get_downloads $SMV$OSX   $SMVURL`
echo FILE=$SMV$WIN   downloads=`get_downloads $SMV$WIN   $SMVURL`
