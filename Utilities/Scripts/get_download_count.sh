#!/bin/bash
URL=https://api.github.com/repos/firemodels/fds/releases
FILE1=FDS_6.6.0-SMV_6.6.0_win64.exe
FILE2=FDS_6.6.0-SMV_6.6.0_linux64.sh
FILE3=FDS_6.6.0-SMV_6.6.0_osx64.sh

get_downloads()
{
FILE=$1
wget -q -O-  $URL | grep -E '\"name\":|download_count'   | grep --no-group-separator -B 1 download_count  | grep $FILE -A 1 | grep download_count | head -1 | awk '{print $2}' | awk -F"," '{print $1}'
}

echo FILE=$FILE1 downloads=`get_downloads $FILE1`
echo FILE=$FILE2 downloads=`get_downloads $FILE2`
echo FILE=$FILE3 downloads=`get_downloads $FILE3`
