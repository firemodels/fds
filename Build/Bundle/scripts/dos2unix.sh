#!/bin/bash
dir=$1
infile=$2

tempfile=/tmp/$infile.$$
cd $dir
tr -d '\r' <  $infile > $tempfile
mv $tempfile $infile
