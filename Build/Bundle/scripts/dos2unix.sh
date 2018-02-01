#!/bin/bash
infile=$1

tempfile=/tmp/$infile.$$
sed -e 's/\r//g' $infile > $tempfile
mv $tempfile $infile
