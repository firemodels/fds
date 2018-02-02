#!/bin/bash
infile=$1

tempfile=/tmp/$infile.$$
tr -d '\r' <  $infile > $tempfile
mv $tempfile $infile
