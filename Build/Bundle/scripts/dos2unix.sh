#!/bin/bash
infile=$1

tempfile=/tmp/$infile.$$
tr -d '\r' <  $HOME/$infile > $tempfile
mv $tempfile $HOME/$infile
