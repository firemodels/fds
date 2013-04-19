#!/bin/bash
DIR=$1
case=$2
cd $DIR
rm -f $case*.s3d
rm -f $case*.sf
rm -f $case*.bf
rm -f $case*.iso
rm -f $case*.sz
rm -f $case*.csv
rm -f $case*.q
rm -f $case*.restart
