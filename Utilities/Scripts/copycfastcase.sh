#!/bin/bash
while getopts 'd:p:' OPTION
do
case $OPTION  in
  d)
   dir="$OPTARG"
   ;;
  p)
   dummy="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

in=$1
infile=${in%.*}

cd $dir
if ! [ -d $OUTDIR/$dir ]; then
  mkdir $OUTDIR/$dir
fi

cp $infile.in $OUTDIR/$dir/.
if [ -e $infile.ini ]; then
  cp $infile.ini $OUTDIR/$dir/.
fi

if [ -e $infile.ssf ]; then
  cp $infile.ssf $OUTDIR/$dir/.
fi
