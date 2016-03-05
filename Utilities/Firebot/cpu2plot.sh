#!/bin/bash
curdir=`pwd`
prefix=fds_

# by default assumes data comes from firebot
# and output goes to the firebot web page

indir=~firebot/.firebot
outdir=/var/www/html/firebot
datadir=$HOME/.firebot
smokebotdir=~smokebot/.smokebot

function usage {
echo "Create a plot from firebot or smokebot timing data"
echo ""
echo "Options:"
echo "-i - input directory [default: $indir]"
echo "-F - force plot creation"
echo "-h - display this message"
echo "-o - output directory [default: $outdir]"
echo "-p - prefix [default: $prefix]"
echo "-s - use smokebot for input [default: $smokebotdir]"
echo "-v - show options used, do not run"
exit
}

FORCE=
SHOW=
while getopts 'd:fFhi:o:p:sv' OPTION
do
case $OPTION  in
  d)
   datadir="$OPTARG"
   ;;
  h)
   usage
   ;;
  F)
   FORCE=-F
   ;;
  i)
   indir="$OPTARG"
   ;;
  o)
   outdir="$OPTARG"
   ;;
  p)
   prefix="$OPTARG"
   ;;
  s)
   indir=$smokebotdir
   datadir=$smokebotdir
   prefix=smv_
   ;;
  v)
   SHOW=1
   ;;
esac
done
shift $(($OPTIND-1))

cd $indir
indir=`pwd`

cd $curdir
if [ ! -d $outdir ]; then
  mkdir -p $outdir
fi
cd $outdir
outdir=`pwd`

cd $curdir

if [ "$SHOW" == "1" ]; then
   echo $0 $FORCE -i $indir -o $outdir   
   exit
fi
date=`date`
cpufrom=$indir/${prefix}times.csv

if [ ! -d $indir ]; then
  echo input directory $indir does not exist
  echo script aborted
  exit
fi
if [ ! -d $outdir ]; then
  echo output directory $outdir does not exist
  echo script aborted
  exit
fi
if [ ! -e $cpufrom ]; then
  echo cpu time file $cpufrom does not exist
  echo script aborted
  exit
fi
touch $outdir/test.$$
if [ ! -e $outdir/test.$$ ]; then
  echo unable to write to output directory $outdir
  echo script aborted
  exit
fi
rm $outdir/test.$$

cpuplot=/tmp/${prefix}times.png.$$
old=$datadir/${prefix}times_trunc_old.csv
cputrunc=$datadir/${prefix}times_trunc.csv

sort -n -k 1 -t , $cpufrom | tail -30 > $cputrunc
if [ "$FORCE" == "" ]; then
  if [ -e $old ]; then
    ndiff=`diff $old $cputrunc|wc -l`
    if [ "$ndiff" == "0" ]; then
      exit
    fi
  fi
fi

cp $cputrunc $old

cat << EOF | gnuplot
set terminal png size 900 600 giant
set xlabel "Days since Jan 1, 2016"
set ylabel "Benchmark Time (s)"
set output "$cpuplot"
set datafile separator ','
set style line 1 lt 1 lw 4 lc rgb "black"
set border ls 1
plot "$cputrunc" using 1:2 title "$date" with lines ls 1
EOF
cp $cpuplot $outdir/${prefix}times.png
rm $cpuplot
