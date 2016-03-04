#!/bin/bash
curdir=`pwd`

indir=~firebot/.firebot
outdir=/var/www/html/firebot

if [ ! -d $indir ]; then
  mkdir -p $indir
fi
cd $indir
firebotdir=~/.firebot
if [ ! -d $firebotdir ]; then
  mkdir -p $firebotdir
fi
cd $curdir

smokebotdir=~smokebot/.smokebot

function usage {
echo "Create a plot from fds timing data"
echo ""
echo "Options:"
echo "-i - input directory [default: $indir]"
echo "-F - force plot creation"
echo "-h - display this message"
echo "-o - output directory [default: $outdir]"
echo "-s - use smokebot for input [default: $smokebotdir]"
echo "-v - show options used, do not run"
exit
}

FORCE=
SHOW=
while getopts 'fFhi:o:sv' OPTION
do
case $OPTION  in
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
  s)
   indir=$smokebotdir
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
cd $outdir
outdir=`pwd`
cd $curdir

if [ "$SHOW" == "1" ]; then
   echo $0 $FORCE -i $indir -o $outdir   
   exit
fi
date=`date`
cpufrom=$indir/fds_times.csv

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
  echo unable to write to outdir $outdir
  echo script aborted
  exit
fi
rm $outdir/test.$$

cpuplot=/tmp/fds_times.png.$$
old=$firebotdir/old
cputo=$firebotdir/fds_times.csv
cputrunc=$firebotdir/fds_times_trunc.csv

tail -30 $cpufrom > $cputrunc
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
cp $cpuplot $outdir/fds_times.png
#rm $cpuplot
