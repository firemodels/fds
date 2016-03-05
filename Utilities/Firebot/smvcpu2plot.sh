#!/bin/bash
curdir=`pwd`

# by default assumes data comes from smokebot
# and output goes to the smokebot web page

indir=~smokebot/.smokebot
outdir=/var/www/html/smokebot
datadir=$HOME/.smokebot

function usage {
echo "Create a plot from smokebot timing data"
echo ""
echo "Options:"
echo "-F - force plot creation"
echo "-i - input directory [default: $indir]"
echo "-h - display this message"
echo "-o - output directory [default: $outdir]"
echo "-v - show options used, do not run"
exit
}

FORCE=
SHOW=
while getopts 'Fhi:o:v' OPTION
do
case $OPTION  in
  F)
   FORCE=-F
   ;;
  h)
   usage
   ;;
  i)
   indir="$OPTARG"
   ;;
  o)
   outdir="$OPTARG"
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
cpufrom=$indir/smv_times.csv

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

cpuplot=/tmp/smv_times.png.$$
old=$datadir/smv_times_trunc_old.csv
cputrunc=$datadir/smv_times_trunc.csv

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
cp $cpuplot $outdir/smv_times.png
rm $cpuplot
