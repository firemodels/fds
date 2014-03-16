#!/bin/bash
if [ $# -lt 1 ] ; then
  echo "Usage: makemovie.sh [-i input_directory] [-d output_directory] [-m movie_name] input_base"
  echo ""
  echo "This script generates a movie from a sequence of "
  echo "png image files.  Each image file has the form basexxxx.png"
  echo "where xxxx is a frame number."
  echo ""
  echo "-i dir - directory where movie frames are located (default: .)"
  echo "-o dir - directory where movie will be placed (default: .)"
  echo "-m movie name - name of movie generated (default: input_base.m1v)"
  echo ""
  exit
fi

indir=.
outdir=.
moviename=
while getopts 'i:o:m:' OPTION
do
case $OPTION in
  i)
  indir="$OPTARG"
  ;;
  o)
  outdir="$OPTARG"
  ;;
  m)
  moviename="$OPTARG"
  ;;
esac
done
shift $((OPTIND-1))

CURDIR=`pwd`

cd $outdir
outdir=`pwd`

cd $CURDIR
cd $indir

base=$1
underscore=_
if [ "$moviename" == "" ] ; then
  moviename=$base.m1v
else
  moviename=$moviename.m1v
fi

echoerr() { echo "$@" 1>&2; }
echoerr Creating the movie file $outdir/$moviename
png2yuv -f 25 -I p -j $base$underscore%04d.png | mpeg2enc -o $outdir/$moviename
#ffmpeg -f image2 -i $base$underscore%04d.png $outdir/$moviename
echoerr The movie file $outdir/$moviename has been created.
