#!/bin/bash
base=$1
outdir=$2
underscore=_
echoerr() { echo "$@" 1>&2; }
echoerr Creating the movie file $outdir/$base.m1v
#jpeg2yuv -f 25 -I p -j $base$underscore%04d.jpg | mpeg2enc -o $outdir/$base.m1v
png2yuv -f 25 -I p -j $base$underscore%04d.png | mpeg2enc -v 0 -o $outdir/$base.m1v
echoerr The movie file $outdir/$base.m1v has been created.
