#!/bin/bash
base=$1
outdir=$2
underscore=_
#jpeg2yuv -f 25 -I p -j $base$underscore%04d.jpg | mpeg2enc -o $outdir/$base.m1v
png2yuv -f 25 -I p -j $base$underscore%04d.png | mpeg2enc -o $outdir/$base.m1v
echo The movie file $outdir/$base.m1v has been created.
