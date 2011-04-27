#!/bin/bash
base=$1
out=$2
#jpeg2yuv -f 25 -I p -j $base%04d.jpg | mpeg2enc -o $out.m1v
png2yuv -f 25 -I p -j $base%04d.png | mpeg2enc -o $out.m1v
