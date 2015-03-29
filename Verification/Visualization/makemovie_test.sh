#!/bin/bash
rm -f movie_test*.png 2>/dev/null
echo making images
smokeview -runscript movie_test > /dev/null 2>/dev/null
movietype=.avi
if [ "`uname`" == "Darwin" ]; then
  movietype=.mp4
fi
echo making movie
ffmpeg -y -i movie_test_%04d.png movie_test$movietype